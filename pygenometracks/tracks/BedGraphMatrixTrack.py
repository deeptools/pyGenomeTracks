from . BedGraphTrack import BedGraphTrack
from . GenomeTrack import GenomeTrack
import numpy as np
import matplotlib.pyplot as plt


class BedGraphMatrixTrack(BedGraphTrack):
    SUPPORTED_ENDINGS = ['.bm', '.bm.gz', '.bedgraphmatrix', '.bm.bgz']
    TRACK_TYPE = 'bedgraph_matrix'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# a bedgraph matrix file is like a bedgraph, except that per bin there
# are more than one value separated by tab: E.g.
# This file type is produced by HiCExplorer tool hicFindTads and contains
# the TAD-separation score at different window sizes
# chrX	18279	40131	0.399113	0.364118	0.320857	0.274307
# chrX	40132	54262	0.479340	0.425471	0.366541	0.324736
#min_value = 0.10
#max_value = 0.70
# if type is set as lines, then the TAD score lines are drawn instead
# of the matrix otherwise a heatmap is plotted
type = lines
# If the type is not lines, you can choose to keep the matrix as not rasterized
# (only used if you use pdf or svg output format) by using:
# rasterize = false
# pos_score_in_bin means 'position of score with respect to bin start and end'
# if the lines option is used, the y values can be put at the
# center of the bin (default) or they can be plot as 'block',
# which mean to plot the values as a line between the start and end of bin
pos_score_in_bin = center
show_data_range = true

# only when type lines is used. Adds horizontal lines
plot_horizontal_lines = false
file_type = {}
    """.format(TRACK_TYPE)
    DEFAULTS_PROPERTIES = {'max_value': None,
                           'min_value': None,
                           'type': 'matrix',
                           'pos_score_in_bin': 'center',
                           'show_data_range': True,
                           'plot_horizontal_lines': False,
                           'orientation': None,
                           'rasterize': True}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'type': ['matrix', 'lines'],
                           'pos_score_in_bin': ['center', 'block'],
                           'orientation': [None, 'inverted']}
    BOOLEAN_PROPERTIES = ['show_data_range', 'plot_horizontal_lines',
                          'rasterize']
    STRING_PROPERTIES = ['file', 'file_type', 'overlay_previous',
                         'type', 'pos_score_in_bin', 'orientation',
                         'title']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}
    # The color cannot be set for the moment

    def set_properties_defaults(self):
        GenomeTrack.set_properties_defaults(self)

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plots a bedgraph matrix file, that instead of having
        a single value per bin, it has several values.
        """
        values_list, start_pos = self.get_scores(chrom_region, start_region, end_region)
        if start_pos == []:
            return
        matrix_rows = []
        for values in values_list:
            values = list(map(float, values))
            matrix_rows.append(values)

        matrix = np.vstack(matrix_rows).T
        if self.properties['orientation'] == 'inverted':
            matrix = np.flipud(matrix)

        if self.properties['type'] == 'lines':
            if self.properties['pos_score_in_bin'] == 'block':
                # convert [(0, 10), (10, 20), (20, 30)] into [0, 10, 10, 20, 20, 30]
                x_values = sum(start_pos, tuple())
            else:
                x_values = [x[0] + (x[1] - x[0]) / 2 for x in start_pos]

            for row in matrix:
                if self.properties['pos_score_in_bin'] == 'block':
                    # convert [1, 2, 3 ...] in [1, 1, 2, 2, 3, 3 ...]
                    row = np.repeat(row, 2)
                ax.plot(x_values, row, color='grey', linewidth=0.5)

            if self.properties['pos_score_in_bin'] == 'block':
                mean_values = np.repeat(matrix.mean(axis=0), 2)
            else:
                mean_values = matrix.mean(axis=0)
            ax.plot(x_values, mean_values, linestyle="--", marker="|")
            ymax = self.properties['max_value']
            ymin = self.properties['min_value']
            ax.set_ylim(ymin, ymax)

            if self.properties['plot_horizontal_lines']:
                ax.grid(True)
                ax.grid(True, axis='y')
                ax.axhline(y=0, color='black', linewidth=1)
                ax.tick_params(axis='y', which='minor', left='on')

        else:
            start_pos = [x[0] for x in start_pos]

            x, y = np.meshgrid(start_pos, np.arange(matrix.shape[0]))
            shading = 'gouraud'
            vmax = self.properties['max_value']
            vmin = self.properties['min_value']
            self.img = ax.pcolormesh(x, y, matrix, vmin=vmin, vmax=vmax, shading=shading)
            if self.properties['rasterize']:
                self.img.set_rasterized(True)

    def plot_y_axis(self, ax, plot_axis):
        if self.properties['type'] == 'lines':
            super(BedGraphMatrixTrack, self).plot_y_axis(ax, plot_axis)
        else:
            try:
                cobar = plt.colorbar(self.img, ax=ax, fraction=0.95)
            except AttributeError:
                return

            cobar.solids.set_edgecolor("face")
            cobar.ax.tick_params(labelsize='smaller')
            cobar.ax.yaxis.set_ticks_position('left')
            # adjust the labels of the colorbar
            ticks = cobar.ax.get_yticks()
            labels = cobar.ax.set_yticklabels(ticks.astype('float32'))
            (vmin, vmax) = cobar.mappable.get_clim()
            for idx in np.where(ticks == vmin)[0]:
                # if the label is at the start of the colobar
                # move it above avoid being cut or overlapping with other track
                labels[idx].set_verticalalignment('bottom')
            for idx in np.where(ticks == vmax)[0]:
                # if the label is at the end of the colobar
                # move it a bit inside to avoid overlapping
                # with other labels
                labels[idx].set_verticalalignment('top')
