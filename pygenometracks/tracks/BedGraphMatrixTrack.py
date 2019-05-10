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
# pos score in bin means 'position of score with respect to bin start and end'
# if the lines option is used, the y values can be put at the
# center of the bin (default) or they can be plot as 'block',
# which mean to plot the values as a line between the start and end of bin
pos score in bin = center
show data range = yes

# only when type lines is used. Adds horizontal lines
plot horizontal lines = no
file_type = {}
    """.format(TRACK_TYPE)

    def set_properties_defaults(self):
        if 'max_value' not in self.properties or self.properties['max_value'] == 'auto':
            self.properties['max_value'] = None

        if 'min_value' not in self.properties or self.properties['min_value'] == 'auto':
            self.properties['min_value'] = None

        if 'type' not in self.properties:
            self.properties['type'] = 'matrix'

        if 'pos score in bin' not in self.properties:
            self.properties['pos score in bin'] = 'center'

        if 'show data range' not in self.properties:
            self.properties['show data range'] = 'yes'

        if 'plot horizontal lines' not in self.properties:
            self.properties['plot horizontal lines'] = 'no'

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plots a bedgraph matrix file, that instead of having
        a single value per bin, it has several values.
        """
        try:
            values_list, start_pos = self.get_scores(chrom_region, start_region, end_region)
        except TypeError:
            return
        matrix_rows = []
        for values in values_list:
            values = list(map(float, values))
            matrix_rows.append(values)

        matrix = np.vstack(matrix_rows).T
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            matrix = np.flipud(matrix)

        if self.properties['type'] == 'lines':
            if self.properties['pos score in bin'] == 'block':
                # convert [(0, 10), (10, 20), (20, 30)] into [0, 10, 10, 20, 20, 30]
                x_values = sum(start_pos, tuple())
            else:
                x_values = [x[0] + (x[1] - x[0]) / 2 for x in start_pos]

            for row in matrix:
                if self.properties['pos score in bin'] == 'block':
                    # convert [1, 2, 3 ...] in [1, 1, 2, 2, 3, 3 ...]
                    row = np.repeat(row, 2)
                ax.plot(x_values, row, color='grey', linewidth=0.5)

            if self.properties['pos score in bin'] == 'block':
                mean_values = np.repeat(matrix.mean(axis=0), 2)
            else:
                mean_values = matrix.mean(axis=0)
            ax.plot(x_values, mean_values, linestyle="--", marker="|")
            ymax = self.properties['max_value']
            ymin = self.properties['min_value']
            ax.set_ylim(ymin, ymax)

            if self.properties['plot horizontal lines'] == 'yes':
                ax.grid(True)
                ax.grid(True, which='y')
                ax.axhline(y=0, color='black', linewidth=1)
                ax.tick_params(axis='y', which='minor', left='on')

        else:
            start_pos = [x[0] for x in start_pos]

            x, y = np.meshgrid(start_pos, np.arange(matrix.shape[0]))
            shading = 'gouraud'
            vmax = self.properties['max_value']
            vmin = self.properties['min_value']
            self.img = ax.pcolormesh(x, y, matrix, vmin=vmin, vmax=vmax, shading=shading)
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
            labels = cobar.ax.get_yticklabels()
            ticks = cobar.ax.get_yticks()
            if ticks[0] == 0:
                # if the label is at the start of the colobar
                # move it above avoid being cut or overlapping with other track
                labels[0].set_verticalalignment('bottom')
            if ticks[-1] == 1:
                # if the label is at the end of the colobar
                # move it a bit inside to avoid overlapping
                # with other labels
                labels[-1].set_verticalalignment('top')
            cobar.ax.set_yticklabels(labels)
