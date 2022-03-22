from . BedGraphTrack import BedGraphTrack
from . GenomeTrack import GenomeTrack
import numpy as np
from matplotlib import cm

DEFAULT_BEDGRAPHMATRIX_COLORMAP = 'viridis'
DEFAULT_BEDGRAPHMATRIX_INDIVIDUAL = 'grey'
DEFAULT_BEDGRAPHMATRIX_SUMMARY = '#1f77b4'  # mpl.rcParams["axes.prop_cycle"].by_key()["color"][0]


class BedGraphMatrixTrack(BedGraphTrack):
    SUPPORTED_ENDINGS = ['.bm', '.bm.gz', '.bedgraphmatrix', '.bm.bgz']
    TRACK_TYPE = 'bedgraph_matrix'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
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
# by default individual lines are colored in grey:
individual_color = grey
# by default the summary line is colored in blue:
summary_color = #1f77b4
# If the type is not lines, you can choose to keep the matrix as not rasterized
# (only used if you use pdf or svg output format) by using:
# rasterize = false
# The different options for color maps can be found here:
# https://matplotlib.org/users/colormaps.html
# the default color map is viridis
# If you want your own colormap you can put the values of the color you want
# For example, colormap = ['blue', 'yellow', 'red']
# or colormap = ['white', (1, 0.88, .66), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
#colormap = Reds
# pos_score_in_bin means 'position of score with respect to bin start and end'
# if the lines option is used, the y values can be put at the
# center of the bin (default) or they can be plot as 'block',
# which mean to plot the values as a line between the start and end of bin
pos_score_in_bin = center
show_data_range = true

# only when type lines is used. Adds horizontal lines
plot_horizontal_lines = false
file_type = {TRACK_TYPE}
    """
    DEFAULTS_PROPERTIES = {'max_value': None,
                           'min_value': None,
                           # In next 1.0 change matrix to lines
                           'type': 'matrix',
                           'pos_score_in_bin': 'center',
                           'show_data_range': True,
                           'plot_horizontal_lines': False,
                           'orientation': None,
                           'rasterize': True,
                           'region': None,  # Cannot be set manually but is set by tracksClass
                           'colormap': DEFAULT_BEDGRAPHMATRIX_COLORMAP,
                           'individual_color': DEFAULT_BEDGRAPHMATRIX_INDIVIDUAL,
                           'summary_color': DEFAULT_BEDGRAPHMATRIX_SUMMARY}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None},
                             'type': {'line': 'lines'}}
    POSSIBLE_PROPERTIES = {'type': ['matrix', 'lines'],
                           'pos_score_in_bin': ['center', 'block'],
                           'orientation': [None, 'inverted']}
    BOOLEAN_PROPERTIES = ['show_data_range', 'plot_horizontal_lines',
                          'rasterize']
    STRING_PROPERTIES = ['file', 'file_type', 'overlay_previous',
                         'type', 'pos_score_in_bin', 'orientation',
                         'title', 'colormap', 'individual_color',
                         'summary_color']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}
    # The colormap can only be a colormap
    # The individual_color and summary_color can only be a color

    def __init__(self, properties_dict):
        GenomeTrack.__init__(self, properties_dict)

        self.load_file()

    def set_properties_defaults(self):
        # To remove in next 1.0
        if 'type' not in self.properties:
            self.log.warning("Deprecated Warning: The section "
                             f"{self.properties['section_name']} did"
                             " not specify the type. For the moment"
                             " the default type is matrix but in the"
                             " next version it will be lines.\n")
        # End to remove
        GenomeTrack.set_properties_defaults(self)
        if self.properties['type'] == 'matrix':
            self.process_color('colormap', colormap_possible=True,
                               colormap_only=True,
                               default_value_is_colormap=True)
            self.cmap = cm.get_cmap(self.properties['colormap'])
        else:
            for param in ['individual_color', 'summary_color']:
                self.process_color(param, colormap_possible=False)

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plots a bedgraph matrix file, that instead of having
        a single value per bin, it has several values.
        """
        values_list, start_pos = self.get_scores(chrom_region, start_region, end_region)
        if start_pos == []:
            self.adjust_ylim(ax)
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
                ax.plot(x_values, row, color=self.properties['individual_color'],
                        linewidth=0.5)

            if self.properties['pos_score_in_bin'] == 'block':
                mean_values = np.repeat(matrix.mean(axis=0), 2)
            else:
                mean_values = matrix.mean(axis=0)
            ax.plot(x_values, mean_values, linestyle="--", marker="|",
                    color=self.properties['summary_color'])
            self.adjust_ylim(ax)

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
            self.last_img_plotted = ax.pcolormesh(x, y, matrix, vmin=vmin, vmax=vmax,
                                                  shading=shading, cmap=self.cmap)
            if self.properties['rasterize']:
                self.last_img_plotted.set_rasterized(True)

    def plot_y_axis(self, ax, plot_axis):
        if self.properties['type'] == 'lines':
            GenomeTrack.plot_y_axis(self, ax, plot_axis)
        else:
            GenomeTrack.plot_custom_cobar(self, ax)

    def __del__(self):
        if self.tbx is not None:
            self.tbx.close()
