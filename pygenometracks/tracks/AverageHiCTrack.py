from hicmatrix import HiCMatrix
import hicmatrix.utilities
import scipy.sparse
from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix, load_npz

from . GenomeTrack import GenomeTrack

DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


class AverageHiCTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.npz']
    TRACK_TYPE = 'average_hic_matrix'
    OPTIONS_TXT = """
title =
# The different options for color maps can be found here: https://matplotlib.org/users/colormaps.html
# the default color map is RdYlBu_r (_r) stands for reverse
#colormap = RdYlBu_r
# depth is the maximum distance that should be plotted.
depth = 100000
# height of track (in cm) can be given. Otherwise, the height is computed such that the proportions of the
# hic matrix are kept (e.g. the image does not appear shrink or extended)
# height = 10
# min_value and max_value refer to the contacts in the matrix.
#min_value = 2.8
#max_value = 3.0
# the matrix can be transformed using the log1 (or log, but zeros could be problematic)
transform = log1p
# show masked bins plots as white lines
# those bins that were not used during the correction
# the default is to extend neighboring bins to
# obtain an aesthetically pleasant output
show_masked_bins = no
# if the track wants to be plotted upside-down:
# orientation = inverted
# optional if the values in the matrix need to be scaled the
# following parameter can be used. This is useful to plot multiple hic-matrices on the same scale
# scale factor = 1
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwargs):
        super(AverageHiCTrack, self).__init__(*args, **kwargs)
        log.debug("correct call of init!")
        self.matrix = load_npz(self.properties['file'])

        self.plot_inverted = False
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            self.plot_inverted = True

        self.norm = None

        if 'colormap' not in self.properties:
            self.properties['colormap'] = DEFAULT_MATRIX_COLORMAP
        self.cmap = cm.get_cmap(self.properties['colormap'])
        self.cmap.set_bad('white')

        self.cmap.set_bad('black')

    def plot(self, ax, chrom_region, region_start, region_end):

        log.debug('type; {}'.format(type(self.matrix)))
        # select only relevant matrix part
        matrix = self.matrix
        
        matrix = np.asarray(matrix.todense().astype(float))
        if 'scale factor' in self.properties:
            matrix = matrix * self.properties['scale factor']

        if 'transform' in self.properties:
            if self.properties['transform'] == 'log1p':
                matrix += 1
                self.norm = colors.LogNorm()

            elif self.properties['transform'] == '-log':
                mask = matrix == 0
                matrix[mask] = matrix[mask is False].min()
                matrix = -1 * np.log(matrix)

            elif self.properties['transform'] == 'log':
                mask = matrix == 0
                matrix[mask] = matrix[mask is False].min()
                matrix = np.log(matrix)

        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            vmax = self.properties['max_value']

        else:
            # try to use a 'aesthetically pleasant' max value
            vmax = np.percentile(matrix.diagonal(1), 80)

        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            vmin = self.properties['min_value']
       
        self.img = self.pcolormesh_45deg(ax, matrix, start_pos, vmax=vmax, vmin=vmin)
        self.img.set_rasterized(True)
        if self.plot_inverted:
            ax.set_ylim(depth, 0)
        else:
            ax.set_ylim(0, depth)

    # def plot_y_axis(self, cbar_ax, plot_ax):

    #     if 'transform' in self.properties and \
    #             self.properties['transform'] in ['log', 'log1p']:
    #         # get a useful log scale
    #         # that looks like [1, 2, 5, 10, 20, 50, 100, ... etc]

    #         # The following code is problematic with some versions of matplotlib.
    #         # Should be uncommented once the problem is clarified
    #         from matplotlib.ticker import LogFormatter
    #         formatter = LogFormatter(10, labelOnlyBase=False)
    #         aa = np.array([1, 2, 5])
    #         tick_values = np.concatenate([aa * 10 ** x for x in range(10)])
    #         cobar = plt.colorbar(self.img, ticks=tick_values, format=formatter, ax=cbar_ax, fraction=0.95)
    #     else:
    #         cobar = plt.colorbar(self.img, ax=cbar_ax, fraction=0.95)

    #     cobar.solids.set_edgecolor("face")
    #     cobar.ax.tick_params(labelsize='smaller')
    #     cobar.ax.yaxis.set_ticks_position('left')

    #     # adjust the labels of the colorbar
    #     labels = cobar.ax.get_yticklabels()
    #     ticks = cobar.ax.get_yticks()

    #     if ticks[0] == 0:
    #         # if the label is at the start of the colobar
    #         # move it above avoid being cut or overlapping with other track
    #         labels[0].set_verticalalignment('bottom')
    #     if ticks[-1] == 1:
    #         # if the label is at the end of the colobar
    #         # move it a bit inside to avoid overlapping
    #         # with other labels
    #         labels[-1].set_verticalalignment('top')

    #     cobar.ax.set_yticklabels(labels)

    def pcolormesh_45deg(self, ax, matrix_c, start_pos_vector, vmin=None, vmax=None):
        """
        Turns the matrix 45 degrees and adjusts the
        bins to match the actual start end positions.
        """
        import itertools
        # code for rotating the image 45 degrees
        n = matrix_c.shape[0]
        # create rotation/scaling matrix
        t = np.array([[1, 0.5], [-1, 0.5]])
        # create coordinate matrix and transform it
        matrix_a = np.dot(np.array([(i[1], i[0])
                                    for i in itertools.product(start_pos_vector[::-1],
                                                               start_pos_vector)]), t)
        # this is to convert the indices into bp ranges
        x = matrix_a[:, 1].reshape(n + 1, n + 1)
        y = matrix_a[:, 0].reshape(n + 1, n + 1)
        # plot
        im = ax.pcolormesh(x, y, np.flipud(matrix_c),
                           vmin=vmin, vmax=vmax, cmap=self.cmap, norm=self.norm)
        return im
