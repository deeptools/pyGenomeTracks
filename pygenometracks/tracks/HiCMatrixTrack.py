from hicmatrix import HiCMatrix
import hicmatrix.utilities
import scipy.sparse
from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

from . GenomeTrack import GenomeTrack

DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


class HiCMatrixTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.h5', '.cool' '.mcool']
    TRACK_TYPE = 'hic_matrix'
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
#min_value =2.8
#max_value = 3.0
# the matrix can be transformed using the log1 (or log, but zeros could be problematic)
transform = log1p
# show masked bins plots as white lines
# those bins that were not used during the correction
# the default is to extend neighboring bins to
# obtain an aesthetically pleasant output
show_masked_bins = no
# if you want to plot the track upside-down:
# orientation = inverted
# optional if the values in the matrix need to be scaled the
# following parameter can be used. This is useful to plot multiple hic-matrices on the same scale
# scale factor = 1
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwargs):
        super(HiCMatrixTrack, self).__init__(*args, **kwargs)

        log.debug('FILE {}'.format(self.properties))
        # log.debug('pRegion {}'.format(pRegion))
        region = None
        if self.properties['region'] is not None:
            if self.properties['region'][2] == 1e15:
                region = [str(self.properties['region'][0])]
            elif len(self.properties['region']) == 3:
                start = int(self.properties['region'][1]) - int(self.properties['depth'])
                if start < 0:
                    start = 0
                end = int(self.properties['region'][2]) + int(self.properties['depth'])

                region = [str(self.properties['region'][0]) + ':' + str(start) + '-' + str(end)]
        # try to open with end region + depth to avoid triangle effect in the plot
        # if it fails open it with given end region.
        try:
            self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'], pChrnameList=region)
        except Exception:
            region = [str(self.properties['region'][0]) + ':' + str(start) + '-' + str(self.properties['region'][2])]
            self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'], pChrnameList=region)

        if len(self.hic_ma.matrix.data) == 0:
            self.log.error("Matrix {} is empty".format(self.properties['file']))
            exit(1)
        if 'show_masked_bins' in self.properties and self.properties['show_masked_bins'] == 'yes':
            pass
        else:
            self.hic_ma.maskBins(self.hic_ma.nan_bins)

        # check that the matrix can be log transformed
        if 'transform' in self.properties:
            if self.properties['transform'] == 'log1p':
                if self.hic_ma.matrix.data.min() + 1 < 0:
                    self.log.error("\n*ERROR*\nMatrix contains negative values.\n"
                                   "log1p transformation can not be applied to \n"
                                   "values in matrix: {}".format(self.properties['file']))
                    exit(1)

            elif self.properties['transform'] == '-log':
                if self.hic_ma.matrix.data.min() < 0:
                    self.log.error("\n*ERROR*\nMatrix contains negative values.\n"
                                   "log(-1 * <values>) transformation can not be applied to \n"
                                   "values in matrix: {}".format(self.properties['file']))
                    exit(1)

            elif self.properties['transform'] == 'log':
                if self.hic_ma.matrix.data.min() < 0:
                    self.log.error("\n*ERROR*\nMatrix contains negative values.\n"
                                   "log transformation can not be applied to \n"
                                   "values in matrix: {}".format(self.properties['file']))
                    exit(1)

        new_intervals = hicmatrix.utilities.enlarge_bins(self.hic_ma.cut_intervals)
        self.hic_ma.interval_trees, self.hic_ma.chrBinBoundaries = \
            self.hic_ma.intervalListToIntervalTree(new_intervals)

        self.hic_ma.cut_intervals = new_intervals
        binsize = self.hic_ma.getBinSize()
        max_depth_in_bins = int(self.properties['depth'] / binsize)

        # work only with the lower matrix
        # and remove all pixels that are beyond
        # 2 * max_depth_in_bis which are not required
        # (this is done by subtracting a second sparse matrix
        # that contains only the lower matrix that wants to be removed.
        limit = 2 * max_depth_in_bins
        self.hic_ma.matrix = scipy.sparse.triu(self.hic_ma.matrix, k=0, format='csr') - \
            scipy.sparse.triu(self.hic_ma.matrix, k=limit, format='csr')
        self.hic_ma.matrix.eliminate_zeros()

        # fill the main diagonal, otherwise it looks
        # not so good. The main diagonal is filled
        # with an array containing the max value found
        # in the matrix
        if sum(self.hic_ma.matrix.diagonal()) == 0:
            self.log.info("Filling main diagonal with max value because it empty and looks bad...\n")
            max_value = self.hic_ma.matrix.data.max()
            main_diagonal = scipy.sparse.dia_matrix(([max_value] * self.hic_ma.matrix.shape[0], [0]),
                                                    shape=self.hic_ma.matrix.shape)
            self.hic_ma.matrix = self.hic_ma.matrix + main_diagonal

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

        log.debug('chrom_region {}, region_start {}, region_end {}'.format(chrom_region, region_start, region_end))
        chrom_sizes = self.hic_ma.get_chromosome_sizes()
        if chrom_region not in chrom_sizes:
            chrom_region_before = chrom_region
            chrom_region = self.change_chrom_names(chrom_region)
            if chrom_region not in chrom_sizes:
                self.log.error("*Error*\nNeither " + chrom_region_before + " "
                               "nor " + chrom_region + " exits as a chromosome"
                               " name on the matrix.\n")
                return

        chrom_region = self.check_chrom_str_bytes(chrom_sizes, chrom_region)
        if region_end > chrom_sizes[chrom_region]:
            self.log.error("*Error*\nThe region to plot extends beyond the chromosome size. Please check.\n")
            self.log.error("{} size: {}. Region to plot {}-{}\n".format(chrom_region, chrom_sizes[chrom_region],
                                                                        region_start, region_end))

        # if self.properties['file'].endswith('.cool'):
        #     # load now the region to be plotted
        #     pass

        # expand region to plus depth on both sides
        # to avoid a 45 degree 'cut' on the edges

        # get bin id of start and end of region in given chromosome
        chr_start_id, chr_end_id = self.hic_ma.getChrBinRange(chrom_region)
        chr_start = self.hic_ma.cut_intervals[chr_start_id][1]
        chr_end = self.hic_ma.cut_intervals[chr_end_id - 1][1]
        start_bp = max(chr_start, region_start - self.properties['depth'])
        end_bp = min(chr_end, region_end + self.properties['depth'])

        idx, start_pos = list(zip(*[(idx, x[1]) for idx, x in
                                    enumerate(self.hic_ma.cut_intervals)
                                    if x[0] == chrom_region and x[1] >= start_bp and x[2] <= end_bp]))

        idx = idx[0:-1]
        # select only relevant matrix part
        matrix = self.hic_ma.matrix[idx, :][:, idx]
        # limit the 'depth' based on the length of the region being viewed

        region_len = region_end - region_start
        depth = min(self.properties['depth'], int(region_len * 1.25))
        depth_in_bins = int(1.5 * region_len / self.hic_ma.getBinSize())

        if depth < self.properties['depth']:
            # remove from matrix all data points that are not visible.
            matrix = matrix - scipy.sparse.triu(matrix, k=depth_in_bins, format='csr')
        matrix = np.asarray(matrix.todense().astype(float))
        if 'scale factor' in self.properties:
            matrix = matrix * self.properties['scale factor']

        if 'transform' in self.properties:
            if self.properties['transform'] == 'log1p':
                matrix += 1
                self.norm = colors.LogNorm()

            elif self.properties['transform'] == '-log':
                mask = matrix == 0
                try:
                    matrix[mask] = matrix[mask == False].min()
                    matrix = -1 * np.log(matrix)
                except ValueError:
                    self.log.info('All values are 0, no log applied.')

            elif self.properties['transform'] == 'log':
                mask = matrix == 0
                try:
                    matrix[mask] = matrix[mask == False].min()
                    matrix = np.log(matrix)
                except ValueError:
                    self.log.info('All values are 0, no log applied.')

        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            vmax = self.properties['max_value']

        else:
            # try to use a 'aesthetically pleasant' max value
            vmax = np.percentile(matrix.diagonal(1), 80)

        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            vmin = self.properties['min_value']
        else:
            if depth_in_bins > matrix.shape[0]:
                depth_in_bins = matrix.shape[0] - 5

            # if the region length is large with respect to the chromosome length, the diagonal may have
            # very few values or none. Thus, the following lines reduce the number of bins until the
            # diagonal is at least length 5
            num_bins_from_diagonal = int(region_len / self.hic_ma.getBinSize())
            for num_bins in range(0, num_bins_from_diagonal)[::-1]:
                distant_diagonal_values = matrix.diagonal(num_bins)
                if len(distant_diagonal_values) > 5:
                    break

            vmin = np.median(distant_diagonal_values)

        self.log.info("setting min, max values for track {} to: {}, {}\n".
                      format(self.properties['section_name'], vmin, vmax))
        self.img = self.pcolormesh_45deg(ax, matrix, start_pos, vmax=vmax, vmin=vmin)
        self.img.set_rasterized(True)
        if self.plot_inverted:
            ax.set_ylim(depth, 0)
        else:
            ax.set_ylim(0, depth)

    def plot_y_axis(self, cbar_ax, plot_ax):

        if 'transform' in self.properties and \
                self.properties['transform'] in ['log', 'log1p']:
            # get a useful log scale
            # that looks like [1, 2, 5, 10, 20, 50, 100, ... etc]

            # The following code is problematic with some versions of matplotlib.
            # Should be uncommented once the problem is clarified
            from matplotlib.ticker import LogFormatter
            formatter = LogFormatter(10, labelOnlyBase=False)
            aa = np.array([1, 2, 5])
            tick_values = np.concatenate([aa * 10 ** x for x in range(10)])
            try:
                cobar = plt.colorbar(self.img, ticks=tick_values, format=formatter, ax=cbar_ax, fraction=0.95)
            except AttributeError:
                return
        else:
            try:
                cobar = plt.colorbar(self.img, ax=cbar_ax, fraction=0.95)
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
