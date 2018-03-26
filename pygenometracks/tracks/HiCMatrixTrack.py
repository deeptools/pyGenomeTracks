import hicexplorer.HiCMatrix as HiCMatrix
import hicexplorer.utilities
import scipy.sparse
import copy
from matplotlib import cm
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
from past.builtins import zip

from . GenomeTrack import GenomeTrack

DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'


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
# if the track wants to be plotted upside-down:
# orientation = inverted
# optional if the values in the matrix need to be scaled the
# following parameter can be used. This is useful to plot multiple hic-matrices on the same scale
# scale factor = 1
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwargs):
        super(HiCMatrixTrack, self).__init__(*args, **kwargs)

        if self.properties['file'].endswith('.cool'):
            # just init the cooler matrix.
            self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'], color_only_init=True)
        else:
            self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'])

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

        new_intervals = hicexplorer.utilities.enlarge_bins(self.hic_ma.cut_intervals)
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

    def plot(self, ax, label_ax, chrom_region, region_start, region_end):
        self.cbar_ax = copy.copy(label_ax)

        chrom_sizes = self.hic_ma.get_chromosome_sizes()
        if chrom_region not in list(chrom_sizes):
            chrom_region = self.change_chrom_names(chrom_region)
            chrom_region = self.check_chrom_str_bytes(chrom_sizes, chrom_region)

        if region_end > chrom_sizes[chrom_region]:
            self.log.error("*Error*\nThe region to plot extends beyond the chromosome size. Please check.\n")
            self.log.error("{} size: {}. Region to plot {}-{}\n".format(chrom_region, chrom_sizes[chrom_region],
                                                                        region_start, region_end))

        if self.properties['file'].endswith('.cool'):
            # load now the region to be plotted
            pass

        # expand region to plus depth on both sides
        # to avoid a 45 degree 'cut' on the edges

        # get bin id of start and end of region in given chromosome
        chr_start_id, chr_end_id = self.hic_ma.getChrBinRange(chrom_region)
        chr_start = self.hic_ma.cut_intervals[chr_start_id][1]
        chr_end = self.hic_ma.cut_intervals[chr_end_id - 1][1]
        start_bp = max(chr_start, region_start - self.properties['depth'])
        end_bp = min(chr_end, region_end + self.properties['depth'])

        idx, start_pos = zip(*[(idx, x[1]) for idx, x in
                               enumerate(self.hic_ma.cut_intervals)
                               if x[0] == chrom_region and x[1] >= start_bp and x[2] <= end_bp])

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
        img = self.pcolormesh_45deg(ax, matrix, start_pos, vmax=vmax, vmin=vmin)
        img.set_rasterized(True)
        if self.plot_inverted:
            ax.set_ylim(depth, 0)
        else:
            ax.set_ylim(0, depth)

        self.cbar_ax.patch.set_alpha(0.0)
        try:
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
                cobar = plt.colorbar(img, ticks=tick_values, format=formatter, ax=self.cbar_ax, fraction=0.95)
            else:
                cobar = plt.colorbar(img, ax=self.cbar_ax, fraction=0.95)
            cobar.solids.set_edgecolor("face")
            cobar.ax.tick_params(labelsize='smaller')
            # cobar.ax.set_ylabel(self.properties['title'])

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

        except ValueError:
            pass

        label_ax.text(0.30, 0.5, self.properties['title'], size='large', verticalalignment='center')

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
