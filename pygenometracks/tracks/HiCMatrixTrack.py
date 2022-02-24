import scipy.sparse
from matplotlib import colors
import numpy as np
from . HiCMatrixLikeTrack import HiCMatrixLikeTrack
import logging
import itertools

DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


class HiCMatrixTrack(HiCMatrixLikeTrack):
    SUPPORTED_ENDINGS = ['.h5', '.cool', '.mcool']
    TRACK_TYPE = 'hic_matrix'
    OPTIONS_TXT = HiCMatrixLikeTrack.OPTIONS_TXT + f"""
# depth is the maximum distance that should be plotted.
# If it is more than 125% of the plotted region, it will
# be adjsted to this maximum value.
depth = 100000
file_type = {TRACK_TYPE}
    """
    DEFAULTS_PROPERTIES = dict({'depth': 100000},
                               **HiCMatrixLikeTrack.DEFAULTS_PROPERTIES)
    INTEGER_PROPERTIES = dict({'depth': [1, np.inf]},
                              **HiCMatrixLikeTrack.INTEGER_PROPERTIES)
    # The colormap can only be a colormap

    def plot(self, ax, chrom_region, region_start, region_end):

        continue_plotting, chrom_region = self.check_before_plotting(chrom_region, region_start, region_end)
        if not continue_plotting:
            return
        # expand region to plus depth on both sides
        # to avoid a 45 degree 'cut' on the edges

        # get bin id of start and end of region in given chromosome
        chr_start_id, chr_end_id = self.hic_ma.getChrBinRange(chrom_region)
        chr_start = self.hic_ma.cut_intervals[chr_start_id][1]
        chr_end = self.hic_ma.cut_intervals[chr_end_id - 1][2]
        start_bp = max(chr_start, region_start - self.properties['depth'])
        end_bp = min(chr_end, region_end + self.properties['depth'])
        idx = [idx for idx, x in enumerate(self.hic_ma.cut_intervals)
               if x[0] == chrom_region and x[1] >= start_bp and x[2] <= end_bp]
        if len(idx) == 0:
            self.log.warning("*Warning*\nThere is no data for the region "
                             "considered on the matrix. "
                             "This will generate an empty track!!\n")
            self.img = None
            return
        start_pos = [x[1] for i, x in enumerate(self.hic_ma.cut_intervals) if i in idx]
        # select only relevant matrix part
        matrix = self.hic_ma.matrix[idx, :][:, idx]
        # update the start_pos to add the last end:
        start_pos = tuple(list(start_pos) + [self.hic_ma.cut_intervals[idx[-1]][2]])
        # limit the 'depth' based on the length of the region being viewed

        region_len = region_end - region_start
        depth = min(self.properties['depth'], int(region_len * 1.25))
        # Need to be sure that you keep at least one bin even if the depth is
        # smaller than the binsize
        depth_in_bins = max(1, int(1.5 * region_len / self.hic_ma.getBinSize()))

        if depth < self.properties['depth']:
            log.warning(f"The depth was set to {self.properties['depth']} which is more than 125%"
                        " of the region plotted. The depth will be set "
                        f"to {depth}.\n")
            # remove from matrix all data points that are not visible.
            matrix = matrix - scipy.sparse.triu(matrix, k=depth_in_bins, format='csr')
        # Using todense will replace all nan values by 0.
        matrix = np.asarray(matrix.todense().astype(float))

        matrix = matrix * self.properties['scale_factor']

        if self.properties['transform'] == 'log1p':
            matrix += 1

        elif self.properties['transform'] in ['-log', 'log']:
            # We first replace 0 values by minimum values after 0
            mask = matrix == 0
            try:
                matrix[mask] = matrix[mask == False].min()
                matrix = np.log(matrix)
            except ValueError:
                self.log.info('All values are 0, no log applied.')
            else:
                if self.properties['transform'] == '-log':
                    matrix = - matrix

        if self.properties['max_value'] is not None:
            vmax = self.properties['max_value']

        else:
            # try to use a 'aesthetically pleasant' max value
            try:
                vmax = np.percentile(matrix.diagonal(1), 80)
            except Exception:
                vmax = None

        if self.properties['min_value'] is not None:
            vmin = self.properties['min_value']
        else:
            if depth_in_bins > matrix.shape[0]:
                # Make sure you keep one bin
                depth_in_bins = max(1, matrix.shape[0] - 5)

            # if the region length is large with respect to the chromosome length, the diagonal may have
            # very few values or none. Thus, the following lines reduce the number of bins until the
            # diagonal is at least length 5 but make sure you have at least one value:
            num_bins_from_diagonal = max(1, int(region_len / self.hic_ma.getBinSize()))
            for num_bins in range(0, num_bins_from_diagonal)[::-1]:
                distant_diagonal_values = matrix.diagonal(num_bins)
                if len(distant_diagonal_values) > 5:
                    break

            vmin = np.median(distant_diagonal_values)

        self.log.info("setting min, max values for track "
                      f"{self.properties['section_name']} to: "
                      f"{vmin}, {vmax}\n")

        if self.properties['transform'] == 'log1p':
            self.current_norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            self.current_norm = colors.Normalize(vmin=vmin, vmax=vmax)

        self.last_img_plotted = self.pcolormesh_45deg(ax, matrix, start_pos)
        if self.properties['rasterize']:
            self.last_img_plotted.set_rasterized(True)
        if self.properties['orientation'] == 'inverted':
            ax.set_ylim(depth, 0)
        else:
            ax.set_ylim(0, depth)

    def pcolormesh_45deg(self, ax, matrix_c, start_pos_vector):
        """
        Turns the matrix 45 degrees and adjusts the
        bins to match the actual start end positions.
        """
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
                           cmap=self.cmap, norm=self.current_norm)
        return im
