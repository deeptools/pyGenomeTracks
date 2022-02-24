from hicmatrix import HiCMatrix
import hicmatrix.utilities
import scipy.sparse
from matplotlib import cm
import numpy as np
from . GenomeTrack import GenomeTrack
from .. utilities import change_chrom_names
import logging
import copy

DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


class HiCMatrixLikeTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = None
    OPTIONS_TXT = """
# The different options for color maps can be found here:
# https://matplotlib.org/users/colormaps.html
# the default color map is RdYlBu_r (_r) stands for reverse
# If you want your own colormap you can put the values of the color you want
# For example, colormap = ['blue', 'yellow', 'red']
# or colormap = ['white', (1, 0.88, .66), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
#colormap = RdYlBu_r
# height of track (in cm) can be given.
# Otherwise, the height is computed such that the proportions of the
# hic matrix are kept (e.g. the image does not appear shrink or extended)
# height = 10
# min_value and max_value refer to the contacts in the matrix.
#min_value =2.8
#max_value = 3.0
# the matrix can be transformed using the log1p (or log or -log, but zeros could be problematic)
transform = log1p
# show masked bins plots as white lines
# those bins that were not used during the correction
# the default is to extend neighboring bins to
# obtain an aesthetically pleasant output
show_masked_bins = false
# optional if the values in the matrix need to be scaled the
# following parameter can be used. This is useful to plot multiple hic-matrices on the same scale
# scale_factor = 1
# You can choose to keep the matrix as not rasterized
# (only used if you use pdf or svg output format) by using:
# rasterize = false
    """
    DEFAULTS_PROPERTIES = {'region': None,  # Cannot be set manually but is set by tracksClass
                           'orientation': None,
                           'show_masked_bins': False,
                           'scale_factor': 1,
                           'transform': 'no',
                           'max_value': None,
                           'min_value': None,
                           'rasterize': True,
                           'colormap': DEFAULT_MATRIX_COLORMAP}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           'transform': ['no', 'log', 'log1p', '-log']}
    BOOLEAN_PROPERTIES = ['show_masked_bins', 'rasterize']
    STRING_PROPERTIES = ['file', 'file_type', 'overlay_previous',
                         'orientation', 'transform',
                         'title', 'colormap']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'scale_factor': [- np.inf, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}
    # The colormap can only be a colormap

    def __init__(self, *args, **kwargs):
        super(HiCMatrixLikeTrack, self).__init__(*args, **kwargs)
        log.debug(f'FILE {self.properties}')

    def set_properties_defaults(self):
        super(HiCMatrixLikeTrack, self).set_properties_defaults()
        # Put default img to None for y axis
        self.last_img_plotted = None
        region = None
        if self.properties['region'] is not None:
            # We need to restrict it to a single region because
            # HiCMatrix does not accept more
            # We check if everything is on a single chrom:
            if len(set([r[0] for r in self.properties['region']])) == 1:
                chrom = self.properties['region'][0][0]
                start = min([r[1] for r in self.properties['region']])
                end = max([r[2] for r in self.properties['region']])
                # I extend of depth to avoid triangle effect in the plot
                if 'depth' in self.properties:
                    start = max(0, start - self.properties['depth'])
                    end += self.properties['depth']
                # I would like to find a way to sligthly above start-end
                # When there is no depth (hic_matrix_square)
                # Like 3 bins each direction but I don't manage
                # To think about a good way.
                region = [f"{chrom}:{start}-{end}"]
        # Cooler and thus HiCMatrix with cool file will raise an error if:
        # - the file is a cool file and:
        #    - the region goes over the chromosome size
        #   or
        #   - the chromosome is not part of the matrix

        # We need to change the log level because we don't want
        # the user to see all the errors raised during the try except
        logging.getLogger('hicmatrix').setLevel(logging.CRITICAL)
        try:
            self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'],
                                              pChrnameList=region)
        except ValueError as ve:
            if region is not None:
                if "Unknown sequence label" in str(ve):
                    rs = region[0].split(':')
                    chrom_region = rs[0]
                    chrom_region_before = chrom_region
                    chrom_region = change_chrom_names(chrom_region)
                    region = [f"{chrom_region}:{rs[1]}"]
                    try:
                        self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'],
                                                          pChrnameList=region)
                    except ValueError as ve2:
                        if "Unknown sequence label" in str(ve2):
                            self.log.warning("*Warning*\nNeither " + chrom_region_before
                                             + " nor " + chrom_region + " exists as a "
                                             "chromosome name on the matrix. "
                                             "This will generate an empty track!!\n")
                            self.hic_ma = HiCMatrix.hiCMatrix()
                            self.hic_ma.matrix = scipy.sparse.csr_matrix((0, 0))
                        elif "Genomic region out of bounds" in str(ve2):
                            region = [chrom_region]
                            self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'],
                                                              pChrnameList=region)
                        else:
                            raise ve2
                elif "Genomic region out of bounds" in str(ve):
                    region = [region[0].split(':')[0]]
                    self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'],
                                                      pChrnameList=region)
                else:
                    raise ve
            else:
                raise ve
        # We put back the log to warning
        logging.getLogger('hicmatrix').setLevel(logging.WARNING)

        if len(self.hic_ma.matrix.data) == 0:
            if region is None:
                # This is not due to a restriction of the matrix
                raise Exception(f"Matrix {self.properties['file']} is empty")
            else:
                return
        # We need to get the size before masking bins because
        # HiCMatrix>=v13 give smaller chromosome_sizes after:
        self.chrom_sizes = self.hic_ma.get_chromosome_sizes()
        if self.properties['show_masked_bins']:
            pass
        else:
            self.hic_ma.maskBins(self.hic_ma.nan_bins)

        # check that the matrix can be log transformed
        if self.properties['transform'] != 'no':
            if self.properties['transform'] == 'log1p':
                if self.hic_ma.matrix.data.min() + 1 <= 0:
                    raise Exception("\n*ERROR*\nMatrix contains values below - 1.\n"
                                    "log1p transformation can not be applied to \n"
                                    f"values in matrix: {self.properties['file']}")

            elif self.properties['transform'] in ['-log', 'log']:
                if self.hic_ma.matrix.data.min() < 0:
                    # For values not filled or equal to zero there will be a
                    # mask, they will be replaced by the minimum value after 0.
                    raise Exception("\n*ERROR*\nMatrix contains negative values.\n"
                                    "log transformation can not be applied to \n"
                                    f"values in matrix: {self.properties['file']}")

        new_intervals = hicmatrix.utilities.enlarge_bins(self.hic_ma.cut_intervals)
        self.hic_ma.interval_trees, self.hic_ma.chrBinBoundaries = \
            self.hic_ma.intervalListToIntervalTree(new_intervals)

        self.hic_ma.cut_intervals = new_intervals
        binsize = self.hic_ma.getBinSize()

        if 'depth' in self.properties:
            max_depth_in_bins = int(self.properties['depth'] / binsize)
            # If the depth is smaller than the binsize. It will display an empty plot
            if max_depth_in_bins < 1:
                self.log.warning(f"*Warning*\nThe depth({self.properties['depth']})"
                                 f" is smaller than binsize({binsize})"
                                 "This will generate an empty track!!\n")
                self.hic_ma.matrix = scipy.sparse.csr_matrix((0, 0))
                return

            self.reduce_matrix(max_depth_in_bins)

        self.process_color('colormap', colormap_possible=True,
                           colormap_only=True, default_value_is_colormap=True)

        self.cmap = copy.copy(cm.get_cmap(self.properties['colormap']))
        self.cmap.set_bad('black')

    def reduce_matrix(self, max_depth_in_bins):
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

    def check_before_plotting(self, chrom_region, region_start, region_end, suffix=''):
        if len(self.hic_ma.matrix.data) == 0:
            self.log.warning("*Warning*\nThere is no data for the region "
                             "considered on the matrix. "
                             "This will generate an empty track!!\n")
            self.last_img_plotted = None
            return False, chrom_region

        log.debug(f'{suffix}chrom_region {chrom_region}, region_start {region_start}, region_end {region_end}')
        if chrom_region not in self.chrom_sizes:
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.chrom_sizes:
                self.log.warning("*Warning*\nNeither " + chrom_region_before
                                 + " nor " + chrom_region + " exists as a "
                                 "chromosome name on the matrix. "
                                 "This will generate an empty track!!\n")
                self.last_img_plotted = None
                return False, chrom_region

        if region_start > self.chrom_sizes[chrom_region]:
            self.log.warning(f"*Warning*\nThe region to plot {suffix}starts beyond the"
                             " chromosome size. This will generate an empty track.\n"
                             f"{chrom_region} size: {self.chrom_sizes[chrom_region]}"
                             f". Region to plot {suffix}{region_start}-{region_end}\n")
            self.last_img_plotted = None
            return False, chrom_region

        if region_end > self.chrom_sizes[chrom_region]:
            self.log.warning(f"*Warning*\nThe region to plot {suffix}extends beyond the"
                             " chromosome size. Please check.\n"
                             f"{chrom_region} size: {self.chrom_sizes[chrom_region]}"
                             f". Region to plot {suffix}{region_start}-{region_end}\n")

        # A chromosome may disappear if it was full of Nan and nan bins were masked:
        if chrom_region not in self.hic_ma.get_chromosome_sizes():
            self.log.warning("*Warning*\nThere is no data for the region "
                             "considered on the matrix. "
                             "This will generate an empty track!!\n")
            self.last_img_plotted = None
            return False, chrom_region
        # Or it may be shortened:
        if region_start > self.hic_ma.get_chromosome_sizes()[chrom_region]:
            self.log.warning(f"*Warning*\nThe region to plot {suffix}starts beyond the"
                             " last bin with data on this chromosome."
                             " This will generate an empty track.\n"
                             f"{chrom_region} last bin: {self.hic_ma.get_chromosome_sizes()[chrom_region]}"
                             f". Region to plot {suffix}{region_start}-{region_end}\n")
            self.last_img_plotted = None
            return False, chrom_region

        return True, chrom_region

    def plot(self, ax, chrom_region, region_start, region_end):
        return

    def plot_y_axis(self, cbar_ax, plot_ax):
        if self.last_img_plotted is None:
            return

        GenomeTrack.plot_custom_cobar(self, cbar_ax)
