from matplotlib import colors
import numpy as np
from . HiCMatrixLikeTrack import HiCMatrixLikeTrack
from .. utilities import change_chrom_names, get_region
import logging

DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


class HiCMatrixSquareTrack(HiCMatrixLikeTrack):
    # The supported endings are hidden to have only one track
    # per extension
    # SUPPORTED_ENDINGS = ['.h5', '.cool', '.mcool']
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = 'hic_matrix_square'
    OPTIONS_TXT = HiCMatrixLikeTrack.OPTIONS_TXT + f"""
# region2 is the region that should be plotted on the y axis.
# Default is the region on the x-axis
# By default the top is the start and the bottom is the end
# If orientation=inverted is used it is the contrary.
#region2 = X:3000000-3500000
file_type = {TRACK_TYPE}
    """
    DEFAULTS_PROPERTIES = dict({'region2': None},
                               **HiCMatrixLikeTrack.DEFAULTS_PROPERTIES)
    STRING_PROPERTIES = HiCMatrixLikeTrack.STRING_PROPERTIES + ['region2']

    def set_properties_defaults(self):
        # First I add region2 to region to get the matrix
        # Containing both region and region2:
        if 'region2' in self.properties:
            region2 = get_region(self.properties['region2'])
            self.properties['region'].append(region2)
        super(HiCMatrixSquareTrack, self).set_properties_defaults()

    def plot(self, ax, chrom_region, region_start, region_end):
        if len(self.hic_ma.matrix.data) == 0:
            self.log.warning("*Warning*\nThere is no data for the region "
                             "considered on the matrix. "
                             "This will generate an empty track!!\n")
            self.img = None
            return

        log.debug(f'chrom_region {chrom_region}, region_start {region_start}, region_end {region_end}')
        if chrom_region not in self.chrom_sizes:
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.chrom_sizes:
                self.log.warning("*Warning*\nNeither " + chrom_region_before
                                 + " nor " + chrom_region + " exists as a "
                                 "chromosome name on the matrix. "
                                 "This will generate an empty track!!\n")
                self.img = None
                return

        chrom_region = self.check_chrom_str_bytes(self.chrom_sizes, chrom_region)
        if region_end > self.chrom_sizes[chrom_region]:
            self.log.warning("*Warning*\nThe region to plot extends beyond the"
                             " chromosome size. Please check.\n"
                             f"{chrom_region} size: {self.chrom_sizes[chrom_region]}"
                             f". Region to plot {region_start}-{region_end}\n")

        # A chromosome may disappear if it was full of Nan and nan bins were masked:
        if chrom_region not in self.hic_ma.get_chromosome_sizes():
            self.log.warning("*Warning*\nThere is no data for the region "
                             "considered on the matrix. "
                             "This will generate an empty track!!\n")
            self.img = None
            return

        # get bin id of start and end of region in given chromosome
        chr_start_id_x, chr_end_id_x = self.hic_ma.getChrBinRange(chrom_region)
        chr_start_x = self.hic_ma.cut_intervals[chr_start_id_x][1]
        chr_end_x = self.hic_ma.cut_intervals[chr_end_id_x - 1][2]
        start_bp_x = max(chr_start_x, region_start - 3 * self.hic_ma.getBinSize())
        end_bp_x = min(chr_end_x, region_end + 3 * self.hic_ma.getBinSize())

        idx, start_pos = list(zip(*[(idx, x[1]) for idx, x in
                                    enumerate(self.hic_ma.cut_intervals)
                                    if x[0] == chrom_region and x[1] >= start_bp_x
                                    and x[2] <= end_bp_x]))
        # Process region2:
        if self.properties['region2'] is None:
            idx_y = idx
            start_pos_y = start_pos
            chrom_region_y, region_start_y, region_end_y = chrom_region, region_start, region_end
        else:
            chrom_region_y, region_start_y, region_end_y = get_region(self.properties['region2'])
            log.debug(f'On y: chrom_region: {chrom_region_y}, region_start {region_start_y}, region_end {region_end_y}')
            if chrom_region_y not in self.chrom_sizes:
                chrom_region_before = chrom_region_y
                chrom_region_y = change_chrom_names(chrom_region_y)
                if chrom_region_y not in self.chrom_sizes:
                    self.log.warning("*Warning*\nNeither " + chrom_region_before
                                     + " nor " + chrom_region_y + " exists as a "
                                     "chromosome name on the matrix. "
                                     "This will generate an empty track!!\n")
                    self.img = None
                    return

            chrom_region_y = self.check_chrom_str_bytes(self.chrom_sizes, chrom_region)
            if region_end_y > self.chrom_sizes[chrom_region_y]:
                self.log.warning("*Warning*\nThe region to plot extends beyond the"
                                 " chromosome size. Please check.\n"
                                 f"{chrom_region_y} size: {self.chrom_sizes[chrom_region_y]}"
                                 f". Region to plot {region_start_y}-{region_end_y}\n")

            # A chromosome may disappear if it was full of Nan and nan bins were masked:
            if chrom_region_y not in self.hic_ma.get_chromosome_sizes():
                self.log.warning("*Warning*\nThere is no data for the region "
                                 "considered on the matrix. "
                                 "This will generate an empty track!!\n")
                self.img = None
                return
            # get bin id of start and end of region2 in given chromosome
            chr_start_id_y, chr_end_id_y = self.hic_ma.getChrBinRange(chrom_region_y)
            chr_start_y = self.hic_ma.cut_intervals[chr_start_id_y][1]
            chr_end_y = self.hic_ma.cut_intervals[chr_end_id_y - 1][2]
            start_bp_y = max(chr_start_y, region_start_y - 3 * self.hic_ma.getBinSize())
            end_bp_y = min(chr_end_y, region_end_y + 3 * self.hic_ma.getBinSize())

            idx_y, start_pos_y = list(zip(*[(idx, x[1]) for idx, x in
                                          enumerate(self.hic_ma.cut_intervals)
                                          if x[0] == chrom_region_y and x[1] >= start_bp_y
                                          and x[2] <= end_bp_y]))
        print(f"{idx[0]} to {idx[-1]}")
        print(f"{idx_y[0]} to {idx_y[-1]}")
        # select only relevant matrix part
        matrix = self.hic_ma.matrix[idx, :][:, idx_y]
        # update the start_pos to add the last end:
        start_pos = tuple(list(start_pos) + [self.hic_ma.cut_intervals[idx[-1]][2]])
        start_pos_y = tuple(list(start_pos_y) + [self.hic_ma.cut_intervals[idx_y[-1]][2]])
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
            vmax = None

        if self.properties['min_value'] is not None:
            vmin = self.properties['min_value']
        else:
            vmin = None

        self.log.info("setting min, max values for track "
                      f"{self.properties['section_name']} to: "
                      f"{vmin}, {vmax}\n")

        if self.properties['transform'] == 'log1p':
            self.norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            self.norm = colors.Normalize(vmin=vmin, vmax=vmax)

        self.img = ax.pcolormesh(start_pos, start_pos_y, np.transpose(matrix),
                                 cmap=self.cmap, norm=self.norm, shading='flat')
        if self.properties['rasterize']:
            self.img.set_rasterized(True)
        if self.properties['orientation'] == 'inverted':
            ax.set_ylim(region_start_y, region_end_y)
        else:
            ax.set_ylim(region_end_y, region_start_y)
