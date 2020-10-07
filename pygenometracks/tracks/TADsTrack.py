from . BedLikeTrack import BedLikeTrack, AROUND_REGION
from .. utilities import change_chrom_names
import numpy as np


class TADsTrack(BedLikeTrack):
    SUPPORTED_ENDINGS = ['.domain', '.domains', '.tad', '.tads']
    TRACK_TYPE = 'domains'
    OPTIONS_TXT = BedLikeTrack.OPTIONS_TXT + f"""
# If the bed file contains a column for color (column 9), then this color can be used by
# setting:
#color = bed_rgb
# if color is a valid colormap name (like RbBlGn), then the score (column 5) is mapped
# to the colormap.
# In this case, the the min_value and max_value for the score can be provided, otherwise
# the maximum score and minimum score found are used.
#color = RdYlBu
#min_value=0
#max_value=100
# If the color is simply a color name, then this color is used and the score is not considered.
color = darkblue
# optional. If not given it is guessed from the file ending.
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = dict({'max_value': None,
                                'min_value': None},
                               **BedLikeTrack.DEFAULTS_PROPERTIES)
    NECESSARY_PROPERTIES = BedLikeTrack.NECESSARY_PROPERTIES
    SYNONYMOUS_PROPERTIES = dict({'max_value': {'auto': None},
                                  'min_value': {'auto': None}},
                                 **BedLikeTrack.SYNONYMOUS_PROPERTIES)
    POSSIBLE_PROPERTIES = BedLikeTrack.POSSIBLE_PROPERTIES
    BOOLEAN_PROPERTIES = BedLikeTrack.BOOLEAN_PROPERTIES
    STRING_PROPERTIES = BedLikeTrack.STRING_PROPERTIES
    FLOAT_PROPERTIES = dict({'max_value': [- np.inf, np.inf],
                             'min_value': [- np.inf, np.inf]},
                            **BedLikeTrack.FLOAT_PROPERTIES)
    INTEGER_PROPERTIES = BedLikeTrack.INTEGER_PROPERTIES
    # The color can be a color or a colormap if bed_type is bed12 or 'bed_rgb'

    def __init__(self, *args, **kwarg):
        super(TADsTrack, self).__init__(*args, **kwarg)
        self.properties['display'] = 'triangles'

    def plot(self, ax, chrom_region, start_region, end_region):
        if chrom_region not in self.interval_tree.keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.interval_tree.keys():
                self.log.warning("*Warning*\nNo interval was found when "
                                 "overlapping with both "
                                 f"{chrom_region_before}:{start_region - AROUND_REGION}-{end_region + AROUND_REGION}"
                                 f" and {chrom_region}:{start_region - AROUND_REGION}-{end_region + AROUND_REGION}"
                                 " inside the bed file. "
                                 "This will generate an empty track!!\n")
                return
        chrom_region = self.check_chrom_str_bytes(self.interval_tree,
                                                  chrom_region)

        genes_overlap = \
            sorted(self.interval_tree[chrom_region][start_region:end_region])

        if self.properties['display'] == 'triangles':
            self.plot_triangles(ax, genes_overlap)
