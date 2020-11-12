from . BedGeneLikeTrack import BedGeneLikeTrack
import numpy as np


class BedTrack(BedGeneLikeTrack):
    SUPPORTED_ENDINGS = ['bed', 'bed3', 'bed4', 'bed5', 'bed6', 'bed8',
                         'bed9', 'bed12',
                         'bed.gz', 'bed3.gz', 'bed4.gz', 'bed5.gz', 'bed6.gz',
                         'bed9.gz', 'bed12.gz']
    TRACK_TYPE = 'bed'
    OPTIONS_TXT = BedGeneLikeTrack.OPTIONS_TXT + f"""
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
# optional. If not given is guessed from the file ending.
file_type = {TRACK_TYPE}
"""

    DEFAULTS_PROPERTIES = dict({'max_value': None,
                                'min_value': None},
                               **BedGeneLikeTrack.DEFAULTS_PROPERTIES)
    NECESSARY_PROPERTIES = BedGeneLikeTrack.NECESSARY_PROPERTIES
    SYNONYMOUS_PROPERTIES = dict({'max_value': {'auto': None},
                                  'min_value': {'auto': None}},
                                 **BedGeneLikeTrack.SYNONYMOUS_PROPERTIES)
    POSSIBLE_PROPERTIES = BedGeneLikeTrack.POSSIBLE_PROPERTIES
    BOOLEAN_PROPERTIES = BedGeneLikeTrack.BOOLEAN_PROPERTIES
    STRING_PROPERTIES = BedGeneLikeTrack.STRING_PROPERTIES
    FLOAT_PROPERTIES = dict({'max_value': [- np.inf, np.inf],
                             'min_value': [- np.inf, np.inf]},
                            **BedGeneLikeTrack.FLOAT_PROPERTIES)
    INTEGER_PROPERTIES = BedGeneLikeTrack.INTEGER_PROPERTIES
