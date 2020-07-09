from . BedTrack import BedTrack
from . GenomeTrack import GenomeTrack
import numpy as np

DEFAULT_BED_COLOR = '#1f78b4'


class TADsTrack(BedTrack):
    SUPPORTED_ENDINGS = ['.domain', '.domains', '.tad', '.tads']
    TRACK_TYPE = 'domains'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
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
# optional: line_width
#line_width = 0.5
# optional, default is black. To remove the border, simply set 'border_color' to none
#border_color = black
# optional. If not given it is guessed from the file ending.
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = {'orientation': None,
                           'color': DEFAULT_BED_COLOR,
                           'border_color': 'black',
                           'line_width': 0.5,
                           # To remove in next 1.0
                           'prefered_name': 'transcript_name',
                           'merge_transcripts': False,
                           # End to remove
                           'max_value': None,
                           'min_value': None,
                           'region': None}  # Cannot be set manually but is set by tracksClass
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted']}
    # In next 1.0: BOOLEAN_PROPERTIES = []
    BOOLEAN_PROPERTIES = ['merge_transcripts']
    STRING_PROPERTIES = ['file', 'file_type',
                         'overlay_previous', 'orientation',
                         'title', 'color', 'border_color',
                         # To remove in next 1.0
                         'prefered_name']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'fontsize': [0, np.inf],
                        'line_width': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}
    # The color can be a color or a colormap if bed_type is bed12 or 'bed_rgb'
    # border_color can only be a color

    def __init__(self, *args, **kwarg):
        super(TADsTrack, self).__init__(*args, **kwarg)
        self.properties['display'] = 'triangles'

    def set_properties_defaults(self):
        self.properties['fontsize'] = 12
        super(TADsTrack, self).set_properties_defaults()
        self.properties['global_max_row'] = False
