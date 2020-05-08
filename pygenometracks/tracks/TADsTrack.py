from . BedTrack import BedTrack
import numpy as np

DEFAULT_BED_COLOR = '#1f78b4'


class TADsTrack(BedTrack):
    SUPPORTED_ENDINGS = ['.domain', '.domains', '.tad', '.tads']
    TRACK_TYPE = 'domains'

    DEFAULTS_PROPERTIES = {'fontsize': 12,
                           'orientation': None,
                           'color': DEFAULT_BED_COLOR,
                           'border_color': 'black',
                           'line_width': 0.5,
                           # To remove in next 1.0
                           'prefered_name': 'transcript_name',
                           'merge_transcripts': False,
                           # End to remove
                           'max_value': None,
                           'min_value': None}
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
