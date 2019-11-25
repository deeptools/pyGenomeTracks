from . BedTrack import BedTrack
import numpy as np

DEFAULT_BED_COLOR = '#1f78b4'


class TADsTrack(BedTrack):
    SUPPORTED_ENDINGS = ['.domain', '.domains', '.tad', '.tads']
    TRACK_TYPE = 'domains'

    DEFAULTS_PROPERTIES = {'fontsize': 12,
                           'orientation': None,
                           'color': DEFAULT_BED_COLOR,
                           'border color': 'black',
                           'interval_height': 100,  # This one is not defined in the documentation
                           'line width': 0.5,
                           'prefered name': 'transcript_name',
                           'merge transcripts': False,
                           'max_value': None,
                           'min_value': None}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted']}
    BOOLEAN_PROPERTIES = ['merge transcripts']
    STRING_PROPERTIES = ['prefered name', 'file', 'file_type',
                         'overlay previous', 'orientation',
                         'title']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'fontsize': [0, np.inf],
                        'interval_height': [0, np.inf],
                        'line width': [0, np.inf]}
    INTEGER_PROPERTIES = {}
    # The color can be a color or a colormap if bed_type is bed12 or 'bed_rgb'
    # border color can only be a color

    def __init__(self, *args, **kwarg):
        super(TADsTrack, self).__init__(*args, **kwarg)
        self.properties['display'] = 'triangles'
