from . BedTrack import BedTrack

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
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted']}
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    BOOLEAN_PROPERTIES = ['merge transcripts']
    STRING_PROPERTIES = ['prefered name']
    STRING_OR_NONE_PROPERTIES = []
    FLOAT_OR_NONE_PROPERTIES = {'max_value': [- np.inf, np.inf],
                                'min_value': [- np.inf, np.inf]}
    FLOAT_CONSTRAINED_PROPERTIES = {'fontsize': [0, np.inf],
                                    'interval_height': [0, np.inf],
                                    'line width': [0, np.inf]}
    INTEGER_OR_NONE_PROPERTIES = {}
    INTEGER_CONSTRAINED_PROPERTIES = {}
    # The color can be a color or a colormap if bed_type is bed12 or 'bed_rgb'
    # border color can only be a color

    def __init__(self, *args, **kwarg):
        super(TADsTrack, self).__init__(*args, **kwarg)
        self.properties['display'] = 'triangles'
