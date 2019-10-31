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

    def __init__(self, *args, **kwarg):
        super(TADsTrack, self).__init__(*args, **kwarg)
        self.properties['display'] = 'triangles'
