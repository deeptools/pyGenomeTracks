from . GenomeTrack import GenomeTrack
from . BedTrack import BedTrack
import numpy as np

DEFAULT_BED_COLOR = '#1f78b4'


class TADsTrack(BedTrack):
    SUPPORTED_ENDINGS = ['.domain', '.domains', '.tad', '.tads']
    TRACK_TYPE = 'domains'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# If a bed file contains the exon
# structure (bed 12) then it is ignored.
# Triangles are only between start and end positions.
# If a gtf is given, it should end with gtf or gtf.gz or
# the type should be defined as gtf:
# type = gtf
# In the case of a gtf file, by default the transcript_name is used.
# If you want to use the gene_name:
# prefered_name = gene_name
# By default, the gtf is transformed to transcripts
# If you want to use see only one structure per gene
# merge_transcripts = true
# If the bed file contains a column for color (column 9), then this color can be used by
# setting:
# color = bed_rgb
#if color is a valid colormap name (like RbBlGn), then the score is mapped
# to the colormap.
# In this case, the the min_value and max_value for the score can be provided, otherwise
# the maximum score and minimum score found are used.
#color = RdYlBu
#min_value=0
#max_value=100
# If the color is simply a color name, then this color is used and the score is not considered.
# If color is set to none the domains will be transparent.
color = darkblue
# height of track in cm
height = 5
# optional: line_width
#line_width = 0.5
#border_color = black
# if you want to plot the track on top of the previous track. Options are 'yes' or 'share-y'. For the 'share-y'
# option the y axis values is shared between this plot and the overlay plot. Otherwise, each plot use its own scale
#overlay_previous = yes
# optional. If not given is guessed from the file ending.
file_type = {}
    """.format(TRACK_TYPE)

    DEFAULTS_PROPERTIES = {'fontsize': 12,
                           'orientation': None,
                           'color': DEFAULT_BED_COLOR,
                           'border_color': 'black',
                           'line_width': 0.5,
                           'prefered_name': 'transcript_name',
                           'merge_transcripts': False,
                           'max_value': None,
                           'min_value': None,
                           'type': 'bed'}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           'type': ['bed', 'gtf']}
    BOOLEAN_PROPERTIES = ['merge_transcripts']
    STRING_PROPERTIES = ['prefered_name', 'file', 'file_type',
                         'overlay_previous', 'orientation',
                         'title', 'color', 'border_color',
                         'type']
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
