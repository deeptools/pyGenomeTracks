from . BedGeneLikeTrack import BedGeneLikeTrack
from . GenomeTrack import GenomeTrack
from . BedLikeTrack import AROUND_REGION
from .. readGtf import ReadGtf
from matplotlib import font_manager
from .. utilities import temp_file_from_intersect


class GtfTrack(BedGeneLikeTrack):
    SUPPORTED_ENDINGS = ['gtf', 'gtf.gz']
    TRACK_TYPE = 'gtf'
    OPTIONS_TXT = BedGeneLikeTrack.OPTIONS_TXT + f"""
# By default the transcript_name is used.
# If you want to use the gene_name:
# prefered_name = gene_name
# By default, the gtf is transformed to transcripts
# If you want to use see only one structure per gene
# merge_transcripts = true
# You can change the color of coding sequences by:
color = darkblue
# optional. If not given is guessed from the file ending.
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = dict({'prefered_name': 'transcript_name',
                                'merge_transcripts': False},
                               **BedGeneLikeTrack.DEFAULTS_PROPERTIES)
    NECESSARY_PROPERTIES = BedGeneLikeTrack.NECESSARY_PROPERTIES
    SYNONYMOUS_PROPERTIES = BedGeneLikeTrack.SYNONYMOUS_PROPERTIES
    POSSIBLE_PROPERTIES = BedGeneLikeTrack.POSSIBLE_PROPERTIES
    BOOLEAN_PROPERTIES = ['merge_transcripts'] + BedGeneLikeTrack.BOOLEAN_PROPERTIES
    STRING_PROPERTIES = ['prefered_name'] + BedGeneLikeTrack.STRING_PROPERTIES
    FLOAT_PROPERTIES = BedGeneLikeTrack.FLOAT_PROPERTIES
    INTEGER_PROPERTIES = BedGeneLikeTrack.INTEGER_PROPERTIES

    def set_properties_defaults(self):
        GenomeTrack.set_properties_defaults(self)
        self.fp = font_manager.FontProperties(size=self.properties['fontsize'])
        # to set the distance between rows
        self.row_scale = 2.3
        # For compatibility with other parent methods
        self.colormap = None
        # Contrary to bed it cannot be a colormap nor bed_rgb
        self.process_color('color', colormap_possible=False,
                           bed_rgb_possible=False,
                           default_value_is_colormap=False)

        # check if border_color and color_utr are colors
        for param in ['border_color', 'color_utr']:
            self.process_color(param, bed_rgb_possible=False)

    def get_bed_handler(self, plot_regions=None):
        if not self.properties['global_max_row']:
            # I do the intersection:
            file_to_open = temp_file_from_intersect(self.properties['file'],
                                                    plot_regions, AROUND_REGION)
        else:
            file_to_open = self.properties['file']

        bed_file_h = ReadGtf(file_to_open,
                             self.properties['prefered_name'],
                             self.properties['merge_transcripts'])
        total_length = bed_file_h.length
        return(bed_file_h, total_length)
