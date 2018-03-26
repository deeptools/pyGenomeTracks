from .. utilities import to_string, to_bytes
import logging


class GenomeTrack(object):
    """
    The TrackPlot object is a holder for all tracks that are to be plotted.
    For example, to plot a bedgraph file a new class that extends TrackPlot
    should be created.

    It is expected that all GenomeTrack objects have a plot method.

    """
    SUPORTED_ENDINGS = []
    TRACK_TYPE = None
    OPTIONS_TXT = """
# title of track (plotted on the right side)
title =
# height of track in cm (ignored if the track is overlay on top the previous track)
height = 2
# if the track wants to be plotted upside-down:
# orientation = inverted
# if the track wants to be plotted on top of the previous track. Options are 'yes' or 'share-y'. For the 'share-y'
# option the y axis values is shared between this plot and the overlay plot. Otherwise, each plot use its own scale
#overlay previous = yes
"""

    def __init__(self, properties_dict):
        self.properties = properties_dict
        self.file_type = 'test'

        FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
        logging.basicConfig(format=FORMAT)
        log = logging.getLogger(__name__)
        log.setLevel(logging.DEBUG)
        self.log = log

    @staticmethod
    def change_chrom_names(chrom):
        """
        Changes UCSC chromosome names to ensembl chromosome names
        and vice versa.
        """
        # TODO: mapping from chromosome names like mithocondria is missing
        if chrom.startswith('chr'):
            # remove the chr part from chromosome name
            chrom = chrom[3:]
        else:
            # prefix with 'chr' the chromosome name
            chrom = 'chr' + chrom

        return chrom

    @staticmethod
    def check_chrom_str_bytes(iteratable_obj, p_obj):
        # determine type
        if isinstance(p_obj, list) and len(p_obj) > 0:
            type_ = type(p_obj[0])
        else:
            type_ = type(p_obj)
        if not isinstance(type(next(iter(iteratable_obj))), type_):
            if type(next(iter(iteratable_obj))) is str:
                p_obj = to_string(p_obj)
            elif type(next(iter(iteratable_obj))) in [bytes, np.bytes_]:
                p_obj = to_bytes(p_obj)
        return p_obj


