from . GenomeTrack import GenomeTrack
from .. utilities import file_to_intervaltree, change_chrom_names
import numpy as np
import pysam

DEFAULT_BEDGRAPH_COLOR = '#a6cee3'


class BedGraphLikeTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = None
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT
    DEFAULTS_PROPERTIES = {}
    NECESSARY_PROPERTIES = []
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = []
    FLOAT_PROPERTIES = {}
    INTEGER_PROPERTIES = {}

    def __init__(self, properties_dict):
        GenomeTrack.__init__(self, properties_dict)
        self.load_file()

    def load_file(self):
        self.tbx = None
        # try to load a tabix file is available
        try:
            # from the tabix file is not possible to know the
            # global min and max
            self.tbx = pysam.TabixFile(self.properties['file'])
        except IOError:
            # load the file as an interval tree
            self.interval_tree, __, __ = file_to_intervaltree(self.properties['file'],
                                                              self.properties['region'])

        self.num_fields = None

    def _get_row_data(self, row, tbx_var='self.tbx'):
        """
        Returns the chrom, start, end and fields from either a tabix or a
        interval tree.
        Args:
            row: if tabix, the row comes from self.tbx.fetch otherwise
            comes from sorted(interval_tree[chrom] ...

        Returns:
            start, end, fields where values is a list

        """
        tbx = eval(tbx_var)
        if tbx is not None:
            fields = row.split("\t")
            values = fields[3:]
            start = int(fields[1])
            end = int(fields[2])

        else:
            values = row.data
            start = row.begin
            end = row.end

        # set the num_fields value
        # it is expected that the number of fields per row
        # is equal. This value is used for regions not covered
        # in the file and that should be represented as nans
        if self.num_fields is None:
            self.num_fields = len(values)
        return start, end, values

    def get_scores(self, chrom_region, start_region, end_region,
                   return_nans=True, tbx_var='self.tbx', inttree_var='self.interval_tree'):
        """
        Retrieves the score (or scores or whatever fields are in a bedgraph like file) and the positions
        for a given region.
        In case there is no item in the region. It returns [], []
        Args:
            chrom_region:
            start_region:
            end_region:
        Returns:
            tuple:
                scores_list, post_list
        """
        score_list = []
        pos_list = []
        tbx = eval(tbx_var)
        if tbx is not None:
            if chrom_region not in tbx.contigs:
                chrom_region_before = chrom_region
                chrom_region = change_chrom_names(chrom_region)
                if chrom_region not in tbx.contigs:
                    self.log.warning("*Warning*\nNeither "
                                     + chrom_region_before + " nor "
                                     + chrom_region + " exists as a "
                                     "chromosome name inside the bedgraph "
                                     "file. This will generate an empty "
                                     "track!!\n")
                    return score_list, pos_list

            chrom_region = self.check_chrom_str_bytes(tbx.contigs,
                                                      chrom_region)
            iterator = tbx.fetch(chrom_region, start_region, end_region)

        else:
            inttree = eval(inttree_var)
            if chrom_region not in list(inttree):
                chrom_region_before = chrom_region
                chrom_region = change_chrom_names(chrom_region)
                if chrom_region not in list(inttree):
                    self.log.warning("*Warning*\nNo interval was found when "
                                     "overlapping with both "
                                     f"{chrom_region_before}:{start_region}-{end_region}"
                                     f" and {chrom_region}:{start_region}-{end_region}"
                                     " inside the bedgraph file. "
                                     "This will generate an empty "
                                     "track!!\n")
                    return score_list, pos_list
            chrom_region = self.check_chrom_str_bytes(inttree, chrom_region)
            iterator = iter(sorted(inttree[chrom_region][start_region - 10000:end_region + 10000]))

        prev_end = start_region
        for row in iterator:
            start, end, values = self._get_row_data(row, tbx_var)
            # if the region is not consecutive with respect to the previous
            # nan values are added.
            if return_nans and prev_end < start:
                score_list.append(np.repeat(np.nan, self.num_fields))
                pos_list.append((prev_end, start))
            prev_end = end
            score_list.append(values)
            pos_list.append((start, end))

        # Add a last value if needed:
        if prev_end < end_region and return_nans:
            score_list.append(np.repeat(np.nan, self.num_fields))
            pos_list.append((prev_end, end_region))

        return score_list, pos_list

    def __del__(self):
        if self.tbx is not None:
            self.tbx.close()
