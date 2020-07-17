# -*- coding: utf-8 -*-
import collections

import gffutils
import warnings
from .utilities import InputError

warnings.filterwarnings("ignore", message="It appears you have a gene feature"
                        " in your GTF file. You may want to use the "
                        "`disable_infer_genes` option to speed up database "
                        "creation")
warnings.filterwarnings("ignore", message="It appears you have a transcript "
                        "feature in your GTF file. You may want to use the "
                        "`disable_infer_transcripts` option to speed up "
                        "database creation")
# In gffutils v0.10 they changed the error message:
warnings.filterwarnings("ignore", message="It appears you have a gene feature"
                        " in your GTF file. You may want to use the "
                        "`disable_infer_genes=True` option to speed up database "
                        "creation")
warnings.filterwarnings("ignore", message="It appears you have a transcript "
                        "feature in your GTF file. You may want to use the "
                        "`disable_infer_transcripts=True` option to speed up "
                        "database creation")


class ReadGtf(object):
    """
    Reads a gtf file.

    Example:
    gtf = ReadGtf("file.gtf")
    for interval in gtf:
        print interval.start

    """

    def __init__(self, file_path, prefered_name="transcript_name",
                 merge_transcripts=True):
        """
        :param file_path: the path of the gtf file
        :return:
        """

        self.file_type = 'bed12'

        # list of bed fields
        self.fields = ['chromosome', 'start', 'end',
                       'name', 'score', 'strand',
                       'thick_start', 'thick_end',
                       'rgb', 'block_count',
                       'block_sizes', 'block_starts']

        self.BedInterval = collections.namedtuple('BedInterval', self.fields)
        # I think the name which should be written
        # should be the transcript_name
        # But we can change it to gene_name
        self.prefered_name = prefered_name
        self.merge_transcripts = merge_transcripts

        # Will process the gtf to get one item per transcript:
        # This will create a database:
        try:
            self.db = gffutils.create_db(file_path, ':memory:')
        except ValueError as ve:
            if "No lines parsed" in str(ve):
                self.length = 0
                self.all_transcripts = open(file_path, 'r')
            else:
                raise InputError("This is not a gtf file.")
        else:
            if self.merge_transcripts:
                self.length = len([i for i in self.db.features_of_type("gene")])
                self.all_transcripts = self.db.features_of_type("gene",
                                                                order_by='start')
            else:
                self.length = len([i for i in self.db.features_of_type("transcript")])
                self.all_transcripts = self.db.features_of_type("transcript",
                                                                order_by='start')

    def __iter__(self):
        return self

    def __next__(self):
        """
        :return: bedInterval object
        """
        bed = self.get_bed_interval()

        return bed

    def get_bed_interval(self):
        """
        Process a transcript from the database,
        retrieve all the values and return
        a namedtuple object
        """
        tr = next(self.all_transcripts)
        # The name would be the prefered_name if exists
        try:
            trName = tr.attributes[self.prefered_name][0]
        except KeyError:
            # Else try to guess the prefered_name from exons:
            try:
                trName = set([e.attributes[self.prefered_name][0]
                              for e in
                              self.db.children(tr,
                                               featuretype='exon',
                                               order_by='start')]).pop()
            except KeyError:
                # Else take the transcript id
                trName = tr.id
        # If the cds is defined in the gtf,
        # use it to define the thick start and end
        # The gtf is 1-based closed intervalls
        # and bed are 0-based half-open so:
        # I need to remove one from each start
        try:
            cds_start = next(self.db.children(tr,
                                              featuretype='CDS',
                                              order_by='start')).start - 1
            cds_end = next(self.db.children(tr,
                                            featuretype='CDS',
                                            order_by='-start')).end
        except StopIteration:
            # If the CDS is not defined, then it is set to the start
            # as proposed here:
            # https://genome.ucsc.edu/FAQ/FAQformat.html#format1
            cds_start = tr.start - 1
            cds_end = tr.start - 1
        # Get all exons starts and end to get lengths
        exons_starts = [e.start - 1
                        for e in self.db.children(tr,
                                                  featuretype='exon',
                                                  order_by='start')]
        exons_ends = [e.end
                      for e in self.db.children(tr,
                                                featuretype='exon',
                                                order_by='start')]
        exons_length = [e - s for s, e in zip(exons_starts, exons_ends)]
        relative_exons_starts = [s - (tr.start - 1) for s in exons_starts]
        line_values = [tr.chrom, tr.start - 1, tr.end, trName, 0, tr.strand,
                       cds_start, cds_end, "0", len(exons_starts),
                       exons_length, relative_exons_starts]
        return self.BedInterval._make(line_values)
