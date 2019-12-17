import argparse
import sys

import gffutils
import warnings

warnings.filterwarnings("ignore", message="It appears you have a gene feature"
                        " in your GTF file. You may want to use the "
                        "`disable_infer_genes` option to speed up database "
                        "creation")
warnings.filterwarnings("ignore", message="It appears you have a transcript "
                        "feature in your GTF file. You may want to use the "
                        "`disable_infer_transcripts` option to speed up "
                        "database creation")



def convert_gtf_to_bed(fn, fo, useGene, mergeTranscripts, ucsc):
    db = gffutils.create_db(fn, ':memory:')
    # For each transcript:
    prefered_name = "transcript_name"
    if useGene:
        prefered_name = "gene_name"
    if mergeTranscripts:
        all_items = db.features_of_type("gene", order_by='start')
    else:
        all_items = db.features_of_type("transcript", order_by='start')
    for tr in all_items:
        # The name would be the name of the transcript/gene if exists
        try:
            # First try to have it directly on the feature
            trName = tr.attributes[prefered_name][0]
        except KeyError:
            # Else try to guess the name of the transcript/gene from exons:
            try:
                trName = set([e.attributes[prefered_name][0]
                              for e in
                              db.children(tr,
                                          featuretype='exon',
                                          order_by='start')]).pop()
            except KeyError:
                # Else take the transcript/gene id
                trName = tr.id
        # If the cds is defined in the gtf,
        # use it to define the thick start and end
        # The gtf is 1-based closed intervalls and
        # bed are 0-based half-open so:
        # I need to remove one from each start
        try:
            # In case of multiple CDS (when there is one entry per gene)
            # I use the first one to get the start
            # and the last one to get the end (order_by=-start)
            cds_start = next(db.children(tr,
                                         featuretype='CDS',
                                         order_by='start')).start - 1
            cds_end = next(db.children(tr,
                                       featuretype='CDS',
                                       order_by='-start')).end
        except StopIteration:
            # If the CDS is not defined, then it is set to the start
            # as proposed here:
            # https://genome.ucsc.edu/FAQ/FAQformat.html#format1
            cds_start = tr.start - 1
            cds_end = tr.start - 1
        # Get all exons starts and lengths
        exons_starts = [e.start - 1
                        for e in
                        db.children(tr, featuretype='exon', order_by='start')]
        exons_length = [len(e)
                        for e in
                        db.children(tr, featuretype='exon', order_by='start')]
        # Rewrite the chromosome name if needed:
        chrom = tr.chrom
        if ucsc and chrom[0:3] != 'chr':
            chrom = 'chr' + chrom
        fo.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n" %
                 (chrom, tr.start - 1, tr.end, trName, 0, tr.strand,
                  cds_start, cds_end, "0", len(exons_starts),
                  ",".join([str(l) for l in exons_length]),
                  ",".join([str(s - (tr.start - 1)) for s in exons_starts])))


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description=("Convert a gtf to a bed12 with one entry"
                     " per transcript/gene"))
    parser.add_argument('input', default=None,
                        help="Input gtf file (can be gzip).")
    parser.add_argument('--output', default=sys.stdout,
                        type=argparse.FileType('w'),
                        help="Output bed12 file.")
    parser.add_argument('--useGene', action="store_true",
                        help="Use the gene name instead of the "
                        "transcript name.")
    parser.add_argument('--mergeTranscripts', action="store_true",
                        help="Merge all transcripts into a single "
                        "entry to have one line per gene.")
    parser.add_argument('--ucscformat', action="store_true",
                        help="If you want that all chromosome names "
                        "begin with 'chr'.")
    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    convert_gtf_to_bed(args.input, args.output, args.useGene,
                       args.mergeTranscripts, args.ucscformat)
