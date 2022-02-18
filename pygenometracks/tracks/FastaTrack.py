from . GenomeTrack import GenomeTrack
from .. utilities import get_optimal_fontsize, change_chrom_names
import numpy as np
import pyfaidx
import os

# Color code for seqs
seq_color = {'A': 'red',
             'T': 'green',
             'G': 'blue',
             'C': 'black',
             'N': 'grey'
             }


class FastaTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.fa', '.fasta']
    TRACK_TYPE = 'fasta'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT

    DEFAULTS_PROPERTIES = {}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['file', 'file_type', 'overlay_previous',
                         'orientation', 'title']
    FLOAT_PROPERTIES = {'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        super(FastaTrack, self).__init__(*args, **kwarg)
        try:
            self.seq = pyfaidx.Fasta(self.properties['file'])
        except pyfaidx.FastaIndexingError:
            self.seq = self.load_fasta(self.properties['file'])

    def load_fasta(self, fastafile):
        """Returns a python dict { id : sequence } for the given .fasta file"""
        with open(os.path.realpath(fastafile), 'r') as filin:
            fasta = filin.read()
            fasta = fasta.split('>')[1:]
            outputdict = {x.split('\n')[0].strip(): "".join(x.split('\n')[1:]) for x in fasta}
        return outputdict

    def plot_y_axis(self, ax, plot_axis):
        pass

    def plot(self, ax, chrom_region, start_region, end_region):

        if chrom_region not in self.seq.keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.seq.keys():
                self.log.warning("*Warning*\nNo record was found for "
                                 f"{chrom_region_before}"
                                 f" nor {chrom_region}"
                                 " inside the fasta file. "
                                 "This will generate an empty track!!\n")
                return

        plotting_figure_width = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).width

        # The first constrain on the fontsize is the width
        ideal_fontsize = 1.4 * get_optimal_fontsize(plotting_figure_width,
                                                    start_region, end_region)
        # The other constraint is the height
        # 1 point = 1/72 inch = height of character
        max_fontsize = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).height * 72

        # Let's take the biggest font possible with these constraints so that the figure is as readable as possible
        fontsize = min(ideal_fontsize, max_fontsize)

        if end_region > len(self.seq[chrom_region]):
            self.log.warning("*Warning*\nPlotting regions goes above"
                             " sequence length")
            end_region = len(self.seq[chrom_region])

        if type(self.seq) == pyfaidx.Fasta:
            seq_overlap = self.seq[chrom_region][start_region:end_region].seq
        else:
            seq_overlap = self.seq[chrom_region][start_region:end_region]

        # If the x-scale is inverted the complement is used:
        xleft, xright = ax.get_xlim()
        if xleft > xright:
            seq_overlap_correct = pyfaidx.complement(seq_overlap)
        else:
            seq_overlap_correct = seq_overlap

        for i, letter in enumerate(seq_overlap_correct):
            ax.text(i + start_region + 0.5, 0.5, letter,
                    color=seq_color[letter.upper()], verticalalignment='center',
                    horizontalalignment='center', fontsize=fontsize)
