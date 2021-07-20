from . GenomeTrack import GenomeTrack
from .. utilities import get_optimal_fontsize
import numpy as np
import os

# Color code for seqs
seq_color = {'A': 'red',
             'T': 'green',
             'G': 'blue',
             'C': 'black'}


def load_fasta(fastafile):
    """Returns a python dict { id : sequence } for the given .fasta file"""
    with open(os.path.realpath(fastafile), 'r') as filin:
        fasta = filin.read()
        fasta = fasta.split('>')[1:]
        outputdict = {x.split('\n')[0].strip(): "".join(x.split('\n')[1:]) for x in fasta}
    return outputdict


class FaTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.fa',".fasta"]
    TRACK_TYPE = 'fasta'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT

    DEFAULTS_PROPERTIES = {'line_width': 0.5,"height":1,"title":"Reference"}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ["file","title"]
    FLOAT_PROPERTIES = {'line_width': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        super(FaTrack, self).__init__(*args, **kwarg)
        self.ref = self.properties['file']
        self.seq = load_fasta(self.ref)

    def plot(self, ax, chrom_region, start_region, end_region):
        plotting_figure_width = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).width

        # The first constrain on the fontsize is the width
        ideal_fontsize = 1.4 * get_optimal_fontsize(plotting_figure_width,
                                              start_region, end_region)
        # The other constraint is the height
        # 1 point = 1/72 inch = height of character
        cm = 2.54
        max_fontsize = self.properties["height"] * 72 / cm

        # Let's take the biggest font possible with these constraints so that the figure is as readable as possible
        fontsize = min(ideal_fontsize,max_fontsize)

        seq_overlap = self.seq[chrom_region][start_region:end_region+1]
        for i in range(len(seq_overlap)):
            ax.text(i+start_region, 0.5, seq_overlap[i],
                    color=seq_color[seq_overlap[i]], verticalalignment='center',
                    horizontalalignment='center',fontsize=fontsize)

