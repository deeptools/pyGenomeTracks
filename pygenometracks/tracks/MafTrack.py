from . GenomeTrack import GenomeTrack
from .. utilities import change_chrom_names
from Bio import AlignIO
import matplotlib
from matplotlib import font_manager
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
from intervaltree import IntervalTree, Interval
import numpy as np
from tqdm import tqdm
import re

# This is to match the database and chromosome from src:
src_re = re.compile(r'^(?P<database>\w+)\.(?P<chromosome>\w+)$')


class MafTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['maf']
    TRACK_TYPE = 'maf'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# In addition to the maf file the reference is required:
# For example
#reference = mm10.chr2
# To speed the access to specific region
# The maf file needs an index.
# The default is the file
# followed by 'index'. Alternatively another
# file can be specified.
# If it does not exists it will be created:
#file_index =
# Set colors
#color_identical = black
#color_mismatch = grey
#color_gap = lightgrey
# not used for the moment: line_width
#line_width = 0.5
# optional: the species order
#species_order = hg18 panTro2
# optional if species_order is specified (don't use space in names)
#species_labels = human chimpanzee
# optional: If not given is guessed from the file ending.
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = {'file_index': None,
                           'orientation': None,
                           'color_identical': 'black',
                           'color_mismatch': 'grey',
                           'color_gap': 'lightgrey',
                           'line_width': 0.5,
                           'species_order': None,
                           'species_labels': None}
    NECESSARY_PROPERTIES = ['file', 'reference']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           }
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['file', 'file_type',
                         'overlay_previous', 'orientation',
                         'title', 'file_index', 'color_identical',
                         'color_mismatch',
                         'color_gap', 'reference',
                         'species_order', 'species_labels']
    FLOAT_PROPERTIES = {'line_width': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        super(MafTrack, self).__init__(*args, **kwarg)
        self.ref = self.properties['reference']
        # If no region is specified use classical
        self.idx = None
        self.database = None
        self.chromosome = None
        if self.properties['region'] is not None:
            # Read or create the index
            self.idx = AlignIO.MafIO.MafIndex(self.properties['file_index'],
                                              self.properties['file'],
                                              self.ref)
            # Process the reference
            match_ref = src_re.match(self.ref)
            if match_ref is not None:
                self.database = match_ref.group('database')
                self.chromosome = match_ref.group('chromosome')
        # Process the species_order and species_labels:
        self.species, self.labels = self.process_species_user()
        # Then get the interval_tree
        self.interval_tree = self.process_maf(self.properties['region'])

    def set_properties_defaults(self):
        super(MafTrack, self).set_properties_defaults()
        # Process colors
        for p in ['color_identical', 'color_mismatch', 'color_gap']:
            self.process_color(p)
        # Process file_index
        if self.properties['file_index'] is None:
            self.properties['file_index'] = self.properties['file'] + 'index'
        # to set the distance between rows
        self.row_scale = 1.3

    def process_species_user(self):
        my_species = None
        my_labels = None
        if self.properties['species_order'] is not None:
            my_species = self.properties['species_order'].split()
        if self.properties['species_labels'] is not None:
            if self.properties['species_order'] is None:
                self.log.warning("species_labels was specified while species_order is not."
                                 " The labels will be ignored.\n")
                my_labels = None
            else:
                my_labels = self.properties['species_labels'].split()
                if len(my_labels) != len(my_species):
                    self.log.warning(f"species_labels was specified but its size ({len(my_labels)}) is different from "
                                     f"what is specified in species_order ({len(my_species)})."
                                     " The labels will be ignored.\n")
                    my_labels = None
        if my_species is not None and my_labels is None:
            my_labels = my_species
        return(my_species, my_labels)
                    
    def process_maf(self, plot_regions=None):

        valid_blocks = 0
        interval_tree = IntervalTree()

        if self.idx is not None:
            # First check that the chromosome plotted
            # is compatible with the ref
            # and extract valid starts/ends
            starts = []
            ends = []
            valid_chrom = None
            if self.chromosome is not None:
                valid_chrom = [self.chromosome, change_chrom_names(self.chromosome)]
            for r in plot_regions:
                if self.chromosome is None or r[0] in valid_chrom:
                    # We only add start and ends to regions to plot in
                    # the correct chromosome
                    starts.append(r[1])
                    ends.append(r[2])
            if len(starts) == 0:
                iterator = []
            else:
                iterator = self.idx.search(starts, ends)
        else:
            iterator = AlignIO.parse(self.properties['file'], "maf")
        # Then go through blocks
        all_ids = {}
        for block in tqdm(iterator):
            # Process the block
            # First get the ref
            ref_alg = None
            for alg in block:
                if alg.id == self.ref:
                    ref_alg = alg
                    break
            # If there is nothing for the ref go to next block
            if ref_alg is None:
                continue
            ref_seq = str(ref_alg.seq.upper())
            gaps = [m.start() for m in re.finditer('(?=-)', ref_seq)]
            ref_seq_no_gap = ''.join([l for i, l in enumerate(ref_seq)
                                      if i not in gaps])
            # First naive:
            for alg in block:
                if alg.id != self.ref:
                    alg_id = alg.id
                    if alg_id not in all_ids:
                        match_alg = src_re.match(alg_id)
                        all_ids[alg_id] = match_alg
                    alg_seq = str(alg.seq.upper())
                    alg_seq_no_gap = ''.join([l for i, l in enumerate(alg_seq)
                                              if i not in gaps])
                    last_status = None
                    last_end = None
                    # This might need some improvement
                    for i, (lr, la) in enumerate(zip(ref_seq_no_gap, alg_seq_no_gap),
                                                 start=ref_alg.annotations['start']):
                        # print(f"{i}:{lr}, {la}")
                        cur_status = self.compare_letters(lr, la)
                        # print(cur_status)
                        if cur_status != last_status:
                            if last_status is not None:
                                interval_tree.add(Interval(last_end, i, {'id':alg_id, 'status':last_status}))
                                valid_blocks += 1
                            last_end = i
                            last_status = cur_status
                    if last_status is not None:
                        interval_tree.add(Interval(last_end, i + 1, {'id':alg_id, 'status':last_status}))
                        valid_blocks += 1

        if valid_blocks == 0:
            self.log.warning("No valid blocks were found in file "
                             f"{self.properties['file']}.\n")
        self.process_ids(all_ids)

        return interval_tree

    def process_ids(self, all_ids):
        ''' 
        Process the all_ids to get a correspondence with the species
        '''
        self.id_to_species = {}
        found_species = []
        for id in all_ids:
            if all_ids[id] is None:
                self.id_to_species[id] = id
                found_species.append(id)
            else:
                genome = all_ids[id].group('database')
                self.id_to_species[id] = genome
                found_species.append(genome)
        if self.species is None:
            self.species = []
            self.labels = []
        for new_species in found_species:
            if new_species not in self.species:
                self.species.append(new_species)
                self.labels.append(new_species)

        # Give row number
        self.species_y = {}
        for y, species in enumerate(self.species):
            self.species_y[species] = y
        self.max_y = y

    def plot(self, ax, chrom_region, start_region, end_region):
        if self.chromosome is not None:
            valid_chrom = [self.chromosome, change_chrom_names(self.chromosome)]
            if chrom_region not in valid_chrom:
                self.log.warning("The plotting chromosome is not the reference."
                                 " Nothing will be plotted.")
        blocks_overlap = \
            sorted(self.interval_tree[start_region:end_region])

        for block in blocks_overlap:
            ypos = self.species_y[self.id_to_species[block.data['id']]] * self.row_scale
            ax.add_patch(Rectangle((block.begin, ypos),
                         block.end - block.begin, 1,
                         edgecolor="none",
                         facecolor=self.properties[f'color_{block.data["status"]}']))

        epsilon = 0.08
        ymax = - epsilon
        ymin = self.row_scale * self.max_y + (1 + epsilon)

        self.log.debug(f"ylim {ymin},{ymax}")
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

    def plot_y_axis(self, ax, plot_axis):
        for y, label in enumerate(self.labels):
            ax.text(0, y * self.row_scale + 0.5, label,
                    verticalalignment='center',
                    horizontalalignment='right', wrap=True)
        ax.set_ylim(*plot_axis.get_ylim())

    @staticmethod
    def compare_letters(lr, la):
        """
            Compare a reference letter lr with another letter lr
            We assume lr cannot be '-' while la can.
        """
        if lr == la:
            return("identical")
        if la == '-':
            return("gap")
        return("mismatch")