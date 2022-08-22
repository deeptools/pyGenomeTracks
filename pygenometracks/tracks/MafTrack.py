from . GenomeTrack import GenomeTrack
from .. utilities import change_chrom_names, InputError, get_optimal_fontsize
import bx.align.maf
import bx.interval_index_file
import bx.seq
import matplotlib as mpl
from matplotlib.patches import Rectangle
import numpy as np
from tqdm import tqdm
import tempfile
import re

# Color code for seqs
seq_color = {'A': 'red',
             'T': 'green',
             'G': 'blue',
             'C': 'black',
             'N': 'grey'}


class MafTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['maf']
    TRACK_TYPE = 'maf'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# In addition to the maf file the reference genome is required:
# For example
#reference = mm10
# To speed the access to specific region
# The maf file needs an index.
# The default is the file
# followed by '.index'. Alternatively another
# file can be specified.
# If it does not exists it will be created:
#file_index =
# Set colors
#color_identical = black
#color_mismatch = grey
#color_gap = lightgrey
# Set line_width
#line_width = 0.5
# optional: the species order
#species_order = hg18 panTro2
# optional if species_order is specified (don't use space in names)
#species_labels = human chimpanzee
# optional if species_order is specified and they are the only one you want:
#species_order_only = true
# optional if you want to see the DNA sequence of the ref
#display_ref_seq = true
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
                           'species_labels': None,
                           'species_order_only': False,
                           'display_ref_seq': False}
    NECESSARY_PROPERTIES = ['file', 'reference']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           }
    BOOLEAN_PROPERTIES = ['species_order_only', 'display_ref_seq']
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
        # Read or create the index
        try:
            self.idx = bx.align.maf.Indexed(self.properties['file'],
                                            self.properties['file_index'],
                                            keep_open=True,
                                            parse_e_rows=True)
        except FileNotFoundError:
            self.idx = None
            self.write_maf_index()
        else:
            if not self.maf_index_has_good_ref():
                new_index = tempfile.NamedTemporaryFile(delete=False)
                self.log.warning(f"The index {self.properties['file_index']}"
                                 f" is not for {self.ref}."
                                 f" Will generate a new one to {new_index.name}")
                self.properties['file_index'] = new_index.name
                self.idx = None
                self.write_maf_index()
        if self.idx is None:
            self.idx = bx.align.maf.Indexed(self.properties['file'],
                                            self.properties['file_index'],
                                            keep_open=True,
                                            parse_e_rows=True)
        # Process the species_order and species_labels:
        self.species, self.labels = self.process_species_user()
        # Initialize current_labels:
        self.current_labels = None
        # Give row number
        self.species_y = {}
        self.max_y = 0
        if self.species is not None:
            for species in self.species:
                self.species_y[species] = self.max_y
                self.max_y += 1
        # self.log.debug(f"self.max_y {self.max_y}")
        # Check incompatibilities
        if self.species is None and self.properties['species_order_only']:
            raise InputError("species_order_only set to true while species_order not specified.")

    def set_properties_defaults(self):
        super(MafTrack, self).set_properties_defaults()
        # Process colors
        for p in ['color_identical', 'color_mismatch', 'color_gap']:
            self.process_color(p)
        # Process file_index
        if self.properties['file_index'] is None:
            self.properties['file_index'] = self.properties['file'] + '.index'
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
        return my_species, my_labels

    def plot(self, ax, chrom_region, start_region, end_region):
        ref_in_index_array = self.ref_chrom_in_maf_index(chrom_region)
        if len(ref_in_index_array) == 0:
            self.log.warning("The plotting chromosome is not the reference."
                             " Nothing will be plotted.")
            return
        ref_in_index = ref_in_index_array[0]
        epsilon = 0.08
        ymax = 0
        valid_blocks = 0
        # Initiate lists and dicts
        if self.species is not None:
            current_species = self.species.copy()
        else:
            current_species = None
        if self.labels is not None:
            self.current_labels = self.labels.copy()
        else:
            self.current_labels = None
        current_species_y = self.species_y.copy()
        current_max_y = self.max_y

        if self.properties['display_ref_seq']:
            ref_seq_dic = {}
        for block in tqdm(self.idx.get_as_iterator(ref_in_index, start_region, end_region)):
            # I need to slice if needed:
            ref = block.get_component_by_src(ref_in_index)
            slice_start = max(start_region, ref.get_forward_strand_start())
            slice_end = min(end_region, ref.get_forward_strand_end())
            if slice_start != ref.get_forward_strand_start() or \
               slice_end != ref.get_forward_strand_end():
                sliced = block.slice_by_component(ref, slice_start, slice_end)
                ref = sliced.get_component_by_src(ref_in_index)
            else:
                sliced = block
            # I need to reverse if the ref is not + strand:
            if ref.strand == "-":
                sliced = sliced.reverse_complement()
                ref = sliced.get_component_by_src(ref_in_index)
            # Get the gap positions in the ref:
            ref_seq = ref.text.upper()
            gaps = [m.start() for m in re.finditer('(?=-)', ref_seq)]
            ref_seq_no_gap = ''.join([l for i, l in enumerate(ref_seq)
                                      if i not in gaps])
            # Store the sequence if required:
            if self.properties['display_ref_seq']:
                # Check the orientation:
                xleft, xright = ax.get_xlim()
                if xleft > xright:
                    # I need to complement (not reversed because I keep the coordinates)
                    ref_seq_no_gap_correct = ref_seq_no_gap.translate(bx.seq.DNA_COMP)
                else:
                    ref_seq_no_gap_correct = ref_seq_no_gap
                for i, lr in enumerate(ref_seq_no_gap_correct,
                                       start=ref.get_forward_strand_start()):
                    ref_seq_dic[i] = lr
            for c in sliced.components:
                # We only plot the non-ref:
                if c.src != ref_in_index:
                    assembly = c.src.split(".")[0]
                    # Check if we plot this assembly:
                    if self.properties["species_order_only"] and assembly not in self.species:
                        continue
                    # Initiate self.current_species if needed:
                    if current_species is None:
                        current_species = []
                        self.current_labels = []
                    # Add this assembly if needed:
                    if assembly not in current_species:
                        current_species.append(assembly)
                        self.current_labels.append(assembly)
                        current_species_y[assembly] = current_max_y
                        current_max_y += 1
                        # self.log.debug(f"self.max_y {self.max_y}")
                    # Get the position:
                    ypos = current_species_y[assembly] * self.row_scale
                    if not c.empty:
                        # Get the sequence to compare with ref:
                        c_seq = c.text.upper()
                        c_seq_no_gap = ''.join([l for i, l in enumerate(c_seq)
                                                if i not in gaps])
                        last_status = None
                        last_end = None
                        # This might need some improvement
                        for i, (lr, la) in enumerate(zip(ref_seq_no_gap, c_seq_no_gap),
                                                     start=ref.get_forward_strand_start()):
                            cur_status = self.compare_letters(lr, la)
                            if cur_status != last_status:
                                if last_status is not None:
                                    # I plot it:
                                    ax.add_patch(Rectangle((last_end, ypos),
                                                           i - last_end, 1,
                                                           edgecolor="none",
                                                           facecolor=self.properties[f'color_{last_status}']))
                                    valid_blocks += 1
                                last_end = i
                                last_status = cur_status
                        if last_status is not None:
                            # I plot it:
                            ax.add_patch(Rectangle((last_end, ypos),
                                                   ref.get_forward_strand_end() - last_end, 1,
                                                   edgecolor="none",
                                                   facecolor=self.properties[f'color_{last_status}']))
                            valid_blocks += 1
                    else:
                        if c.synteny_empty == "C":
                            # the sequence before and after is contiguous
                            # implying that this region was either deleted
                            # in the source or inserted in the reference sequence.
                            # The browser draws a single line or a "-" in base mode in these blocks.
                            ax.plot([ref.get_forward_strand_start(), ref.get_forward_strand_end()],
                                    [ypos + 0.5, ypos + 0.5], color="black", linewidth=self.properties['line_width'])
                        elif c.synteny_empty == "I":
                            # there are non-aligning bases in the source species
                            # between chained alignment blocks before and
                            # after this block.
                            # The browser shows a double line or "=" in base mode.
                            ax.plot([ref.get_forward_strand_start(), ref.get_forward_strand_end()],
                                    [ypos + 0.3, ypos + 0.3], color="black", linewidth=self.properties['line_width'])
                            ax.plot([ref.get_forward_strand_start(), ref.get_forward_strand_end()],
                                    [ypos + 0.7, ypos + 0.7], color="black", linewidth=self.properties['line_width'])
                        elif c.synteny_empty == "M":
                            # there are non-aligning bases in the source and
                            # more than 90% of them are Ns in the source.
                            # The browser shows a pale yellow bar.
                            ax.add_patch(Rectangle((ref.get_forward_strand_start(), ypos),
                                                   ref.size, 1,
                                                   edgecolor="none",
                                                   facecolor="lightyellow"))
                        elif c.synteny_empty == "n":
                            # there are non-aligning bases in the source
                            # and the next aligning block starts
                            # in a new chromosome or scaffold
                            # that is bridged by a chain between
                            # still other blocks.
                            # The browser shows either a single line
                            # or a double line based on how many bases
                            # are in the gap between the bridging alignments.
                            # LD: My observation is that nothing is plotted
                            continue
                        else:
                            self.log.warning(f"Unknown synteny empty code: {c.synteny_empty}"
                                             f" for {assembly}. Nothing is plotted.")
        if valid_blocks == 0:
            self.log.warning("No valid blocks were found in file "
                             f"{self.properties['file']} for region"
                             f"{chrom_region}:{start_region}-{end_region}.\n")
        ymax -= epsilon
        ymin = self.row_scale * current_max_y + epsilon
        # I need to know how many species are plotted before plotting the sequence:
        if self.properties['display_ref_seq']:
            plotting_figure_width = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).width

            # The first constrain on the fontsize is the width
            ideal_fontsize = 1.4 * get_optimal_fontsize(plotting_figure_width,
                                                        start_region, end_region)
            # The other constraint is the height
            real_height_in = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).height
            # We cannot have the sequence that take more than the height of an alignment:
            max_height_in = real_height_in / (current_max_y + 1)

            # 1 point = 1/72 inch = height of character
            max_fontsize = max_height_in * 72

            # Let's take the biggest font possible with these constraints so that the figure is as readable as possible
            fontsize = min(ideal_fontsize, max_fontsize)

            # Evaluate the dedicated space for letters:
            target_height_in = fontsize / 72
            relative_height_letters = target_height_in / real_height_in

            height_seq = (ymin - ymax) / (1 / relative_height_letters - 1)
            ymax -= height_seq
            for i in range(start_region, end_region):
                if i in ref_seq_dic:
                    ax.text(i + 0.5, - height_seq / 2, ref_seq_dic[i],
                            color=seq_color[ref_seq_dic[i]], verticalalignment='center',
                            horizontalalignment='center', fontsize=fontsize)

        self.log.debug(f"ylim {ymin},{ymax}")
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

    def plot_y_axis(self, ax, plot_axis):
        if self.current_labels is not None:
            if int(mpl.__version__.split(".")[1]) < 3:
                wrap = False
            else:
                wrap = True
            for y, label in enumerate(self.current_labels):
                ax.text(0, y * self.row_scale + 0.5, label,
                        verticalalignment='center',
                        horizontalalignment='right', wrap=wrap)
            ax.set_ylim(*plot_axis.get_ylim())

    @staticmethod
    def compare_letters(lr, la):
        """
            Compare a reference letter lr with another letter lr
            We assume lr cannot be '-' while la can.
        """
        if lr == la:
            return "identical"
        if la == '-':
            return "gap"
        return "mismatch"

    # This is inspired from galaxy tools util maf_utilities
    # Except that parse_e_rows=True
    def build_maf_index_species_chromosomes(self, filename, index_species=None):
        species = []
        species_chromosomes = {}
        indexes = bx.interval_index_file.Indexes()
        blocks = 0
        try:
            maf_reader = bx.align.maf.Reader(open(filename), parse_e_rows=True)
            while True:
                pos = maf_reader.file.tell()
                block = next(maf_reader)
                if block is None:
                    break
                blocks += 1
                for c in block.components:
                    spec = c.src
                    chrom = None
                    if "." in spec:
                        spec, chrom = spec.split(".", 1)
                    if spec not in species:
                        species.append(spec)
                        species_chromosomes[spec] = []
                    if chrom and chrom not in species_chromosomes[spec]:
                        species_chromosomes[spec].append(chrom)
                    if (index_species is None or spec in index_species) and \
                       c.size > 0 and not c.empty:
                        # We only index species if the size is positive
                        # And if the it is not a e line because multiple
                        # blocks may correspond to the same region in the src
                        forward_strand_start = c.forward_strand_start
                        forward_strand_end = c.forward_strand_end
                        try:
                            forward_strand_start = int(forward_strand_start)
                            forward_strand_end = int(forward_strand_end)
                        except ValueError:
                            continue
                            # I am not sure this can happen
                            # Original comment from @blankenberg:
                            # start and end are not integers, can't add component to index, goto next component
                            # this likely only occurs when parse_e_rows is True?
                            # could a species exist as only e rows? should the
                        indexes.add(c.src, forward_strand_start, forward_strand_end, pos, max=c.src_size)
        except Exception as e:
            # most likely a bad MAF
            self.log.warning(f'Building MAF index on {filename} failed: {e}')
            return (None, [], {}, 0)
        return (indexes, species, species_chromosomes, blocks)

    def write_maf_index(self):
        indexes, *_ = self.build_maf_index_species_chromosomes(self.properties['file'],
                                                               self.ref)
        if indexes is not None:
            with open(self.properties['file_index'], 'wb') as index:
                indexes.write(index)
        else:
            raise InputError("Unable to generate index for " + self.properties['file'])

    def maf_index_has_good_ref(self):
        ref_species = [ref.split(".")[0] for ref in self.idx.indexes.indexes.keys()]
        return self.ref in ref_species

    def ref_chrom_in_maf_index(self, chrom_name):
        compatible_ref = [ref for ref in self.idx.indexes.indexes.keys()
                          if ref in [f"{self.ref}.{chrom_name}", f"{self.ref}.{change_chrom_names(chrom_name)}"]]
        return compatible_ref
