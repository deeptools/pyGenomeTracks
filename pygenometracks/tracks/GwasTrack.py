import numpy as np
from intervaltree import Interval, IntervalTree
from tqdm import tqdm

from ..readGwas import ReadGwas
from ..utilities import change_chrom_names, transform
from .GenomeTrack import GenomeTrack

DEFAULT_GWAS_COLOR = '#ff7f00'


class GwasTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.gwas', '.linear', '.logistic', '.assoc', '.qassoc']  # this is used by make_tracks_file to guess the type of track based on file name
    TRACK_TYPE = 'gwas'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# File containing the data.
# We expect either:
# - a tabular without header, with the first four columns:
#   CHR, BP, SNP and P.
#   Optionally, extra annotation columns can be added.
# - a tabular with a header indicating the column values.
#   Required columns are:
#   - chromosome (can also be labelled chr or chrom)
#   - position (can also be bp, pos or base_pair_location)
#   - pvalue (can also be p, pval, p-value, p_value, p.value)
#   The header can start with '#'
file =
# Indicate if your file has a header:
file_has_header = false
# Each SNP will be plotted as a 'o' and you can control color/size etc...
# Inside color
#color = red
# Border color
#border_color = black
# Line width
#line_width = 0.5
# Size
#marker_size = 45
# To log transform your PVALUE you can use transform and log_pseudocount:
# For the transform values:
# 'no': do not transform the values
# 'log1p': transformed_values = log(1 + initial_values)
# 'log': transformed_values = log(log_pseudocount + initial_values)
# 'log2': transformed_values = log2(log_pseudocount + initial_values)
# 'log10': transformed_values = log10(log_pseudocount + initial_values)
# '-log': transformed_values = - log(log_pseudocount + initial_values)
# '-log10': transformed_values = - log10(log_pseudocount + initial_values)
# The default is:
# transform = -log10
# log_pseudocount = 0
# When a transformation is applied, by default the y axis
# gives the transformed values, if you prefer to see
# the original values:
#y_axis_values = original
# If you want to have a grid on the y-axis
#grid = true
# set show_data_range to false to hide the text on the left showing the data range
show_data_range = true
# the default for min_value and max_value is 'auto' which means that the scale will go
# roughly from the minimum value found in the region plotted to the maximum value found.
# To change set min_value and max_value before transformation.
# Use for example:
min_value = 1
max_value = 1e-15
# If your gwas file is large and you are plotting large regions
# This can lead to very large pdf/svg files.
# A way to decrease the size of your file
# is to rasterize the dots by using:
# rasterize = true
# Optional. If not given is guessed from the file ending.
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = {'max_value': None,
                           'min_value': None,
                           'show_data_range': True,
                           'orientation': None,
                           'color': DEFAULT_GWAS_COLOR,
                           'border_color': 'black',
                           'line_width': 0.5,
                           'marker_size': 45,
                           'transform': '-log10',
                           'log_pseudocount': 0,
                           'y_axis_values': 'transformed',
                           'file_has_header': False,
                           'rasterize': False,
                           'grid': False}

    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'transform': ['no', 'log', 'log1p', '-log', 'log2',
                                         'log10', '-log10'],
                           'y_axis_values': ['original', 'transformed']}
    BOOLEAN_PROPERTIES = ['file_has_header', 'show_data_range',
                          'rasterize', 'grid']
    STRING_PROPERTIES = ['title', 'file_type', 'file', 'overlay_previous',
                         'orientation', 'color',
                         'border_color', 'transform', 'y_axis_values']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'log_pseudocount': [- np.inf, np.inf],
                        'height': [0, np.inf],
                        'marker_size': [0, np.inf],
                        'line_width': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        super(GwasTrack, self).__init__(*args, **kwarg)
        self.interval_tree = self.process_gwas(self.properties['region'])

    def set_properties_defaults(self):
        super(GwasTrack, self).set_properties_defaults()
        self.process_color('color')
        self.process_color('border_color')

    def process_gwas(self, plot_regions=None):
        """Read the gwas file and store values in a IntervalTree

        :param list plot_regions: list of plotted regions (like [(chrom1, start1, end1), (chrom2, start2, end2)]), defaults to None
        :return None
        """

        gwas_file_h = ReadGwas(self.properties['file'],
                               has_header=self.properties['file_has_header'])

        valid_intervals = 0
        interval_tree = {}

        if plot_regions is not None:
            chroms_to_plot = set([v[0] for v in plot_regions])
        else:
            chroms_to_plot = None

        for record in tqdm(gwas_file_h, total=gwas_file_h.length):

            if plot_regions is not None and record.chromosome not in chroms_to_plot:
                continue

            if record.chromosome not in interval_tree:
                interval_tree[record.chromosome] = IntervalTree()

            interval_tree[record.chromosome].add(Interval(record.position,
                                                          record.position + 1, record))
            valid_intervals += 1

        try:
            gwas_file_h.file_handle.close()
        except AttributeError:
            pass

        if valid_intervals == 0:
            self.log.warning("No valid intervals were found in file "
                             f"{self.properties['file']} for regions"
                             f"{plot_regions}.\n")

        return interval_tree

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plot a scatter plot for the GWAS data.

        :param ax: matplotlib axis
        :param chrom_region: chromosome name
        :param start_region: start position of the region
        :param end_region: end position of the region
        :return: None
        """
        if chrom_region not in self.interval_tree.keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.interval_tree.keys():
                self.log.warning("*Warning*\nNo interval was found when "
                                 "overlapping with both "
                                 f"{chrom_region_before}:{start_region}-{end_region}"
                                 f" and {chrom_region}:{start_region}-{end_region}"
                                 " inside the gwas file. "
                                 "This will generate an empty track!!\n")
                self.adjust_ylim(ax)
                return

        gwas_overlap = \
            self.interval_tree[chrom_region][start_region:end_region]

        # Fill in the position and pvalues lists with data from the GWAS file
        position = [region.begin for region in gwas_overlap]
        score_list = [region.data.pvalue
                      for region in gwas_overlap]
        
        transformed_scores = transform(np.array(score_list),
                                       self.properties['transform'],
                                       self.properties['log_pseudocount'],
                                       self.properties['file'])

        # Plot the scatterplot
        ax.scatter(position, transformed_scores,
                   s=self.properties['marker_size'],
                   color=self.properties['color'], marker='o',
                   edgecolors=self.properties['border_color'],
                   linewidths=self.properties['line_width'],
                   rasterized=self.properties['rasterize'])

        if self.properties['grid']:
            ax.grid(axis='y', zorder=0)

        self.adjust_ylim(ax)

    def plot_y_axis(self, ax, plot_axis):
        super(GwasTrack, self).plot_y_axis(ax, plot_axis,
                                           self.properties['transform'],
                                           self.properties['log_pseudocount'],
                                           self.properties['y_axis_values'],
                                           self.properties['grid'])
