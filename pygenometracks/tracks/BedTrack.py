from . GenomeTrack import GenomeTrack
from .. readBed import ReadBed
# To remove next 1.0
from .. readGtf import ReadGtf
# End to remove
from .. utilities import opener, get_length_w, count_lines, temp_file_from_intersect, change_chrom_names
import matplotlib
from matplotlib import font_manager
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
from intervaltree import IntervalTree, Interval
import numpy as np
from tqdm import tqdm

DEFAULT_BED_COLOR = '#1f78b4'
DISPLAY_BED_VALID = ['collapsed', 'triangles', 'interleaved', 'stacked', 'squares']
DISPLAY_BED_SYNONYMOUS = {'interlaced': 'interleaved', 'domain': 'interleaved'}
DEFAULT_DISPLAY_BED = 'stacked'
AROUND_REGION = 100000


class BedTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['bed', 'bed3', 'bed4', 'bed5', 'bed6', 'bed8',
                         'bed9', 'bed12',
                         'bed.gz', 'bed3.gz', 'bed4.gz', 'bed5.gz', 'bed6.gz',
                         'bed9.gz', 'bed12.gz']
    TRACK_TYPE = 'bed'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# If the bed file contains a column for color (column 9), then this color can be used by
# setting:
#color = bed_rgb
# if color is a valid colormap name (like RbBlGn), then the score (column 5) is mapped
# to the colormap.
# In this case, the the min_value and max_value for the score can be provided, otherwise
# the maximum score and minimum score found are used.
#color = RdYlBu
#min_value=0
#max_value=100
# If the color is simply a color name, then this color is used and the score is not considered.
color = darkblue
# optional: line_width
#line_width = 0.5
# optional: border_color
# default is black.
# To remove the border, simply set 'border_color' to none
# Not used in tssarrow style
#border_color = black
# the display parameter defines how the bed file is plotted.
# Default is 'stacked' where regions are plotted on different lines so
# we can see all regions and all labels.
# The other options are ['collapsed', 'interleaved', 'triangles', 'squares']
# These 2 options assume that the regions do not overlap.
# `collapsed`: The bed regions are plotted one after the other in one line.
# `interleaved`: The bed regions are plotted in two lines, first up, then down, then up etc.
# If the bed file contains the exon
# structure (bed 12) then this is plotted. Otherwise
# a region **with direction** is plotted.
# style to plot the genes when the display is 'stacked', 'collapsed' or 'interleaved'
#style = UCSC
#style = flybase
#style = tssarrow
# maximum number of gene rows to be plotted. This
# field is useful to limit large number of close genes
# to be printed over many rows. When several images want
# to be combined this must be set to get equal size
# otherwise, on each image the height of each gene changes
#gene_rows = 10
# by default the ymax is the number of
# rows occupied by the genes in the region plotted. However,
# by setting this option, the global maximum is used instead.
# This is useful to combine images that are all consistent and
# have the same number of rows.
#global_max_row = true
# whether printing the labels
labels = false
# optional:
# by default the labels are not printed if you have more than 60 features.
# to change it, just increase the value:
#max_labels = 60
# optional: font size can be given to override the default size
fontsize = 10
# If you want to plot all labels inside the plotting region:
#all_labels_inside = true
# If you want to display the name of the gene which goes over the plotted
# region in the right margin put:
#labels_in_margin = true
# If you want to use italic for your labels:
#fontstyle = italic
# if you use UCSC style, you can set the relative distance between 2 arrows on introns
# default is 2
#arrow_interval = 2
# if you use tssarrow style, you can choose the length of the arrow in bp
# (default is 4% of the plotted region)
#arrow_length = 5000
# if you use flybase or tssarrow style, you can choose the color of non-coding intervals:
#color_utr = grey
# as well as the proportion between their height and the one of coding
# (by default they are the same height):
#height_utr = 1
# if you use flybase or UCSC style, you can choose the color of the backbone
#color_backbone = red
# By default, for oriented intervals in flybase style,
# or bed files with less than 12 columns, the arrowhead is added
# outside of the interval.
# If you want that the tip of the arrow correspond to
# the extremity of the interval use:
#arrowhead_included = true
# By default the size of this arrow is 0.4% of the plotted region.
# This size is also used to put space between the bed regions and
# their labels.
# To increase it:
#arrowhead_fraction = 0.01
# The two other display options are really different and no label can be display:
# `triangles` display each region as a triangle, can be useful to overlay with a hic_matrix
# `squares` display each region as a square along the diagonal, can be useful to overlay with a hic_matrix_square
# optional. If not given is guessed from the file ending.
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = {'fontsize': 12,
                           'orientation': None,
                           'color': DEFAULT_BED_COLOR,
                           'border_color': 'black',
                           'labels': True,
                           'style': 'flybase',
                           'display': DEFAULT_DISPLAY_BED,
                           'line_width': 0.5,
                           'max_labels': 60,
                           # To remove in next 1.0
                           'prefered_name': 'transcript_name',
                           'merge_transcripts': False,
                           'merge_overlapping_exons': False,
                           # end to remove
                           'global_max_row': False,
                           'gene_rows': None,
                           'max_value': None,
                           'min_value': None,
                           'arrow_interval': 2,
                           'arrowhead_included': False,
                           'arrowhead_fraction': 0.004,
                           'color_utr': 'grey',
                           'color_backbone': 'black',
                           'height_utr': 1,
                           'region': None,  # Cannot be set manually but is set by tracksClass
                           'arrow_length': None,
                           'all_labels_inside': False,
                           'labels_in_margin': False,
                           'fontstyle': 'normal'}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None},
                             'display': DISPLAY_BED_SYNONYMOUS}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           'style': ['flybase', 'UCSC', 'tssarrow'],
                           'display': DISPLAY_BED_VALID,
                           'fontstyle': ['normal', 'italic', 'oblique']}
    BOOLEAN_PROPERTIES = ['labels', 'global_max_row',
                          'arrowhead_included', 'all_labels_inside',
                          'labels_in_margin',
                          # To remove in next 1.0
                          'merge_transcripts', 'merge_overlapping_exons']
    STRING_PROPERTIES = ['file', 'file_type',
                         'overlay_previous', 'orientation',
                         'title', 'style', 'color', 'border_color',
                         'color_utr', 'display', 'fontstyle',
                         'color_backbone',
                         # To remove in next 1.0
                         'prefered_name']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'fontsize': [0, np.inf],
                        'line_width': [0, np.inf],
                        'height': [0, np.inf],
                        'height_utr': [0, 1],
                        'arrowhead_fraction': [0, np.inf]}
    INTEGER_PROPERTIES = {'gene_rows': [0, np.inf],
                          'max_labels': [0, np.inf],
                          'arrow_interval': [1, np.inf],
                          'arrow_length': [0, np.inf]}

    def __init__(self, *args, **kwarg):
        super(BedTrack, self).__init__(*args, **kwarg)
        self.bed_type = None  # once the bed file is processed,
        # this is bed3, bed4, bed5, bed6, bed8, bed9 or bed12
        self.current_len_w = None  # this is the length of the letter 'w' given the font size
        self.interval_tree = {}  # interval tree of the bed regions
        self.interval_tree, min_score, max_score = self.process_bed(self.properties['region'])
        if self.colormap is not None:
            if self.properties['min_value'] is not None:
                min_score = self.properties['min_value']
            if self.properties['max_value'] is not None:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.colormap)
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def set_properties_defaults(self):
        super(BedTrack, self).set_properties_defaults()
        self.fp = font_manager.FontProperties(size=self.properties['fontsize'])
        self.colormap = None
        self.parametersUsingColormap = []
        # check if the color given is a color map
        is_colormap = self.process_color('color', colormap_possible=True,
                                         bed_rgb_possible=True,
                                         default_value_is_colormap=False)
        if is_colormap:
            self.colormap = self.properties['color']
            self.parametersUsingColormap.append('color')

        # check if border_color and color_utr and color_backbone are colors
        # if they are part of self.properties
        # (for example, TADsTracks do not have color_utr nor color_backbone)
        for param in [p for p in ['border_color', 'color_utr', 'color_backbone']
                      if p in self.properties]:
            is_colormap = self.process_color(param, colormap_possible=True,
                                             bed_rgb_possible=True)
            if is_colormap:
                if self.colormap is None:
                    self.colormap = self.properties[param]
                    self.parametersUsingColormap.append(param)
                else:
                    if self.colormap == self.properties[param]:
                        self.parametersUsingColormap.append(param)
                    else:
                        self.log.warning("*WARNING* section "
                                         f"{self.properties['section_name']}: "
                                         f"{param} was set to "
                                         f"{self.properties[param]}, but "
                                         f"{self.parametersUsingColormap[0]}"
                                         f" was set to {self.colormap}. "
                                         "It is not possible to have multiple"
                                         f" colormap. {param} set to "
                                         f"{self.DEFAULTS_PROPERTIES[param]}.\n")
                        self.properties[param] = self.DEFAULTS_PROPERTIES[param]

        # to set the distance between rows
        self.row_scale = 2.3

    def get_bed_handler(self, plot_regions=None):
        if not self.properties['global_max_row']:
            # I do the intersection:
            file_to_open = temp_file_from_intersect(self.properties['file'],
                                                    plot_regions, AROUND_REGION)
        else:
            file_to_open = self.properties['file']
        # To remove in next 1.0
        if self.properties['file'].endswith('gtf') or \
           self.properties['file'].endswith('gtf.gz'):
            self.log.warning("Deprecation Warning: "
                             f"In section {self.properties['section_name']},"
                             f" file_type was set to {self.TRACK_TYPE}"
                             " whereas it is a gtf file. In the future"
                             " only bed files will be accepted, please"
                             " use file_type = gtf.\n")
            bed_file_h = ReadGtf(file_to_open,
                                 self.properties['prefered_name'],
                                 self.properties['merge_transcripts'],
                                 self.properties['merge_overlapping_exons'])
            total_length = bed_file_h.length
        else:
            # end of remove
            total_length = count_lines(opener(file_to_open),
                                       asBed=True)
            bed_file_h = ReadBed(opener(file_to_open))

        return bed_file_h, total_length

    def process_bed(self, plot_regions=None):

        bed_file_h, total_length = self.get_bed_handler(plot_regions)
        self.bed_type = bed_file_h.file_type

        if self.properties['color'] == 'bed_rgb' and \
           self.bed_type not in ['bed12', 'bed9']:
            self.log.warning("*WARNING* Color set to 'bed_rgb', "
                             "but bed file does not have the rgb field. "
                             f"The color has been set to {DEFAULT_BED_COLOR}.\n")
            self.properties['color'] = DEFAULT_BED_COLOR

        valid_intervals = 0
        interval_tree = {}

        max_score = float('-inf')
        min_score = float('inf')
        for bed in tqdm(bed_file_h, total=total_length):
            if bed.score < min_score:
                min_score = bed.score
            if bed.score > max_score:
                max_score = bed.score

            if bed.chromosome not in interval_tree:
                interval_tree[bed.chromosome] = IntervalTree()

            interval_tree[bed.chromosome].add(Interval(bed.start,
                                                       bed.end, bed))
            valid_intervals += 1

        try:
            bed_file_h.file_handle.close()
        except AttributeError:
            pass

        if valid_intervals == 0:
            self.log.warning("No valid intervals were found in file "
                             f"{self.properties['file']}.\n")

        return interval_tree, min_score, max_score

    def get_max_num_row(self, len_w, small_relative):
        ''' Process the whole bed regions at the given figure length
        and font size to
        determine the maximum number of rows required.
        :return:
        '''

        self.max_num_row = {}
        for chrom in self.interval_tree:
            row_last_position = []  # each entry in this list contains the end position
            self.max_num_row[chrom] = 0
            for region in sorted(self.interval_tree[chrom][0:500000000]):
                bed = region.data
                if self.properties['labels']:
                    bed_extended_end = int(bed.end + (len(bed.name) * len_w))
                    # To uniformize the label position and max_row calc should be:
                    # bed_extended_end = int(bed.end + small_relative + ((len(bed.name) + 2) * len_w))
                else:
                    bed_extended_end = (bed.end + 2 * small_relative)

                # get smallest free row
                if len(row_last_position) == 0:
                    free_row = 0
                    row_last_position.append(bed_extended_end)
                else:
                    # get list of rows that are less than bed.start, then take the min
                    idx_list = [idx for idx, value in enumerate(row_last_position) if value < bed.start]
                    if len(idx_list):
                        free_row = min(idx_list)
                        row_last_position[free_row] = bed_extended_end
                    else:
                        free_row = len(row_last_position)
                        row_last_position.append(bed_extended_end)

                if free_row > self.max_num_row[bed.chromosome]:
                    self.max_num_row[bed.chromosome] = free_row

        self.log.debug(f"max number of rows set to {self.max_num_row}")
        return self.max_num_row

    def get_y_pos(self, free_row):
        """
        The y_pos is set such that regions to be plotted
        do not overlap (stacked). To override this
        the properties['collapsed'] needs to be set.

        The algorithm uses a interval tree (self.region_interval)
        to check the overlaps
        and a sort of coverage vector 'rows used'
        to identify the row in which to plot
        :return: int y position
        """

        # if the interleaved directive is given,
        # ypos simply oscilates between 0 and 1
        if self.properties['display'] == 'interleaved':
            ypos = 1 \
                if self.counter % 2 == 0 \
                else 0
        # if the collapsed directive is given
        # ypos is always 0
        elif self.properties['display'] == 'collapsed':
            ypos = 0
        # if it is stacked
        # it will got the the free_row
        else:
            ypos = free_row * self.row_scale
        return ypos

    def plot(self, ax, chrom_region, start_region, end_region):
        if chrom_region not in self.interval_tree.keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.interval_tree.keys():
                self.log.warning("*Warning*\nNo interval was found when "
                                 "overlapping with both "
                                 f"{chrom_region_before}:{start_region - AROUND_REGION}-{end_region + AROUND_REGION}"
                                 f" and {chrom_region}:{start_region - AROUND_REGION}-{end_region + AROUND_REGION}"
                                 " inside the bed file. "
                                 "This will generate an empty track!!\n")
                return

        genes_overlap = \
            sorted(self.interval_tree[chrom_region][start_region:end_region])

        if self.properties['display'] == 'triangles':
            self.plot_triangles(ax, genes_overlap)
        elif self.properties['display'] == 'squares':
            self.plot_squares(ax, genes_overlap)
            ax.set_ylim(end_region, start_region)
        else:
            self.counter = 0
            self.current_small_relative = self.properties['arrowhead_fraction'] * (end_region - start_region)
            if self.properties['labels']:
                self.current_len_w = get_length_w(ax.get_figure().get_figwidth(),
                                                  start_region, end_region,
                                                  self.properties['fontsize'])
            else:
                self.current_len_w = 1

            if self.properties['global_max_row']:
                self.get_max_num_row(self.current_len_w, self.current_small_relative)

            # do not print labels when too many intervals are visible.
            display_labels = self.properties['labels']
            if self.properties['labels'] and \
               len(genes_overlap) > self.properties['max_labels']:
                display_labels = False

            linewidth = self.properties['line_width']
            max_num_row_local = 1
            max_ypos = 0
            # check for the number of other intervals that overlap
            #    with the given interval
            #            1         2
            #  012345678901234567890123456
            #  1=========       4=========
            #       2=========
            #         3============
            #
            # for 1 row_last_position = [9]
            # for 2 row_last_position = [9, 14]
            # for 3 row_last_position = [9, 14, 19]
            # for 4 row_last_position = [26, 14, 19]

            row_last_position = []
            # each entry in this list contains the end position
            # of genomic interval. The list index is the row
            # in which the genomic interval was plotted.
            # Any new genomic interval that wants to be plotted,
            # knows the row to use by finding the list index that
            # is larger than its start

            # check for overlapping genes including
            # label size (if plotted)

            if ax.get_xlim()[0] > ax.get_xlim()[1]:
                genes_overlap = reversed(genes_overlap)
            for region in genes_overlap:
                """
                BED12 gene format with exon locations at the end
                chrX    20850   23076   CG17636-RA      0       -       20850   23017   0       3       946,765,64,     0,1031,2162,

                BED9
                bed with rgb at end
                chr2L   0       70000   ID_5    0.26864549832   .       0       70000   51,160,44

                BED6
                bed without rgb
                chr2L   0       70000   ID_5    0.26864549832   .
                """
                self.counter += 1
                bed = region.data

                if ax.get_xlim()[0] < ax.get_xlim()[1]:
                    bed_left = bed.start
                    bed_right = bed.end

                    def add_to_right(a, b):
                        return a + b

                    def add_to_left(a, b):
                        return a - b

                    def is_left_to(a, b):
                        return a < b

                    def is_right_to(a, b):
                        return a > b

                else:
                    bed_left = bed.end
                    bed_right = bed.start

                    def add_to_right(a, b):
                        return a - b

                    def add_to_left(a, b):
                        return a + b

                    def is_left_to(a, b):
                        return a > b

                    def is_right_to(a, b):
                        return a < b

                if display_labels:
                    num_name_characters = len(bed.name) + 2
                    # +2 to account for a space before and after the name
                    bed_extended_right = int(add_to_right(bed_right, (num_name_characters * self.current_len_w)))
                    # To uniformize the label position and max_row calc should be:
                    # bed_extended_right = int(add_to_right(bed_right, (num_name_characters * self.current_len_w + self.current_small_relative)))
                else:
                    bed_extended_right = add_to_right(bed_right, 2 * self.current_small_relative)

                bed_extended_left = bed_left
                # get smallest free row
                if len(row_last_position) == 0:
                    free_row = 0
                    row_last_position.append(bed_extended_right)
                else:
                    # If all_labels_inside = True
                    # genes which goes over will have their labels inside
                    if self.properties['all_labels_inside'] and display_labels \
                       and is_right_to(bed_extended_right, ax.get_xlim()[1]):
                        bed_extended_left = int(add_to_left(bed_left, (num_name_characters * self.current_len_w)))
                        # To uniformize the label position and max_row calc should be:
                        # bed_extended_left = int(add_to_left(bed_left, (num_name_characters * self.current_len_w + self.current_small_relative)))

                        # Check that the start position is not outside:
                        if is_left_to(bed_extended_left, ax.get_xlim()[0]):
                            # If it would be outside, we use the default right label
                            bed_extended_left = bed_left
                        else:
                            # If we keep the label to the left, we update the right extended
                            bed_extended_right = add_to_right(bed_right, 2 * self.current_small_relative)

                    # get list of rows that are left to bed_extended_left, then take the min
                    idx_list = [idx for idx, value in enumerate(row_last_position)
                                if is_left_to(value, bed_extended_left)]
                    if len(idx_list):
                        free_row = min(idx_list)
                        row_last_position[free_row] = bed_extended_right
                    else:
                        free_row = len(row_last_position)
                        row_last_position.append(bed_extended_right)

                rgb = self.get_rgb(bed)
                edgecolor = self.get_rgb(bed, param='border_color', default=rgb)

                ypos = self.get_y_pos(free_row)

                # do not plot if the maximum interval rows to plot is reached
                if self.properties['gene_rows'] is not None and \
                   free_row >= self.properties['gene_rows']:
                    continue

                if free_row > max_num_row_local:
                    max_num_row_local = free_row
                if ypos > max_ypos:
                    max_ypos = ypos

                if self.properties['style'] == 'tssarrow':
                    self.draw_gene_tssarrow_style(ax, bed, ypos, rgb,
                                                  linewidth)
                elif self.bed_type == 'bed12':
                    if self.properties['style'] == 'flybase':
                        self.draw_gene_with_introns_flybase_style(ax, bed, ypos,
                                                                  rgb, edgecolor,
                                                                  linewidth)
                    else:
                        self.draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor,
                                                    linewidth)
                else:
                    self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor, linewidth)

                if not display_labels:
                    pass
                elif bed_extended_left != bed_left:
                    # The label will be plotted before
                    ax.text(add_to_left(bed_left, self.current_small_relative),
                            ypos + (1 / 2),
                            bed.name, horizontalalignment='right',
                            verticalalignment='center', fontproperties=self.fp,
                            fontstyle=self.properties['fontstyle'])
                    # To uniformize the label position and max_row calc should be:
                    # ax.text(add_to_left(bed_left, self.current_small_relative + self.current_len_w),
                elif bed_right > start_region and bed_right < end_region:
                    ax.text(add_to_right(bed_right, self.current_small_relative),
                            ypos + 0.5,
                            bed.name, horizontalalignment='left',
                            verticalalignment='center', fontproperties=self.fp,
                            fontstyle=self.properties['fontstyle'])
                    # To uniformize the label position and max_row calc should be:
                    # ax.text(add_to_right(bed_right, self.current_small_relative + self.current_len_w),
                elif self.properties['labels_in_margin'] \
                        and (bed_right == end_region or is_right_to(bed_right, end_region)):
                    ax.text(add_to_right(ax.get_xlim()[1], self.current_small_relative),
                            ypos + (1 / 2),
                            bed.name, horizontalalignment='left',
                            verticalalignment='center', fontproperties=self.fp,
                            fontstyle=self.properties['fontstyle'])
                    # To uniformize the label position and max_row calc should be:
                    # ax.text(add_to_right(ax.get_xlim()[1], self.current_small_relative + self.current_len_w),

            if self.counter == 0:
                self.log.warning("*Warning* No intervals were found for file"
                                 f" {self.properties['file']} in "
                                 f"section '{self.properties['section_name']}'"
                                 " for the interval plotted"
                                 f" ({chrom_region}:{start_region}-{end_region}).\n")

            epsilon = 0.08
            ymax = - epsilon

            # We set ymin and ymax to have genes centered epsilon from the border

            if self.properties['global_max_row']:
                max_ypos = self.max_num_row[chrom_region] * self.row_scale

            elif self.properties['gene_rows'] is not None:
                max_ypos = self.properties['gene_rows'] * self.row_scale

            ymin = max_ypos + (1 + epsilon)

            self.log.debug(f"ylim {ymin},{ymax}")
            # the axis is inverted (thus, ymax < ymin)
            ax.set_ylim(ymin, ymax)

            if self.properties['display'] == 'interleaved':
                ax.set_ylim(2 + epsilon, ymax)
            elif self.properties['display'] == 'collapsed':
                ax.set_ylim(1 + epsilon, ymax)

        if self.properties['orientation'] == 'inverted':
            ylims = ax.get_ylim()
            ax.set_ylim(ylims[1], ylims[0])

    def plot_label(self, label_ax, width_dpi, h_align='left'):
        if h_align == 'left':
            label_ax.text(0.05, 1, self.properties['title'],
                          horizontalalignment='left', size='large',
                          verticalalignment='top',
                          transform=label_ax.transAxes,
                          wrap=True)
        elif h_align == 'right':
            txt = label_ax.text(1, 1, self.properties['title'],
                                horizontalalignment='right', size='large',
                                verticalalignment='top',
                                transform=label_ax.transAxes,
                                wrap=True)
            # To be able to wrap to the left:
            txt._get_wrap_line_width = lambda: width_dpi
        else:
            txt = label_ax.text(0.5, 1, self.properties['title'],
                                horizontalalignment='center', size='large',
                                verticalalignment='top',
                                transform=label_ax.transAxes,
                                wrap=True)
            # To be able to wrap to the left:
            txt._get_wrap_line_width = lambda: width_dpi

    def plot_y_axis(self, ax, plot_axis):
        if self.colormap is not None:
            self.colormap.set_array([])
            GenomeTrack.plot_custom_cobar(self, ax, fraction=1)

    def get_rgb(self, bed, param='color', default=DEFAULT_BED_COLOR):
        """
        get the rgb value for the bed and the param given:
        :param bed:
        :param param:
        :param default: the default value if it fails
        :return: color
        """
        rgb = self.properties[param]

        if self.colormap is not None and param in self.parametersUsingColormap:
            # translate value field (in the example above is 0 or 0.2686...)
            # into a color
            rgb = self.colormap.to_rgba(bed.score)
        elif self.properties[param] == 'bed_rgb':
            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if self.bed_type in ['bed9', 'bed12'] and len(bed.rgb) == 3:
                try:
                    rgb = [float(x) / 255 for x in bed.rgb]
                except IndexError:
                    rgb = default
            else:
                rgb = default
        return rgb

    def draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor, linewidth):
        """
        draws an interval with direction (if given)
        """

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos),
                         bed.end - bed.start, 1,
                         edgecolor=edgecolor, facecolor=rgb,
                         linewidth=linewidth))
        else:
            vertices = self._draw_arrow(bed.start, bed.end, bed.strand, ypos)
            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=rgb,
                                 linewidth=linewidth))

    def draw_gene_with_introns_flybase_style(self, ax, bed, ypos, rgb,
                                             edgecolor, linewidth):
        """
        draws a gene like in flybase gbrowse.
        """
        if bed.block_count == 0 and bed.thick_start == bed.start and \
           bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor, linewidth)
            return
        half_height = 1 / 2
        # draw 'backbone', a line from the start until the end of the gene
        rgb_backbone = self.get_rgb(bed, param='color_backbone', default='black')
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height],
                color=rgb_backbone, linewidth=linewidth, zorder=-1)

        # get start, end of all the blocks
        positions = self._split_bed_to_blocks(bed)

        if bed.strand != '.':
            # plot all blocks as rectangles except the last if the strand is + or
            # the first is the strand is -, which are drawn as arrows.
            if bed.strand == '-':
                positions = positions[::-1]

            first_pos = positions.pop()
            if first_pos[2] == 'UTR':
                _rgb = self.get_rgb(bed, param='color_utr', default=rgb)
                # The arrow will be centered on
                # ypos + 1 / 2
                # The total height will be
                # self.properties['height_utr']
                y0 = ypos + (1 - self.properties['height_utr']) / 2
                half_height = self.properties['height_utr'] / 2
            else:
                _rgb = rgb
                y0 = ypos
                half_height = 1 / 2

            vertices = self._draw_arrow(first_pos[0], first_pos[1], bed.strand,
                                        y0, half_height)

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=linewidth))

        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = self.get_rgb(bed, param='color_utr', default=rgb)
                y0 = ypos + (1 - self.properties['height_utr']) / 2
                height = self.properties['height_utr']
            else:
                _rgb = rgb
                y0 = ypos
                height = 1

            vertices = [(start_pos, y0), (start_pos, y0 + height),
                        (end_pos, y0 + height), (end_pos, y0)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=linewidth))

    def _draw_arrow(self, start, end, strand, ypos, half_height=None):
        """
        Draws a filled arrow.
        :param start:
        :param end:
        :param strand:
        :param ypos:
        :param half_height:
        :return: None
        """
        if half_height is None:
            half_height = 1 / 2
        # The y values are common to both strands:
        y0 = ypos
        y1 = ypos + 2 * half_height
        if strand == '+':
            x0 = start
            if self.properties['arrowhead_included']:
                x1 = max(start, end - self.current_small_relative)
                x2 = end
            else:
                x1 = end
                x2 = end + self.current_small_relative
            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            =================>
            x0             x1 x2

            """

            vertices = [(x0, y0), (x0, y1), (x1, y1),
                        (x2, y0 + half_height), (x1, y0)]

        else:
            if self.properties['arrowhead_included']:
                x0 = min(end, start + self.current_small_relative)
                xb = start
            else:
                x0 = start
                xb = start - self.current_small_relative
            x1 = end
            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

              <=================
            xb x0              x1
            """
            vertices = [(x0, y0), (xb, y0 + half_height), (x0, y1),
                        (x1, y1), (x1, y0)]

        return vertices

    def _split_bed_to_blocks(self, bed):
        """
        Split a bed entry into blocks to plot
        :param bed: a namedtuple with at least 6 fields
        :return: a list of tuple (start, end, type) with type in ['UTR', 'coding']
        """
        if self.bed_type != 'bed12':
            # No thick_start, thick_end, block_count:
            return [(bed.start, bed.end, 'coding')]

        # get start, end of all the blocks
        positions = []
        for idx in range(0, bed.block_count):
            # x0 and x1 are the start/end of the current block
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            # We deal with the special case where
            # there is no coding independently
            if bed.thick_start == bed.thick_end:
                positions.append((x0, x1, 'UTR'))
                continue
            # If the beginning of the coding region
            # is withing the current block
            if x0 < bed.thick_start < x1:
                # What is before is UTR
                positions.append((x0, bed.thick_start, 'UTR'))
                # The start of the interval is updated
                x0 = bed.thick_start

            # If the end of the coding region
            # is withing the current block
            if x0 < bed.thick_end < x1:
                # What is before is coding
                positions.append((x0, bed.thick_end, 'coding'))
                # The start of the interval is updated
                x0 = bed.thick_end

            if x1 < bed.thick_start or x0 >= bed.thick_end:
                type = 'UTR'
            else:
                type = 'coding'

            positions.append((x0, x1, type))

        return positions

    def draw_gene_with_introns(self, ax, bed, ypos, rgb, edgecolor, linewidth):
        """
        draws a gene like in UCSC
        Except that for the moment no arrow are plotted
        on the coding part
        """

        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor, linewidth)
            return

        # draw 'backbone', a line from the start until the end of the gene
        rgb_backbone = self.get_rgb(bed, param='color_backbone', default='black')
        ax.plot([bed.start, bed.end], [ypos + 1 / 2, ypos + 1 / 2], color=rgb_backbone, linewidth=linewidth, zorder=-1)

        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x1 < bed.thick_start or x0 > bed.thick_end or \
               bed.thick_start == bed.thick_end:
                y0 = ypos + 1 / 4
                y1 = ypos + 3 / 4
            else:
                y0 = ypos
                y1 = ypos + 1

            if x0 < bed.thick_start < x1 and x0 < bed.thick_end < x1:
                vertices = ([(x0, ypos + 1 / 4),
                             (x0, ypos + 3 / 4),
                             (bed.thick_start, ypos + 3 / 4),
                             (bed.thick_start, ypos + 1),
                             (bed.thick_end, ypos + 1),
                             (bed.thick_end, ypos + 3 / 4),
                             (x1, ypos + 3 / 4),
                             (x1, ypos + 1 / 4),
                             (bed.thick_end, ypos + 1 / 4),
                             (bed.thick_end, ypos),
                             (bed.thick_start, ypos),
                             (bed.thick_start, ypos + 1 / 4)])
            elif x0 < bed.thick_start < x1:
                vertices = ([(x0, ypos + 1 / 4),
                             (x0, ypos + 3 / 4),
                             (bed.thick_start, ypos + 3 / 4),
                             (bed.thick_start, ypos + 1),
                             (x1, ypos + 1),
                             (x1, ypos),
                             (bed.thick_start, ypos),
                             (bed.thick_start, ypos + 1 / 4)])
            elif x0 < bed.thick_end < x1:
                vertices = ([(x0, ypos),
                             (x0, ypos + 1),
                             (bed.thick_end, ypos + 1),
                             (bed.thick_end, ypos + 3 / 4),
                             (x1, ypos + 3 / 4),
                             (x1, ypos + 1 / 4),
                             (bed.thick_end, ypos + 1 / 4),
                             (bed.thick_end, ypos)])
            else:
                vertices = ([(x0, y0), (x0, y1), (x1, y1), (x1, y0)])

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 linewidth=linewidth,
                                 edgecolor='none',
                                 facecolor=rgb))

            if idx < bed.block_count - 1:
                # plot small arrows over the back bone
                intron_length = bed.block_starts[idx + 1] - (bed.block_starts[idx] + bed.block_sizes[idx])
                arrow_interval = self.properties['arrow_interval']
                if intron_length > self.current_small_relative:
                    intron_center = x1 + int(intron_length) / 2
                    pos = np.arange(x1 + 1 * self.current_small_relative,
                                    x1 + intron_length + self.current_small_relative,
                                    int(arrow_interval * self.current_small_relative))
                    # center them
                    pos = pos + intron_center - pos.mean()
                    # plot them
                    for xpos in pos:
                        self._plot_small_arrow(ax, xpos, ypos, bed.strand, bed)

    def draw_gene_tssarrow_style(self, ax, bed, ypos, rgb, linewidth):
        """
        Draw genes like this:
          -->
          |
          ----------- ^ ---
          |         |   | |
          -----------   ---
        """
        # get start, end of all the blocks
        positions = self._split_bed_to_blocks(bed)

        y_bottom = ypos + 1
        y_top_intron = ypos + 1 / 4

        if bed.strand in ["+", "-"]:
            if self.properties['arrow_length'] is None:
                arrow_length = 10 * self.current_small_relative
            else:
                arrow_length = self.properties['arrow_length']
            y_arrow = ypos + 1 / 8
            head_width = 1 / 4
            head_length = self.current_small_relative * 3
            # plot the arrow to indicate tss
            if bed.strand == "+":
                x = bed.start
                dx = arrow_length
            else:
                x = bed.end
                dx = - arrow_length
            # First plot the vertical line:
            ax.add_line(Line2D((x, x), (y_bottom, y_arrow), color=rgb, linewidth=linewidth))
            # Then the arrow
            ax.arrow(x, y_arrow, dx, 0, overhang=1, width=0,
                     head_width=head_width,
                     head_length=head_length,
                     length_includes_head=True,
                     color=rgb, linewidth=linewidth)

        # plot all blocks as rectangles like in the flybase mode but
        # with half the height and no border
        # as well as the junction between exons:
        last_corner = None
        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = self.get_rgb(bed, param='color_utr', default=rgb)
                y0 = y_bottom - 1 / 2 * (1 - self.properties['height_utr']) / 2
                height = 1 / 2 * self.properties['height_utr']
            else:
                _rgb = rgb
                y0 = y_bottom
                height = 1 / 2

            vertices = [(start_pos, y0), (start_pos, y0 - height),
                        (end_pos, y0 - height), (end_pos, y0)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor="none",
                                 facecolor=_rgb,
                                 linewidth=linewidth))
            if last_corner is not None:
                if last_corner[0] < start_pos:
                    xdata = (last_corner[0], (last_corner[0] + start_pos) / 2,
                             start_pos)
                    ydata = (last_corner[1], y_top_intron,
                             y0 - height)
                    ax.add_line(Line2D(xdata, ydata, color=last_corner[2],
                                       linewidth=linewidth))

            last_corner = (end_pos, y0 - height, _rgb)

    def plot_triangles(self, ax, genes_overlap):
        """
        Plots the boundaries as triangles in the given ax.
        """
        ymax = 0.001
        valid_regions = 0
        for region in genes_overlap:
            """
                   ______ y2
                  ""
                 "  "
                "    "
               "      "_____ y1
            _____________________
               x1 x2 x3

            """
            x1 = region.begin
            x2 = x1 + float(region.end - region.begin) / 2
            x3 = region.end
            y1 = 0
            y2 = (region.end - region.begin)

            rgb = self.get_rgb(region.data)
            edgecolor = self.get_rgb(region.data, param='border_color', default=rgb)

            triangle = Polygon([[x1, y1], [x2, y2], [x3, y1]], closed=True,
                               facecolor=rgb, edgecolor=edgecolor, linewidth=self.properties['line_width'])
            ax.add_artist(triangle)
            valid_regions += 1

            if y2 > ymax:
                ymax = y2

        if valid_regions == 0:
            self.log.warning(f"No regions found for section {self.properties['section_name']}.\n")

        ax.set_ylim(0, ymax)

    def plot_squares(self, ax, genes_overlap):
        """
        Plots the boundaries as squares along the diagonal in the given ax.
        """
        valid_regions = 0
        for region in genes_overlap:
            """
        2->  "  " <- 1
        3->  "  " <- 0
            """
            x0 = region.end
            x2 = region.begin

            y0 = x0
            x1 = x0
            y2 = x2
            y1 = y2
            x3 = x2
            y3 = y0
            rgb = self.get_rgb(region.data)
            edgecolor = self.get_rgb(region.data, param='border_color', default=rgb)

            rectangle = Polygon(np.array([[x0, y0], [x1, y1], [x2, y2], [x3, y3]]),
                                facecolor=rgb, edgecolor=edgecolor,
                                linewidth=self.properties['line_width'])
            ax.add_artist(rectangle)
            valid_regions += 1

        if valid_regions == 0:
            self.log.warning(f"No regions found for section {self.properties['section_name']}.\n")

    def _plot_small_arrow(self, ax, xpos, ypos, strand, bed):
        """
        Draws a broken line with 2 parts:
        For strand = +:  > For strand = -: <
        :param xpos:
        :param ypos:
        :param strand:
        :
        :return: None
        """
        if strand == '.':
            return
        if strand == '+':
            xdata = [xpos - self.current_small_relative / 4,
                     xpos + self.current_small_relative / 4,
                     xpos - self.current_small_relative / 4]
        else:
            xdata = [xpos + self.current_small_relative / 4,
                     xpos - self.current_small_relative / 4,
                     xpos + self.current_small_relative / 4]
        ydata = [ypos + 1 / 4,
                 ypos + 1 / 2,
                 ypos + 3 / 4]

        rgb_backbone = self.get_rgb(bed, param='color_backbone', default='black')
        ax.add_line(Line2D(xdata, ydata, color=rgb_backbone, linewidth=self.properties['line_width']))
