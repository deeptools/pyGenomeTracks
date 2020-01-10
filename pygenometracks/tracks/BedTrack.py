from . GenomeTrack import GenomeTrack
from .. readBed import ReadBed
from .. readGtf import ReadGtf
from .. utilities import opener
import matplotlib
from matplotlib import font_manager
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from intervaltree import IntervalTree, Interval
import numpy as np

DEFAULT_BED_COLOR = '#1f78b4'
DISPLAY_BED_VALID = ['collapsed', 'triangles', 'interleaved', 'stacked']
DISPLAY_BED_SYNONYMOUS = {'interlaced': 'interleaved', 'domain': 'interleaved'}
DEFAULT_DISPLAY_BED = 'stacked'


class BedTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['bed', 'bed3', 'bed4', 'bed5', 'bed6', 'bed8',
                         'bed9', 'bed12',
                         'bed.gz', 'bed3.gz', 'bed4.gz', 'bed5.gz', 'bed6.gz',
                         'bed9.gz', 'bed12.gz',
                         'gtf', 'gtf.gz']
    TRACK_TYPE = 'bed'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# If the bed file contains the exon
# structure (bed 12) then this is plotted. Otherwise
# a region **with direction** is plotted.
# If a gtf is given, it should end with gtf or gtf.gz or
# the type should be defined as gtf:
# type = gtf
# In the case of a gtf file, by default the transcript_name is used.
# If you want to use the gene_name:
# prefered_name = gene_name
# By default, the gtf is transformed to transcripts
# If you want to use see only one structure per gene
# merge_transcripts = true
# If the bed file contains a column for color (column 9), then this color can be used by
# setting:
# color = bed_rgb
#if color is a valid colormap name (like RbBlGn), then the score is mapped
# to the colormap.
# In this case, the the min_value and max_value for the score can be provided, otherwise
# the maximum score and minimum score found are used.
#color = RdYlBu
#min_value=0
#max_value=100
# If the color is simply a color name, then this color is used and the score is not considered.
color = darkblue
# height of track in cm
height = 5
# whether printing the labels
labels = false
# optional:
# by default the labels are not printed if you have more than 60 features.
# to change it, just increase the value:
#max_labels = 60
# optional: font size can be given to override the default size
fontsize = 10
# optional: line_width
#line_width = 0.5
# the display parameter defines how the bed file is plotted.
# The options are ['collapsed', 'interleaved', 'triangles'] These options asume that the regions do not overlap.
# `collapsed`: The bed regions are plotted one after the other in one line.
# `interleaved`: The bed regions are plotted in two lines, first up, then down, then up etc.
# if display is not given, then each region is plotted using the gene style
#optional, default is black. To remove the background color, simply set 'color' and 'background color' to the
# same value
#border_color = black
# style to plot the genes when they have exon information
#style = UCSC
#style = flybase
# maximum number of gene rows to be plotted. This
# field is useful to limit large number of close genes
# to be printed over many rows. When several images want
# to be combined this must be set to get equal size, otherwise, on each image the height of each gene changes
#gene_rows = 10
# by default the ymax is the number of
# rows occupied by the genes in the region plotted. However,
# by setting this option, the global maximum is used instead.
# This is useful to combine images that are all consistent and
# have the same number of rows.
#global_max_row = true
# if you use UCSC style, you can set the relative distance between 2 arrows on introns
# default is 2
#arrow_interval = 2
# if you use flybase style, you can choose the color of non-coding intervals:
#color_utr = grey
# as well as the proportion between their height and the one of coding
# (by default they are the same height):
#height_utr = 1
# By default, for oriented intervals the arrowhead is added
# outside of the interval.
# If you want that the tip of the arrow correspond to
# the extremity of the interval use:
# arrowhead_included = true
# if you want to plot the track on top of the previous track. Options are 'yes' or 'share-y'. For the 'share-y'
# option the y axis values is shared between this plot and the overlay plot. Otherwise, each plot use its own scale
#overlay_previous = yes
# optional. If not given is guessed from the file ending.
file_type = {}
    """.format(TRACK_TYPE)

    DEFAULTS_PROPERTIES = {'fontsize': 12,
                           'orientation': None,
                           'color': DEFAULT_BED_COLOR,
                           'border_color': 'black',
                           'labels': True,
                           'style': 'flybase',
                           'display': DEFAULT_DISPLAY_BED,
                           'interval_height': 100,  # This one is not defined in the documentation
                           'line_width': 0.5,
                           'max_labels': 60,
                           'prefered_name': 'transcript_name',
                           'merge_transcripts': False,
                           'global_max_row': False,
                           'gene_rows': None,
                           'max_value': None,
                           'min_value': None,
                           'arrow_interval': 2,
                           'arrowhead_included': False,
                           'color_utr': 'grey',
                           'height_utr': 1}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None},
                             'display': DISPLAY_BED_SYNONYMOUS}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           'style': ['flybase', 'UCSC'],
                           'display': DISPLAY_BED_VALID}
    BOOLEAN_PROPERTIES = ['labels', 'merge_transcripts', 'global_max_row',
                          'arrowhead_included']
    STRING_PROPERTIES = ['prefered_name', 'file', 'file_type',
                         'overlay_previous', 'orientation',
                         'title', 'style', 'color', 'border_color',
                         'color_utr', 'display']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'fontsize': [0, np.inf],
                        'interval_height': [0, np.inf],
                        'line_width': [0, np.inf],
                        'height': [0, np.inf],
                        'height_utr': [0, 1]}
    INTEGER_PROPERTIES = {'gene_rows': [0, np.inf],
                          'max_labels': [0, np.inf],
                          'arrow_interval': [1, np.inf]}
    # The color can be a color or a colormap or 'bed_rgb'
    # border_color, color_utr can only be a color

    def __init__(self, *args, **kwarg):
        super(BedTrack, self).__init__(*args, **kwarg)
        self.bed_type = None  # once the bed file is read,
        # this is bed3, bed4, bed5, bed6, bed8, bed9 or bed12
        self.len_w = None  # this is the length of the letter 'w' given the font size
        self.interval_tree = {}  # interval tree of the bed regions
        self.interval_tree, min_score, max_score = self.process_bed()
        if self.colormap is not None:
            if self.properties['min_value'] is not None:
                min_score = self.properties['min_value']
            if self.properties['max_value'] is not None:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def set_properties_defaults(self):
        super(BedTrack, self).set_properties_defaults()
        self.fp = font_manager.FontProperties(size=self.properties['fontsize'])
        self.colormap = None

        # check if the color given is a color map
        if not matplotlib.colors.is_color_like(self.properties['color']) \
           and self.properties['color'] != 'bed_rgb':
            # check if the color is a valid colormap name
            if self.properties['color'] not in matplotlib.cm.datad:
                self.log.warning("*WARNING* color: '{}' for section {}"
                                 " is not valid. Color has "
                                 "been set to "
                                 "{}".format(self.properties['color'],
                                             self.properties['section_name'],
                                             DEFAULT_BED_COLOR))
                self.properties['color'] = DEFAULT_BED_COLOR
            else:
                self.colormap = self.properties['color']

        # check if border_color and color_utr are colors
        # if they are part of self.properties
        # (for example, TADsTracks do not have color_utr)
        for param in [p for p in ['border_color', 'color_utr']
                      if p in self.properties]:
            if not matplotlib.colors.is_color_like(self.properties[param]):
                self.log.warning("*WARNING* {}: '{}' for section {}"
                                 " is not valid. Color has "
                                 "been set to "
                                 "{}".format(param,
                                             self.properties[param],
                                             self.properties['section_name'],
                                             self.DEFAULTS_PROPERTIES[param]))
                self.properties[param] = self.DEFAULTS_PROPERTIES[param]

        # to set the distance between rows
        self.row_scale = self.properties['interval_height'] * 2.3

    def get_length_w(self, fig_width, region_start, region_end):
        """
        to improve the visualization of the genes
        it is good to have an estimation of the label
        length. In the following code I try to get the
        length of a 'W' in base pairs.
        """
        if self.properties['labels']:
            # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
            inches_per_pt = 1.0 / 72.27
            font_in_inches = self.properties['fontsize'] * inches_per_pt
            region_len = region_end - region_start
            bp_per_inch = region_len / fig_width
            font_in_bp = font_in_inches * bp_per_inch
            self.len_w = font_in_bp
        else:
            self.len_w = 1

        return self.len_w

    def process_bed(self):

        if self.properties['file'].endswith('gtf') or \
           self.properties['file'].endswith('gtf.gz'):
            bed_file_h = ReadGtf(self.properties['file'],
                                 self.properties['prefered_name'],
                                 self.properties['merge_transcripts'])
        else:
            bed_file_h = ReadBed(opener(self.properties['file']))
        self.bed_type = bed_file_h.file_type

        if self.properties['color'] == 'bed_rgb' and \
           self.bed_type not in ['bed12', 'bed9']:
            self.log.warning("*WARNING* Color set to 'bed_rgb', "
                             "but bed file does not have the rgb field. "
                             "The color has been set to {}".format(DEFAULT_BED_COLOR))
            self.properties['color'] = DEFAULT_BED_COLOR

        valid_intervals = 0
        interval_tree = {}

        max_score = float('-inf')
        min_score = float('inf')
        for bed in bed_file_h:
            if bed.score < min_score:
                min_score = bed.score
            if bed.score > max_score:
                max_score = bed.score

            if bed.chromosome not in interval_tree:
                interval_tree[bed.chromosome] = IntervalTree()

            interval_tree[bed.chromosome].add(Interval(bed.start,
                                                       bed.end, bed))
            valid_intervals += 1

        if valid_intervals == 0:
            self.log.warning("No valid intervals were found in file "
                             "{}".format(self.properties['file_name']))

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

        self.log.debug("max number of rows set to {}".format(self.max_num_row))
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
        # ypos simply oscilates between 0 and 100
        if self.properties['display'] == 'interleaved':
            ypos = self.properties['interval_height'] \
                if self.counter % 2 == 0 \
                else 1
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
            possible_chrom_names = self.get_alternative_chrom_names(chrom_region)
            compatible_chrom_names = [c for c in possible_chrom_names
                                      if c in self.interval_tree]
            if len(compatible_chrom_names) == 0:
                self.log.warning("*Warning*\nNeither " + chrom_region_before
                                 + " nor " + str(possible_chrom_names)
                                 + " exists as a "
                                 "chromosome name inside the bed file. "
                                 "This will generate an empty track!!\n")
                return
            else:
                chrom_region = compatible_chrom_names[0]
        chrom_region = self.check_chrom_str_bytes(self.interval_tree,
                                                  chrom_region)

        genes_overlap = \
            sorted(self.interval_tree[chrom_region][start_region:end_region])

        if self.properties['display'] == 'triangles':
            self.plot_triangles(ax, genes_overlap)
        else:
            self.counter = 0
            self.small_relative = 0.004 * (end_region - start_region)
            self.get_length_w(ax.get_figure().get_figwidth(), start_region,
                              end_region)
            if self.properties['global_max_row']:
                self.get_max_num_row(self.len_w, self.small_relative)

            # do not print labels when too many intervals are visible.
            if self.properties['labels'] and \
               len(genes_overlap) > self.properties['max_labels']:
                self.properties['labels'] = False

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

                if self.properties['labels']:
                    num_name_characters = len(bed.name) + 2
                    # +2 to account for a space before and after the name
                    bed_extended_end = int(bed.end + (num_name_characters * self.len_w))
                else:
                    bed_extended_end = (bed.end + 2 * self.small_relative)

                # get smallest free row
                if len(row_last_position) == 0:
                    free_row = 0
                    row_last_position.append(bed_extended_end)
                else:
                    # get list of rows that are less than bed.start, then take the min
                    idx_list = [idx for idx, value in enumerate(row_last_position)
                                if value < bed.start]
                    if len(idx_list):
                        free_row = min(idx_list)
                        row_last_position[free_row] = bed_extended_end
                    else:
                        free_row = len(row_last_position)
                        row_last_position.append(bed_extended_end)

                rgb, edgecolor = self.get_rgb_and_edge_color(bed)

                ypos = self.get_y_pos(free_row)

                # do not plot if the maximum interval rows to plot is reached
                if self.properties['gene_rows'] is not None and \
                   free_row >= self.properties['gene_rows']:
                    continue

                if free_row > max_num_row_local:
                    max_num_row_local = free_row
                if ypos > max_ypos:
                    max_ypos = ypos

                if self.bed_type == 'bed12':
                    if self.properties['style'] == 'flybase':
                        self.draw_gene_with_introns_flybase_style(ax, bed, ypos,
                                                                  rgb, edgecolor,
                                                                  linewidth)
                    else:
                        self.draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor,
                                                    linewidth)
                else:
                    self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor, linewidth)

                if not self.properties['labels']:
                    pass
                elif bed.end > start_region and bed.end < end_region:
                    ax.text(bed.end + self.small_relative,
                            ypos + (self.properties['interval_height'] / 2),
                            bed.name, horizontalalignment='left',
                            verticalalignment='center', fontproperties=self.fp)

            if self.counter == 0:
                self.log.warning("*Warning* No intervals were found for file {} "
                                 "in section '{}' for the interval plotted"
                                 " ({}:{}-{}).\n".
                                 format(self.properties['file'],
                                        self.properties['section_name'],
                                        chrom_region, start_region, end_region))
            ymax = 0

            if self.properties['global_max_row']:
                ymin = self.max_num_row[chrom_region] * self.row_scale

            elif self.properties['gene_rows'] is not None:
                ymin = self.properties['gene_rows'] * self.row_scale
            else:
                ymin = max_ypos + self.properties['interval_height']

            self.log.debug("ylim {},{}".format(ymin, ymax))
            # the axis is inverted (thus, ymax < ymin)
            ax.set_ylim(ymin, ymax)

            if self.properties['display'] == 'interleaved':
                ax.set_ylim(-5, 205)
            elif self.properties['display'] == 'collapsed':
                ax.set_ylim(-5, 105)

    def plot_label(self, label_ax):
        label_ax.text(0.05, 1, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='top',
                      transform=label_ax.transAxes,
                      wrap=True)

    def plot_y_axis(self, ax, plot_axis):
        if self.colormap is not None:
            self.colormap.set_array([])

            cobar = plt.colorbar(self.colormap, ax=ax, fraction=1,
                                 orientation='vertical')

            cobar.solids.set_edgecolor("face")
            cobar.ax.tick_params(labelsize='smaller')
            cobar.ax.yaxis.set_ticks_position('left')
            # adjust the labels of the colorbar
            ticks = cobar.ax.get_yticks()
            labels = cobar.ax.set_yticklabels(ticks.astype('float32'))
            (vmin, vmax) = cobar.mappable.get_clim()
            for idx in np.where(ticks == vmin)[0]:
                # if the label is at the start of the colobar
                # move it above avoid being cut or overlapping with other track
                labels[idx].set_verticalalignment('bottom')
            for idx in np.where(ticks == vmax)[0]:
                # if the label is at the end of the colobar
                # move it a bit inside to avoid overlapping
                # with other labels
                labels[idx].set_verticalalignment('top')

    def get_rgb_and_edge_color(self, bed):
        rgb = self.properties['color']
        edgecolor = self.properties['border_color']

        if self.colormap:
            # translate value field (in the example above is 0 or 0.2686...)
            # into a color
            rgb = self.colormap.to_rgba(bed.score)
            self.rgb = rgb
        if self.properties['color'] == 'bed_rgb':
            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if self.bed_type in ['bed9', 'bed12'] and len(bed.rgb) == 3:
                try:
                    rgb = [float(x) / 255 for x in bed.rgb]
                except IndexError:
                    rgb = DEFAULT_BED_COLOR
            else:
                rgb = DEFAULT_BED_COLOR
        return rgb, edgecolor

    def draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor, linewidth):
        """
        draws an interval with direction (if given)
        """

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos),
                         bed.end - bed.start,
                         self.properties['interval_height'],
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
        draws a gene using different styles
        """
        if bed.block_count == 0 and bed.thick_start == bed.start and \
           bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = self.properties['interval_height'] / 2
        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height],
                'black', linewidth=linewidth, zorder=-1)

        # get start, end of all the blocks
        positions = []
        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if bed.thick_start == bed.thick_end:
                positions.append((x0, x1, 'UTR'))

            elif x0 < bed.thick_start < x1:
                positions.append((x0, bed.thick_start, 'UTR'))
                positions.append((bed.thick_start, x1, 'coding'))

            elif x0 < bed.thick_end < x1:
                positions.append((x0, bed.thick_end, 'coding'))
                positions.append((bed.thick_end, x1, 'UTR'))

            else:
                if x1 < bed.thick_start or x0 > bed.thick_end:
                    type = 'UTR'
                else:
                    type = 'coding'

                positions.append((x0, x1, type))

        # plot all blocks as rectangles except the last if the strand is + or
        # the first is the strand is -, which are drawn as arrows.
        if bed.strand == '-':
            positions = positions[::-1]

        first_pos = positions.pop()
        if first_pos[2] == 'UTR':
            _rgb = self.properties['color_utr']
            # The arrow will be centered on
            # ypos + self.properties['interval_height'] /2
            # The total height will be
            # self.properties['interval_height'] * self.properties['height_utr']
            y0 = ypos + self.properties['interval_height'] * \
                (1 - self.properties['height_utr']) / 2
            half_height = self.properties['interval_height'] * \
                self.properties['height_utr'] / 2
        else:
            _rgb = rgb
            y0 = ypos
            half_height = self.properties['interval_height'] / 2

        vertices = self._draw_arrow(first_pos[0], first_pos[1], bed.strand,
                                    y0, half_height)

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=_rgb,
                             linewidth=linewidth))

        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = self.properties['color_utr']
                y0 = ypos + self.properties['interval_height'] * \
                    (1 - self.properties['height_utr']) / 2
                height = self.properties['interval_height'] * \
                    self.properties['height_utr']
            else:
                _rgb = rgb
                y0 = ypos
                height = self.properties['interval_height']

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
            half_height = self.properties['interval_height'] / 2
        # The y values are common to both strands:
        y0 = ypos
        y1 = ypos + 2 * half_height
        if strand == '+':
            x0 = start
            if self.properties['arrowhead_included']:
                x1 = max(start, end - self.small_relative)
                x2 = end
            else:
                x1 = end
                x2 = end + self.small_relative
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
                x0 = min(end, start + self.small_relative)
                xb = start
            else:
                x0 = start
                xb = start - self.small_relative
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

    def draw_gene_with_introns(self, ax, bed, ypos, rgb, edgecolor, linewidth):
        """
        draws a gene like in flybase gbrowse.
        """

        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor, linewidth)
            return
        half_height = self.properties['interval_height'] / 2
        quarter_height = self.properties['interval_height'] / 4
        three_quarter_height = quarter_height * 3

        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=linewidth, zorder=-1)

        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x1 < bed.thick_start or x0 > bed.thick_end or \
               bed.thick_start == bed.thick_end:
                y0 = ypos + quarter_height
                y1 = ypos + three_quarter_height
            else:
                y0 = ypos
                y1 = ypos + self.properties['interval_height']

            if x0 < bed.thick_start < x1:
                vertices = ([(x0, ypos + quarter_height), (x0, ypos + three_quarter_height),
                             (bed.thick_start, ypos + three_quarter_height),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (x1, ypos + self.properties['interval_height']), (x1, ypos),
                             (bed.thick_start, ypos), (bed.thick_start, ypos + quarter_height)])

            elif x0 < bed.thick_end < x1:
                vertices = ([(x0, ypos),
                             (x0, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + three_quarter_height),
                             (x1, ypos + three_quarter_height),
                             (x1, ypos + quarter_height),
                             (bed.thick_end, ypos + quarter_height),
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
                if intron_length > self.small_relative:
                    intron_center = x1 + int(intron_length) / 2
                    pos = np.arange(x1 + 1 * self.small_relative,
                                    x1 + intron_length + self.small_relative,
                                    int(arrow_interval * self.small_relative))
                    # center them
                    pos = pos + intron_center - pos.mean()
                    # plot them
                    for xpos in pos:
                        self._plot_small_arrow(ax, xpos, ypos, bed.strand)

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

            rgb, edgecolor = self.get_rgb_and_edge_color(region.data)

            triangle = Polygon([[x1, y1], [x2, y2], [x3, y1]], closed=True,
                               facecolor=rgb, edgecolor=edgecolor, linewidth=self.properties['line_width'])
            ax.add_artist(triangle)
            valid_regions += 1

            if y2 > ymax:
                ymax = y2

        if valid_regions == 0:
            self.log.warning("No regions found for section {}.".format(self.properties['section_name']))

        if self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, 0)
        else:
            ax.set_ylim(0, ymax)

    def _plot_small_arrow(self, ax, xpos, ypos, strand):
        """
        Draws a broken line with 2 parts:
        For strand = +:  > For strand = -: <
        :param xpos:
        :param ypos:
        :param strand:
        :
        :return: None
        """
        if strand == '+':
            xdata = [xpos - self.small_relative / 4,
                     xpos + self.small_relative / 4,
                     xpos - self.small_relative / 4]
        else:
            xdata = [xpos + self.small_relative / 4,
                     xpos - self.small_relative / 4,
                     xpos + self.small_relative / 4]
        ydata = [ypos + self.properties['interval_height'] / 4,
                 ypos + self.properties['interval_height'] / 2,
                 ypos + self.properties['interval_height'] * 3 / 4]
        ax.add_line(Line2D(xdata, ydata, color='black', linewidth=self.properties['line_width']))
