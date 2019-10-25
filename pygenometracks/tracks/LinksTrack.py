from . GenomeTrack import GenomeTrack
from intervaltree import IntervalTree, Interval
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Polygon

DEFAULT_LINKS_COLOR = 'blue'


class LinksTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.arcs', '.arc', '.link', '.links', '.bedpe']
    TRACK_TYPE = 'links'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# the file format for links is (tab separated)
#   chr1 start1 end1 chr2 start2 end2 (score ...)
# The score field is optional
# The fields after the score field will be ignored
# for example:
#   chr1 100 200 chr1 250 300 0.5
# depending on the links type either an 'arc' or a 'triangle' or a 'loop' can be plotted.
# If an arc, a line will be drawn from the center of the first region (chr1: 150),
# to the center of the other region (chr1: 275).
# if a triangle, the vertix of the triangle will be drawn at the center between the two points (also the center of
# each position is used)
# if a loop, a rectangle highlighting the intersection between the 2 regions will be shown
# the triangles, and loops options are convenient to overlay over a
# Hi-C matrix to highlight the matrix pixel of the highlighted link
# For these tracks do not hesitate to put large line width like 5 or 10.
links type = arcs
# color of the lines
# if color is a valid colormap name (like RdYlGn),
# then the score is mapped to the colormap.
color = red
# if line width is not given, the score is used to set the line width
# using the following formula (0.5 * square root(score)
# line width = 0.5
# options for line style are 'solid', 'dashed', 'dotted' etc. The full list of
# styles can be found here: https://matplotlib.org/gallery/lines_bars_and_markers/linestyles.html
line style = solid
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwarg):
        super(LinksTrack, self).__init__(*args, **kwarg)
        # the file format expected is similar to file format of links in
        # circos:
        # chr1 100 200 chr1 250 300 0.5
        # where the last value is a score.
        if 'line style' not in self.properties:
            self.properties['line style'] = 'solid'
        if 'links type' not in self.properties:
            self.properties['links type'] = 'arcs'
        self.max_height = None
        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        has_score = True
        max_score = float('-inf')
        min_score = float('inf')
        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                    continue
                try:
                    chrom1, start1, end1, chrom2, start2, end2 = line.strip().split('\t')[:6]
                except Exception as detail:
                    self.log.error('File not valid. The format is chrom1 start1, end1, '
                                   'chrom2, start2, end2\nError: {}\n in line\n {}'.format(detail, line))
                    exit(1)
                try:
                    score = line.strip().split('\t')[6]
                except IndexError:
                    has_score = False
                    score = np.nan

                try:
                    start1 = int(start1)
                    end1 = int(end1)
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError as detail:
                    self.log.error("Error reading line: {}. One of the fields is not "
                                   "an integer.\nError message: {}".format(line_number, detail))
                    exit(1)

                assert start1 <= end1, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)
                assert start2 <= end2, "Error in line #{}, end2 larger than start2 in {}".format(line_number, line)

                if has_score:
                    try:
                        score = float(score)
                    except ValueError as detail:
                        self.log.warning("Warning: reading line: {}. The score is not valid {} will not be used. "
                                         "\nError message: {}".format(line_number, score, detail))
                        score = np.nan
                        has_score = False
                    else:
                        if score < min_score:
                            min_score = score
                        if score > max_score:
                            max_score = score

                if chrom1 != chrom2:
                    self.log.warn("Only links in same chromosome are used. Skipping line\n{}\n".format(line))
                    continue

                if chrom1 not in interval_tree:
                    interval_tree[chrom1] = IntervalTree()

                if start2 < start1:
                    start1, start2 = start2, start1
                    end1, end2 = end2, end1

                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(Interval(start1, end2, [start1, end1, start2, end2, score]))
                valid_intervals += 1

        if valid_intervals == 0:
            self.log.warn("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

        if 'line width' not in self.properties and not has_score:
            self.log.warning("*WARNING* for section {}"
                             " no line width has been set but some "
                             "lines do not have scores."
                             "line width has been set to "
                             "0.5".format(self.properties['section_name']))
            self.properties['line width'] = 0.5

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_LINKS_COLOR
        self.colormap = None
        # check if the color given is a color map
        if not matplotlib.colors.is_color_like(self.properties['color']):
            # check if the color is a valid colormap name
            if self.properties['color'] not in matplotlib.cm.datad:
                self.log.warning("*WARNING* color: '{}' for section {}"
                                 " is not valid. Color has "
                                 "been set to "
                                 "{}".format(self.properties['color'],
                                             self.properties['section_name'],
                                             DEFAULT_LINKS_COLOR))
                self.properties['color'] = DEFAULT_LINKS_COLOR
            else:
                if not has_score:
                    self.log.warning("*WARNING* for section {}"
                                     " a colormap was chosen but some "
                                     "lines do not have scores."
                                     "Color has been set to "
                                     "{}".format(self.properties['section_name'],
                                                 DEFAULT_LINKS_COLOR))
                    self.properties['color'] = DEFAULT_LINKS_COLOR
                else:
                    self.colormap = self.properties['color']

        if self.colormap is not None:
            if 'min_value' in self.properties:
                min_score = self.properties['min_value']
            if 'max_value' in self.properties:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def plot(self, ax, chrom_region, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between Hi-C bins.
        :param ax: matplotlib axis
        :param label_ax: matplotlib axis for labels
        """
        self.max_height = 0
        count = 0

        if chrom_region not in list(self.interval_tree):
            chrom_region_before = chrom_region
            chrom_region = self.change_chrom_names(chrom_region)
            if chrom_region not in list(self.interval_tree):
                self.log.error("*Error*\nNeither " + chrom_region_before + " "
                               "nor " + chrom_region + " exits as a chromosome"
                               " name inside the link file.\n")
                return

        chrom_region = self.check_chrom_str_bytes(self.interval_tree, chrom_region)

        arcs_in_region = sorted(self.interval_tree[chrom_region][region_start:region_end])

        for idx, interval in enumerate(arcs_in_region):
            # skip intervals whose start and end are outside the plotted region
            if interval.begin < region_start and interval.end > region_end:
                continue

            if 'line width' in self.properties:
                self.line_width = float(self.properties['line width'])
            else:
                self.line_width = 0.5 * np.sqrt(interval.data[4])

            if self.properties['links type'] == 'triangles':
                self.plot_triangles(ax, interval)
            elif self.properties['links type'] == 'loops':
                self.plot_loops(ax, interval.data)
            else:
                self.plot_arcs(ax, interval)

            count += 1

        # the arc height is equal to the radius, the track height is the largest
        # radius plotted plus an small increase to avoid cropping of the arcs
        self.max_height += self.max_height * 0.1
        self.log.debug("{} were links plotted".format(count))
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(self.max_height, -1)
        else:
            ax.set_ylim(-1, self.max_height)

        self.log.debug('title is {}'.format(self.properties['title']))

    def plot_y_axis(self, ax, plot_ax):
        if self.colormap is not None and self.properties['overlay previous'] == 'no':
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

    def plot_arcs(self, ax, interval):

        diameter = (interval.end - interval.begin)
        radius = float(diameter) / 2
        center = interval.begin + float(diameter) / 2
        if radius > self.max_height:
            self.max_height = radius
        ax.plot([center], [diameter])
        if self.colormap:
            # translate score field
            # into a color
            rgb = self.colormap.to_rgba(interval.data[4])
        else:
            rgb = self.properties['color']
        ax.add_patch(Arc((center, 0), diameter,
                         diameter, 0, 0, 180, color=rgb,
                         linewidth=self.line_width, ls=self.properties['line style']))

    def plot_triangles(self, ax, interval):
        x1 = interval.begin
        x2 = x1 + float(interval.end - interval.begin) / 2
        x3 = interval.end
        y1 = 0
        y2 = (interval.end - interval.begin)
        if self.colormap:
            # translate score field
            # into a color
            rgb = self.colormap.to_rgba(interval.data[4])
        else:
            rgb = self.properties['color']

        triangle = Polygon(np.array([[x1, y1], [x2, y2], [x3, y1]]),
                           closed=False,
                           facecolor='none', edgecolor=rgb,
                           linewidth=self.line_width,
                           ls=self.properties['line style'])
        ax.add_artist(triangle)
        if y2 > self.max_height:
            self.max_height = y2

    def plot_loops(self, ax, loop):
        """
              " <- 2
        3->  "  " <- 1
               " <- 0
            """
        width1 = loop[1] - loop[0]
        width2 = loop[3] - loop[2]
        x0 = (loop[1] + loop[2]) / 2
        y0 = loop[2] - loop[1]

        x1 = x0 + width2 / 2
        y1 = y0 + width2

        x2 = (loop[0] + loop[3]) / 2
        y2 = loop[3] - loop[0]

        x3 = x0 - width1 / 2
        y3 = y0 + width1

        if self.colormap:
            # translate score field
            # into a color
            rgb = self.colormap.to_rgba(loop[4])
        else:
            rgb = self.properties['color']

        rectangle = Polygon(np.array([[x0, y0], [x1, y1], [x2, y2], [x3, y3]]),
                            facecolor='none', edgecolor=rgb,
                            linewidth=self.line_width,
                            ls=self.properties['line style'])
        ax.add_artist(rectangle)
        if y2 > self.max_height:
            self.max_height = y2
