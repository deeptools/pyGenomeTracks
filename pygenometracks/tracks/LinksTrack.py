from . GenomeTrack import GenomeTrack
from intervaltree import IntervalTree, Interval
import numpy as np


class LinksTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.arcs', '.arc' '.link', '.links']
    TRACK_TYPE = 'links'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# the file format for links is (tab separated)
#   chr1 start1 end1 chr2 start2 end2 score
# for example:
#   chr1 100 200 chr1 250 300 0.5
# depending on the links type either and arc or a 'triangle' can be plotted. If an arc,
# a line will be drawn from the center of the first region (chr1: 150, tot the center of the other region (chr1:275).
# if a triangle, the vertix of the triangle will be drawn at the center between the two points (also the center of
# each position is used)
# links whose start or end is not in the region plotted are not shown.
# color of the lines
color = red
# for the links type, the options are arcs and triangles, the triangles option is convenient to overlay over a
# Hi-C matrix to highlight the matrix pixel of the highlighted link
links type = arcs
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
        if 'line width' not in self.properties:
            self.properties['line width'] = 0.5
        if 'line style' not in self.properties:
            self.properties['line style'] = 'solid'
        if 'links type' not in self.properties:
            self.properties['links type'] = 'arcs'
        self.max_height = None
        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                    continue
                try:
                    chrom1, start1, end1, chrom2, start2, end2, score = line.strip().split('\t')
                except Exception as detail:
                    self.log.error('File not valid. The format is chrom1 start1, end1, '
                                   'chrom2, start2, end2, score\nError: {}\n in line\n {}'.format(detail, line))
                    exit(1)

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
                try:
                    score = float(score)
                except ValueError as detail:
                    self.log.error("Error reading line: {}. The score is not valid {}. "
                                   "\nError message: {}".format(line_number, detail))
                    exit(1)

                if chrom1 != chrom2:
                    self.log.warn("Only links in same chromosome are used. Skipping line\n{}\n".format(line))
                    continue

                if chrom1 not in interval_tree:
                    interval_tree[chrom1] = IntervalTree()

                if start2 < start1:
                    start1, start2 = start2, start1
                    end1, end2 = end2, end1

                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(Interval(start1, end2, score))
                valid_intervals += 1

        if valid_intervals == 0:
            self.log.warn("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

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
                self.line_width = 0.5 * np.sqrt(interval.data)

            if self.properties['links type'] == 'triangles':
                self.plot_triangles(ax, interval)
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
        pass

    def plot_arcs(self, ax, interval):
        from matplotlib.patches import Arc

        diameter = (interval.end - interval.begin)
        radius = float(diameter) / 2
        center = interval.begin + float(diameter) / 2
        if radius > self.max_height:
            self.max_height = radius
        ax.plot([center], [diameter])
        ax.add_patch(Arc((center, 0), diameter,
                         diameter, 0, 0, 180, color=self.properties['color'],
                         linewidth=self.line_width, ls=self.properties['line style']))

    def plot_triangles(self, ax, interval):
        from matplotlib.patches import Polygon
        x1 = interval.begin
        x2 = x1 + float(interval.end - interval.begin) / 2
        x3 = interval.end
        y1 = 0
        y2 = (interval.end - interval.begin)

        triangle = Polygon(np.array([[x1, y1], [x2, y2], [x3, y1]]), closed=False,
                           facecolor='none', edgecolor=self.properties['color'],
                           linewidth=self.line_width,
                           ls=self.properties['line style'])
        ax.add_artist(triangle)
        if y2 > self.max_height:
            self.max_height = y2
