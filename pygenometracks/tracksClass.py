# -*- coding: utf-8 -*-
from __future__ import division
from past.builtins import map
from past.builtins import zip

import numpy as np
import logging
from configparser import ConfigParser

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.textpath
import matplotlib.colors
import matplotlib.gridspec
import matplotlib.cm
import mpl_toolkits.axisartist as axisartist
import textwrap
from . utilities import file_to_intervaltree

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ndarray size changed")
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=ImportWarning)

# import warnings
# warnings.filterwarnings('error')

from collections import OrderedDict
from intervaltree import IntervalTree, Interval

from pygenometracks.tracks import *

FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
logging.basicConfig(format=FORMAT)
log = logging.getLogger(__name__)
# logging.basicConfig()
# log = logging.getLogger("tracksClass")
log.setLevel(logging.DEBUG)

DEFAULT_BEDGRAPH_COLOR = '#a6cee3'
DEFAULT_TRACK_HEIGHT = 0.5  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
# proportion of width dedicated to (figure, legends)
DEFAULT_WIDTH_RATIOS = (0.90, 0.1)
DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0.03, 'top': 0.97}


class MultiDict(OrderedDict):
    """
    Class to allow identically named
    sections in configuration file
    by appending the section number as
    for example:
    1. section name
    """
    _unique = 0

    def __setitem__(self, key, val):
        if isinstance(val, OrderedDict):
            self._unique += 1
            key = "{}. [{}]".format(str(self._unique), key)
        OrderedDict.__setitem__(self, key, val)


class PlotTracks(object):

    def __init__(self, tracks_file, fig_width=DEFAULT_FIGURE_WIDTH,
                 fig_height=None, fontsize=None, dpi=None, track_label_width=None):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.dpi = dpi
        self.vlines_intval_tree = None
        self.vlines_properties = None
        self.track_list = None
        start = self.print_elapsed(None)
        self.available_tracks = self.get_available_tracks()
        self.parse_tracks(tracks_file)
        if fontsize:
            fontsize = fontsize
        else:
            fontsize = float(fig_width) * 0.3
        # the track label width is the fraction of the figure width that is used
        # for the track 'title' or label.
        if track_label_width is None:
            self.width_ratios = DEFAULT_WIDTH_RATIOS
        else:
            self.width_ratios = (1 - track_label_width, track_label_width)

        font = {'size': fontsize}
        matplotlib.rc('font', **font)
        # initialize each track
        self.track_obj_list = []
        for idx, properties in enumerate(self.track_list):
            if 'spacer' in properties:
                self.track_obj_list.append(SpacerTrack(properties))
                continue
            elif 'x-axis' in properties:
                self.track_obj_list.append(XAxisTrack(properties))
                continue
            else:
                # for all other tracks that are not axis or spacer
                # the track_class is obtained from the available tracks
                track_class = self.available_tracks[properties['file_type']]
                self.track_obj_list.append(track_class(properties))

            if 'title' in properties:
                # adjust titles that are too long
                # if the track label space is small
                if track_label_width < 0.1:
                    if sys.version_info[0] == 2:
                        properties['title'] = textwrap.fill(properties['title'].encode("UTF-8"), 12)
                    else:
                        properties['title'] = textwrap.fill(properties['title'], 12)
                else:
                    if sys.version_info[0] == 2:
                        properties['title'] = textwrap.fill(properties['title'].encode("UTF-8"), 30)
                    else:
                        properties['title'] = textwrap.fill(properties['title'], 30)

        log.info("time initializing track(s):")
        self.print_elapsed(start)

    @staticmethod
    def get_available_tracks():
        avail_tracks = {}
        work = [GenomeTrack]
        while work:
            parent = work.pop()
            for child in parent.__subclasses__():
                if child not in avail_tracks:
                    track_type = child.TRACK_TYPE
                    avail_tracks[track_type] = child
                    work.append(child)
        return avail_tracks

    def get_tracks_height(self, start_region=None, end_region=None):
        """
        The main purpose of the following loop is
        to get the height of each of the tracks
        because for the Hi-C the height is variable with respect
        to the range being plotted, the function is called
        when each plot is going to be printed.

        Args:
            start_region: start of the region to plot. Only used in case the plot is a Hi-C matrix
            end_region: end of the region to plot. Only used in case the plot is a Hi-C matrix

        Returns:

        """
        track_height = []
        for track_dict in self.track_list:
            # if overlay previous is set to a value other than no
            # then, skip this track height
            if track_dict['overlay previous'] != 'no':
                continue
            elif 'x-axis' in track_dict and track_dict['x-axis'] is True:
                height = track_dict['fontsize'] / 10
            elif 'height' in track_dict:
                height = track_dict['height']
            # compute the height of a Hi-C track
            # based on the depth such that the
            # resulting plot appears proportional
            #
            #      /|\
            #     / | \
            #    /  |d \   d is the depth that we want to be proportional
            #   /   |   \  when plotted in the figure
            # ------------------
            #   region len
            #
            # d (in cm) =  depth (in bp) * width (in cm) / region len (in bp)

            elif 'depth' in track_dict and track_dict['file_type'] == 'hic_matrix':
                # to compute the actual width of the figure the margins and the region
                # set for the legends have to be considered
                # DEFAULT_MARGINS[1] - DEFAULT_MARGINS[0] is the proportion of plotting area

                hic_width = \
                    self.fig_width * (DEFAULT_MARGINS['right'] - DEFAULT_MARGINS['left']) * self.width_ratios[0]
                scale_factor = 0.6  # the scale factor is to obtain a 'pleasing' result.
                depth = min(track_dict['depth'], (end_region - start_region))

                height = scale_factor * depth * hic_width / (end_region - start_region)
            else:
                height = DEFAULT_TRACK_HEIGHT

            track_height.append(height)

        return track_height

    def plot(self, file_name, chrom, start, end, title=None):
        track_height = self.get_tracks_height(start_region=start, end_region=end)

        if self.fig_height:
            fig_height = self.fig_height
        else:
            fig_height = sum(track_height)

        log.debug("Figure size in cm is {} x {}. Dpi is set to {}\n".format(self.fig_width,
                                                                            fig_height, self.dpi))
        fig = plt.figure(figsize=self.cm2inch(self.fig_width, fig_height))
        if title:
            fig.suptitle(title)

        grids = matplotlib.gridspec.GridSpec(len(track_height), 2,
                                             height_ratios=track_height,
                                             width_ratios=self.width_ratios)
        axis_list = []
        # skipped_tracks is the count of tracks that have the
        # 'overlay previous' parameter and should be skipped
        skipped_tracks = 0
        axis = None
        for idx, track in enumerate(self.track_obj_list):
            if idx == 0 and track.properties['overlay previous'] != 'no':
                log.warn("First track can not have the `overlay previous` option")
                track.properties['overlay previous'] = 'no'

            if track.properties['overlay previous'] in ['yes', 'share-y']:
                overlay = True
                skipped_tracks += 1
            else:
                overlay = False

            if track.properties['overlay previous'] == 'share-y':
                ylim = axis.get_ylim()
            else:
                idx -= skipped_tracks
                axis = axisartist.Subplot(fig, grids[idx, 0])
                fig.add_subplot(axis)
                # turns off the lines around the tracks
                axis.axis[:].set_visible(False)
                # to make the background transparent
                axis.patch.set_visible(False)
                label_axis = plt.subplot(grids[idx, 1])
                label_axis.set_axis_off()

            axis.set_xlim(start, end)
            track.plot(axis, label_axis, chrom, start, end)

            if track.properties['overlay previous'] == 'share-y':
                axis.set_ylim(ylim)

            if not overlay:
                axis_list.append(axis)

        if self.vlines_intval_tree:
            self.plot_vlines(axis_list, chrom, start, end)

        fig.subplots_adjust(wspace=0, hspace=0.0,
                            left=DEFAULT_MARGINS['left'],
                            right=DEFAULT_MARGINS['right'],
                            bottom=DEFAULT_MARGINS['bottom'],
                            top=DEFAULT_MARGINS['top'])

        fig.savefig(file_name, dpi=self.dpi, transparent=False)
        return fig.get_size_inches()

    def plot_vlines(self, axis_list, chrom_region, start_region, end_region):
        """
        Plots dotted lines from the top of the first plot to the bottom
        of the last plot at the specified positions.

        :param axis_list: list of plotted axis
        :param chrom_region chromosome name
        :param start_region start position
        :param end_region end position

        :return: None
        """
        vlines_list = []
        if 'line width' in self.vlines_properties:
            line_width = self.vlines_properties['line width']
        else:
            line_width = 0.5

        if chrom_region not in list(self.vlines_intval_tree):
            chrom_region = GenomeTrack.change_chrom_names(chrom_region)
        chrom_region = GenomeTrack.check_chrom_str_bytes(self.vlines_intval_tree, chrom_region)

        for region in sorted(self.vlines_intval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            vlines_list.append(region.begin)

        for ax in axis_list:
            ymin, ymax = ax.get_ylim()

            ax.vlines(vlines_list, ymin, ymax, linestyle='dashed', zorder=10, linewidth=line_width,
                      color=(0, 0, 0, 0.7), alpha=0.5)

        return

    def parse_tracks(self, tracks_file):
        """
        Parses a configuration file

        :param tracks_file: file path containing the track configuration
        :return: array of dictionaries and vlines_file. One dictionary per track
        """
        from ast import literal_eval
        parser = ConfigParser(dict_type=MultiDict, strict=False)
        parser.read_file(open(tracks_file, 'r'))

        tracks_file_path = os.path.dirname(tracks_file)

        track_list = []
        for section_name in parser.sections():
            track_options = dict({"section_name": section_name})
            if section_name.endswith('[spacer]'):
                track_options['spacer'] = True
            elif section_name.endswith('[x-axis]'):
                track_options['x-axis'] = True
            for name, value in parser.items(section_name):
                if name in ['max_value', 'min_value', 'depth', 'height', 'line width',
                            'fontsize', 'scale factor', 'number of bins'] and value != 'auto':
                    track_options[name] = literal_eval(value)
                else:
                    track_options[name] = value

            if 'type' in track_options and track_options['type'] == 'vlines':
                self.vlines_properties = self.check_file_exists(track_options, tracks_file_path)
            elif 'skip' in track_options and track_options['skip'] != 'no':
                pass
            else:
                track_list.append(track_options)

        updated_track_list = []
        for track_dict in track_list:
            warn = None
            if 'file' in track_dict and track_dict['file'] != '':
                track_dict = self.check_file_exists(track_dict, tracks_file_path)
                if 'file_type' not in track_dict:
                    track_dict['file_type'] = self.guess_filetype(track_dict, self.available_tracks)

            if 'overlay previous' not in track_dict:
                track_dict['overlay previous'] = 'no'
            #  set some default values
            if 'title' not in track_dict:
                track_dict['title'] = ''
                if track_dict['overlay previous'] != 'no' or track_dict['section_name'].endswith('[x-axis]') \
                        or track_dict['section_name'].endswith('[spacer]'):
                    pass
                else:
                    warn = "\ntitle not set for 'section {}'\n".format(track_dict['section_name'])
            if warn:
                sys.stderr.write(warn)
            updated_track_list.append(track_dict)
        self.track_list = updated_track_list
        if self.vlines_properties:
            self.vlines_intval_tree, __, __ = file_to_intervaltree(self.vlines_properties['file'])

    @staticmethod
    def check_file_exists(track_dict, tracks_path):
        """
        Checks if a file or list of files exists. If the file does not exists
        tries to check if the file may be relative to the track_file path, in such case
        the path is updated.
        :param track_dict: dictionary of track values. Should contain
                            a 'file' key containing the path of the file
                            or files to be checked separated by space
                            For example: file1 file2 file3
        :param tracks_path: path of the tracks file
        :return: None
        """
        for file_type in ['boundaries_file', 'file']:
            if file_type in track_dict:
                file_names = [x for x in track_dict[file_type].split(" ") if x != '']
                full_path_file_names = []
                for file_name in file_names:
                    try:
                        open(file_name, 'r').close()
                        full_path_file_names.append(file_name)
                    except IOError:
                        try:
                            # try to find the file in the same path as the
                            # the path of the
                            name_with_tracks_path = tracks_path + "/" + file_name
                            open(name_with_tracks_path, 'r').close()
                            full_path_file_names.append(name_with_tracks_path)
                        except IOError:
                            sys.stderr.write("\n*ERROR*\nFile in section [{}] "
                                             "not found:\n{}\n\n".format(track_dict['section_name'],
                                                                         file_name))
                            sys.exit(1)

                track_dict[file_type] = " ".join(full_path_file_names)
        return track_dict

    @staticmethod
    def guess_filetype(track_dict, available_tracks):
        """

        :param track_dict: dictionary of track values with the 'file' key
                    containing a string path of the file or files. Only the ending
                     of the last file is used in case when there are more files
        :param: available_tracks: list of available tracks

        :return: string file type detected
        """
        file_ = track_dict['file'].strip()
        file_type = None
        for track_type, track_class in available_tracks.items():
            for ending in track_class.SUPPORTED_ENDINGS:
                if file_.endswith(ending):
                    if file_type == track_class.TRACK_TYPE:
                        log.error("file_type already defined in other GenomeTrack")
                        exit()
                    else:
                        file_type = track_class.TRACK_TYPE

        if file_type is None:
            sys.exit("Section {}: can not identify file type. Please specify "
                     "the file_type for '{}'".format(track_dict['section_name'], file_))

        return file_type

    @staticmethod
    def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

    @staticmethod
    def print_elapsed(start):
        import time
        if start:
            log.info(time.time() - start)
        return time.time()


class SpacerTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = None

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        pass


class XAxisTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = None

    def __init__(self, *args, **kwargs):
        super(XAxisTrack, self).__init__(*args, **kwargs)
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 15

    def plot(self, ax, label_axis, chrom_region, region_start, region_end):
        ticks = ax.get_xticks()
        if ticks[-1] - ticks[1] <= 1e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"

        elif 1e5 < ticks[-1] - ticks[1] < 4e6:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"
        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in ticks]
            labels[-2] += " Mbp"

        if 'where' in self.properties and self.properties['where'] == 'top':
            ax.axis["x"] = ax.new_floating_axis(0, 0.2)
            ax.axis["x"].set_axis_direction("top")
            label_y_pos = 0.99
            vert_align = 'top'
        else:
            ax.axis["x"] = ax.new_floating_axis(0, 0.9)
            label_y_pos = 0.01
            vert_align = 'bottom'
        ax.text(0.5, label_y_pos, chrom_region, horizontalalignment='center',
                fontsize=int(self.properties['fontsize']), verticalalignment=vert_align, transform=ax.transAxes)

        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')

        ax.axis["x"].major_ticklabels.set(size=int(self.properties['fontsize']))
        label_axis.text(0.15, 0.5, self.properties['title'],
                        horizontalalignment='left', size='large',
                        verticalalignment='center')


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
                    msg = 'File not valid. The format is chrom1 start1, end1, ' \
                          'chrom2, start2, end2, score\nError: {}\n in line\n {}'.format(detail, line)
                    sys.exit(msg)

                try:
                    start1 = int(start1)
                    end1 = int(end1)
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError as detail:
                    msg = "Error reading line: {}. One of the fields is not " \
                          "an integer.\nError message: {}".format(line_number, detail)
                    sys.exit(msg)

                assert start1 <= end1, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)
                assert start2 <= end2, "Error in line #{}, end2 larger than start2 in {}".format(line_number, line)
                try:
                    score = float(score)
                except ValueError as detail:
                    msg = "Error reading line: {}. The score is not valid {}. " \
                          "\nError message: {}".format(line_number, detail)
                    sys.exit(msg)

                if chrom1 != chrom2:
                    sys.stderr.write("Only links in same chromosome are used. Skipping line\n{}\n".format(line))
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
            sys.stderr.write("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

    def plot(self, ax, label_ax, chrom_region, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between Hi-C bins.
        :param ax: matplotlib axis
        :param label_ax: matplotlib axis for labels
        """
        self.max_height = 0
        count = 0

        if chrom_region not in list(self.interval_tree):
            chrom_region = self.change_chrom_names(chrom_region)
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
        log.debug("{} were links plotted".format(count))
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(self.max_height, -1)
        else:
            ax.set_ylim(-1, self.max_height)

        log.debug('title is {}'.format(self.properties['title']))
        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')

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
