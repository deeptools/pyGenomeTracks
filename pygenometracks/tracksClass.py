# -*- coding: utf-8 -*-

import logging
import os
from configparser import ConfigParser
import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import matplotlib.textpath
import matplotlib.colors
import matplotlib.gridspec
import matplotlib.cm
import mpl_toolkits.axisartist as axisartist
from . utilities import file_to_intervaltree, change_chrom_names, MyBasePairFormatter, get_region
from collections import OrderedDict
from pygenometracks.tracks.GenomeTrack import GenomeTrack
from pygenometracks.utilities import InputError

import warnings

matplotlib.use('Agg')

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ndarray size changed")
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=ImportWarning)

# import warnings
# warnings.filterwarnings('error')

FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
logging.basicConfig(format=FORMAT)
log = logging.getLogger(__name__)
# logging.basicConfig()
# log = logging.getLogger("tracksClass")
log.setLevel(logging.DEBUG)

DEFAULT_TRACK_HEIGHT = 0.5  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0.03, 'top': 0.97}


class VType(object):
    """The VType object is a holder for all vertical types
    which will be plotted above other tracks.

    It is expected that all VType has a plot method
    """
    DEFAULTS_PROPERTIES = {}
    NECESSARY_PROPERTIES = []
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['type']
    FLOAT_PROPERTIES = {}
    INTEGER_PROPERTIES = {}

    def __init__(self, properties_dict):
        FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
        logging.basicConfig(format=FORMAT)
        log = logging.getLogger(__name__)
        log.setLevel(logging.DEBUG)
        self.log = log
        self.properties = properties_dict
        self.set_properties_defaults()
        self.type = 'test'

    def set_properties_defaults(self):
        # put to default all properties which are not set:
        for prop in self.DEFAULTS_PROPERTIES:
            if prop not in self.properties:
                self.properties[prop] = self.DEFAULTS_PROPERTIES[prop]
        # check if properties are possible:
        for prop in self.POSSIBLE_PROPERTIES:
            possibles = self.POSSIBLE_PROPERTIES[prop]
            if self.properties[prop] not in possibles:
                default_value = self.DEFAULTS_PROPERTIES[prop]
                self.log.warning(f"*WARNING* {prop}: '{self.properties[prop]}'"
                                 " for section "
                                 f"{self.properties['section_name']}"
                                 f" is not valid. {prop} has "
                                 "been set to "
                                 f"{default_value}.\n")
                self.properties[prop] = default_value

    def process_color(self, param):
        """
        Put a valid color in self.properties[param]
        Args:
            param: param to check/update the value

        Returns:
            None

        """
        if not matplotlib.colors.is_color_like(self.properties[param]):
            default_value = self.DEFAULTS_PROPERTIES[param]
            log.warning("*Warning*\nThe color "
                        f"{self.properties[param]} specified"
                        " is not a valid color."
                        f" {default_value} will be used."
                        "\n")
            self.properties[param] = default_value
        return


class VlineType(VType):
    TYPE = 'vlines'
    OPTIONS_TXT = f"""
# Mandatory:
# bed file with positions to plot (only start are used)
file = tads.bed
# type:
type = {TYPE}
# Optional:
# line width
#line_width = 0.5
# color of the line:
#color = black
# transparency:
#alpha = 0.7
# zorder, choose positive value (above) to put on top
# and null or negative values to put behind
#zorder = 10
# line style choose among 'solid', 'dashed', 'dashdot' or 'dotted'
# line_style = 'dashed'
    """
    DEFAULTS_PROPERTIES = {'line_width': 0.5,
                           'color': 'black',
                           'alpha': 0.7,
                           'zorder': 10,
                           'line_style': 'dashed'}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'line_style': ['solid', 'dashed', 'dotted', 'dashdot']
                           }
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['file', 'type',
                         'color', 'line_style']
    FLOAT_PROPERTIES = {'line_width': [0, np.inf],
                        'zorder': [-np.inf, np.inf],
                        'alpha': [0, 1]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        super(VlineType, self).__init__(*args, **kwarg)
        # load the file as an interval tree
        self.interval_tree, __, __ = file_to_intervaltree(self.properties['file'],
                                                          self.properties['region'])

    def set_properties_defaults(self):
        super(VlineType, self).set_properties_defaults()
        self.process_color('color')

    def plot(self, axis_list, figure, chrom_region, start_region, end_region):
        """
        If the zorder is positive:
        Plots lines from the top of the first plot to the bottom
        of the last plot at the specified positions.
        Else:
        Plots vlines on each track

        :param axis_list: list of plotted axis
        :param chrom_region chromosome name
        :param start_region start position
        :param end_region end position

        :return: None
        """
        if chrom_region not in list(self.interval_tree):
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in list(self.interval_tree):
                log.warning("*Warning*\nNo interval was found when "
                            f"overlapping with both {chrom_region_before}:{start_region}-{end_region}"
                            f" and {chrom_region}:{start_region}-{end_region} inside the "
                            "file with vertical lines. "
                            "No vertical lines will be "
                            "plotted!!\n")
                return
        vlines_list = []
        for region in sorted(self.interval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            vlines_list.append(region.begin)

        for pos in vlines_list:
            con = ConnectionPatch(xyA=(pos, axis_list[0].get_ylim()[1]), xyB=(pos, axis_list[-1].get_ylim()[0]),
                                  coordsA=axis_list[0].transData, coordsB=axis_list[-1].transData,
                                  axesA=axis_list[0],
                                  axesB=axis_list[-1], color=self.properties['color'],
                                  lw=self.properties['line_width'], alpha=self.properties['alpha'],
                                  zorder=self.properties['zorder'], ls=self.properties['line_style'])
            figure.add_artist(con)


class VhighlightType(VType):
    TYPE = 'vhighlight'
    OPTIONS_TXT = f"""
# Mandatory:
# bed file with regions to highlight
file = regions.bed
# type:
type = {TYPE}
# Optional:
# line width
#line_width = 0.5
# border_color of the line (put none to have no border):
#border_color = black
# color of the rectangle:
#color = yellow
# transparency:
#alpha = 0.3
# zorder, choose high values (10 or 100) to put on top
# and low values to put behind (0 or negative values)
#zorder = -100
    """
    DEFAULTS_PROPERTIES = {'line_width': 0.5,
                           'border_color': 'none',
                           'color': 'yellow',
                           'alpha': 0.5,
                           'zorder': -100}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['file', 'type',
                         'color', 'border_color']
    FLOAT_PROPERTIES = {'line_width': [0, np.inf],
                        'zorder': [-np.inf, np.inf],
                        'alpha': [0, 1]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        super(VhighlightType, self).__init__(*args, **kwarg)
        # load the file as an interval tree
        self.interval_tree, __, __ = file_to_intervaltree(self.properties['file'],
                                                          self.properties['region'])

    def set_properties_defaults(self):
        super(VhighlightType, self).set_properties_defaults()
        self.process_color('color')
        self.process_color('border_color')

    def plot(self, axis_list, figure, chrom_region, start_region, end_region):
        """
        Plots rectangle on each track,
        at the specified positions.

        :param axis_list: list of plotted axis
        :param chrom_region chromosome name
        :param start_region start position
        :param end_region end position

        :return: None
        """
        if chrom_region not in list(self.interval_tree):
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in list(self.interval_tree):
                log.warning("*Warning*\nNo interval was found when "
                            f"overlapping with both {chrom_region_before}:{start_region}-{end_region}"
                            f" and {chrom_region}:{start_region}-{end_region} inside the "
                            "file with vhighlight. "
                            "No vhighlight will be "
                            "plotted!!\n")
                return

        for ax in axis_list:
            for region in sorted(self.interval_tree[chrom_region][start_region - 10000:end_region + 10000]):
                ax.axvspan(region.begin, region.end,
                           facecolor=self.properties['color'],
                           alpha=self.properties['alpha'],
                           zorder=self.properties['zorder'],
                           lw=self.properties['line_width'],
                           edgecolor=self.properties['border_color'])


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
            key = f"{str(self._unique)}. [{key}]"
        OrderedDict.__setitem__(self, key, val)


class PlotTracks(object):

    def __init__(self, tracks_file, fig_width=DEFAULT_FIGURE_WIDTH,
                 fig_height=None, fontsize=None, dpi=None,
                 track_label_width=0.1,
                 plot_regions=None, plot_width=None):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.dpi = dpi
        self.type_list = None
        self.track_list = None
        start = self.print_elapsed(None)
        self.available_tracks = self.get_available_tracks()
        self.available_types = self.get_available_types()
        self.parse_tracks(tracks_file, plot_regions=plot_regions)
        if fontsize:
            fontsize = fontsize
        else:
            fontsize = float(fig_width) * 0.3
        # the track label width is the fraction of
        # the figure width that is used
        # for the track 'title' or label.
        self.width_ratios = (0.01,
                             1 - track_label_width,
                             track_label_width)

        # Process the width:
        if plot_width is not None:
            self.fig_width = plot_width \
                / (DEFAULT_MARGINS['right'] - DEFAULT_MARGINS['left']) \
                * (1 + 2 / 3 * 0.01) \
                / (self.width_ratios[1] / sum(self.width_ratios))

        font = {'size': fontsize}
        matplotlib.rc('font', **font)
        # initialize each track
        self.track_obj_list = []
        for idx, properties in enumerate(self.track_list):
            log.info(f"initialize {properties['section_name']}")
            track_class = properties['track_class']
            properties['region'] = plot_regions.copy()
            self.track_obj_list.append(track_class(properties))

        # initialize each type
        self.type_obj_list = []
        for idx, properties in enumerate(self.type_list):
            log.info(f"initialize {properties['section_name']}")
            track_class = properties['track_class']
            properties['region'] = plot_regions.copy()
            self.type_obj_list.append(track_class(properties))

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

    @staticmethod
    def get_available_types():
        avail_types = {}
        work = [VType]
        while work:
            parent = work.pop()
            for child in parent.__subclasses__():
                if child not in avail_types:
                    track_type = child.TYPE
                    avail_types[track_type] = child
                    work.append(child)
        return avail_types

    def get_tracks_height(self, start_region=None, end_region=None):
        """
        The main purpose of the following loop is
        to get the height of each of the tracks
        because for the Hi-C the height is variable with respect
        to the range being plotted, the function is called
        when each plot is going to be printed.

        Args:
            start_region: start of the region to plot.
                          Only used in case the plot is a Hi-C matrix
            end_region: end of the region to plot.
                        Only used in case the plot is a Hi-C matrix

        Returns:

        """
        track_height = []
        for i, track_dict in enumerate(self.track_list):
            if i == 0 and track_dict['overlay_previous'] != 'no':
                log.warning("First track can not have the `overlay_previous` option.\n")
                self.track_list[i]['overlay_previous'] = 'no'
            # if overlay_previous is set to a value other than no
            # then, skip this track height
            if track_dict['overlay_previous'] != 'no':
                continue
            elif 'height' in track_dict:
                height = track_dict['height']
            elif track_dict['file_type'] == 'x_axis':
                height = track_dict['fontsize'] / 8
            elif ('depth' in track_dict
                  and track_dict['file_type'] == 'hic_matrix') or \
                    track_dict['file_type'] == 'hic_matrix_square':
                # compute the height of a Hi-C track
                # based on the depth such that the
                # resulting plot appears proportional
                # For hic_matrix:
                #
                #      /|\
                #     / | \
                #    /  |d \   d is the depth that we want to be proportional
                #   /   |   \  when plotted in the figure
                # ------------------
                #   region len
                #
                # d (in cm) =  depth (in bp) * 0.5 *
                #              width (in cm) / region len (in bp)
                # For hic_matrix_square:
                # d (in cm) = length of region2 (in bp) *
                #              width (in cm) / region len (in bp)
                #
                # to compute the actual width of the figure the margins
                # and the region
                # set for the legends have to be considered
                # DEFAULT_MARGINS['right'] - DEFAULT_MARGINS['left']
                # is the proportion of plotting area
                # This plotting area is divided in three part as specified in
                # self.width_ratios (normalized to 1)
                # And as wspace is specified 0.01,
                # 0.01 of the mean of the 3 regions is not occupied.
                # 1 / (1 + 2 / 3 * 0.01) is used to plot.

                hic_width = \
                    self.fig_width * \
                    (DEFAULT_MARGINS['right'] - DEFAULT_MARGINS['left']) / \
                    (1 + 2 / 3 * 0.01) * \
                    self.width_ratios[1] / sum(self.width_ratios)
                # the scale factor is to obtain each bin as a square
                # (a 45 degree rotated matrix)
                if track_dict['file_type'] == 'hic_matrix':
                    scale_factor = 0.5
                    depth = min(track_dict['depth'],
                                int((end_region - start_region) * 1.25))
                elif track_dict['file_type'] == 'hic_matrix_square':
                    scale_factor = 1  # No rotation
                    if track_dict['region2'] is None:
                        depth = end_region - start_region
                    else:
                        region2 = get_region(track_dict['region2'])
                        depth = region2[2] - region2[1]
                height = scale_factor * depth * hic_width / \
                    (end_region - start_region)
            else:
                height = DEFAULT_TRACK_HEIGHT
                self.track_list[i]['height'] = height

            track_height.append(height)

        return track_height

    def plot(self, file_name, chrom, start, end, title=None,
             h_align_titles='left', decreasing_x_axis=False):
        track_height = self.get_tracks_height(start_region=start,
                                              end_region=end)

        if self.fig_height:
            fig_height = self.fig_height
        else:
            fig_height = sum(track_height) / \
                (DEFAULT_MARGINS['top'] - DEFAULT_MARGINS['bottom'])

        log.debug(f"Figure size in cm is {self.fig_width} x {fig_height}."
                  f" Dpi is set to {self.dpi}\n")
        fig = plt.figure(figsize=self.cm2inch(self.fig_width, fig_height))

        fig.subplots_adjust(wspace=0, hspace=0.0,
                            left=DEFAULT_MARGINS['left'],
                            right=DEFAULT_MARGINS['right'],
                            bottom=DEFAULT_MARGINS['bottom'],
                            top=DEFAULT_MARGINS['top'])

        if title:
            fig.suptitle(title)

        grids = matplotlib.gridspec.GridSpec(len(track_height), 3,
                                             height_ratios=track_height,
                                             width_ratios=self.width_ratios,
                                             wspace=0.01)
        axis_list = []
        # skipped_tracks is the count of tracks that have the
        # 'overlay_previous' parameter and should be skipped
        skipped_tracks = 0
        plot_axis = None
        for idx, track in enumerate(self.track_obj_list):
            log.info(f"plotting {track.properties['section_name']}")

            if track.properties['overlay_previous'] in ['yes', 'share-y']:
                overlay = True
                skipped_tracks += 1
            else:
                overlay = False

            if track.properties['overlay_previous'] == 'share-y':
                ylim = plot_axis.get_ylim()
            else:
                idx -= skipped_tracks
                plot_axis = axisartist.Subplot(fig, grids[idx, 1])
                fig.add_subplot(plot_axis)
                # turns off the lines around the tracks
                plot_axis.axis[:].set_visible(False)
                # to make the background transparent
                plot_axis.patch.set_visible(False)
                if not overlay:
                    y_axis = plt.subplot(grids[idx, 0])
                    y_axis.set_axis_off()

                    label_axis = plt.subplot(grids[idx, 2])
                    label_axis.set_axis_off()
                    # I get the width of the label_axis to be able to wrap the
                    # labels when right or center aligned.
                    width_inch = label_axis.get_window_extent().width
                    width_dpi = width_inch * self.dpi / fig.dpi

            if decreasing_x_axis:
                plot_axis.set_xlim(end, start)
            else:
                plot_axis.set_xlim(start, end)
            track.plot(plot_axis, chrom, start, end)
            track.plot_y_axis(y_axis, plot_axis)
            track.plot_label(label_axis, width_dpi=width_dpi,
                             h_align=h_align_titles)

            if track.properties['overlay_previous'] == 'share-y':
                plot_axis.set_ylim(ylim)

            if not overlay:
                axis_list.append(plot_axis)

        for current_type in self.type_obj_list:
            current_type.plot(axis_list, fig, chrom, start, end)

        fig.savefig(file_name, dpi=self.dpi, transparent=False)
        return fig

    def parse_tracks(self, tracks_file, plot_regions=None):
        """
        Parses a configuration file

        :param tracks_file: file path containing the track configuration
        :param plot_regions: a list of tuple [(chrom1, start1, end1), (chrom2, start2, end2)]
                             on which the data should be loaded
                             here the vlines and the vhightlight
        :return: None
        """
        parser = ConfigParser(dict_type=MultiDict, strict=False)
        parser.read_file(open(tracks_file, 'r'))

        tracks_file_path = os.path.dirname(tracks_file)

        track_list = []
        type_list = []
        for section_name in parser.sections():
            # track_options is what will become the self.properties
            track_options = dict({"section_name": section_name})
            all_keywords = [i[0] for i in parser.items(section_name)]
            # First we check if there is a skip set to true:
            if 'skip' in all_keywords and \
               parser.getboolean(section_name, 'skip'):
                # In this case we just do not explore the section
                continue
            # We will append properties dictionnaries
            # to the track_list or type_list
            # Then the types are treated slightly differently:
            if 'type' in all_keywords and parser.get(section_name, 'type') in self.available_types:
                track_options['type'] = parser.get(section_name, 'type')
                track_options['track_class'] = \
                    self.available_types[track_options['type']]
            # If the sections are spacer or x-axis we fill the file_type:
            # (They are special sections where the title defines the track type)
            elif section_name.endswith('[spacer]'):
                track_options['file_type'] = 'spacer'
                track_options['track_class'] = SpacerTrack
            elif section_name.endswith('[x-axis]'):
                track_options['file_type'] = 'x_axis'
                track_options['track_class'] = XAxisTrack
            # For the others we need to have a 'file_type'
            # Either the file_type is part of the keywords
            elif 'file_type' in all_keywords:
                track_options['file_type'] = parser.get(section_name,
                                                        'file_type')
                if track_options['file_type'] not in self.available_tracks:
                    raise InputError(f"Section {section_name}: the file_type "
                                     f"{track_options['file_type']} does not"
                                     " exists.\npossible file_type are:"
                                     f"{self.available_tracks.keys()}.")
                track_options['track_class'] = \
                    self.available_tracks[track_options['file_type']]
            # Or we guess it from the file:
            elif 'file' in all_keywords:
                track_options['file'] = parser.get(section_name,
                                                   'file')
                track_options['file_type'] = \
                    self.guess_filetype(track_options,
                                        self.available_tracks)
                track_options['track_class'] = \
                    self.available_tracks[track_options['file_type']]
            else:
                raise InputError(f"Section {section_name}: there is no file_type nor file "
                                 "specified and it is not a [spacer] nor a "
                                 "[x-axis] section. This is not a valid "
                                 "section.")
            # Now we should have a 'track_class' set.
            # We can get for it all the necessary and possible keywords
            track_class = track_options['track_class']
            NECESSARY_PROPERTIES = track_class.NECESSARY_PROPERTIES
            for necessary_name in NECESSARY_PROPERTIES:
                if necessary_name not in all_keywords:
                    raise InputError(f"The section {section_name} is "
                                     "describing a object of"
                                     f" type {track_class} but the necessary "
                                     f"property {necessary_name}"
                                     " is not part of the config file.")
            unused_keys = []
            # Now we can proceed with the keywords:
            for name, value in parser.items(section_name):
                # To be removed in the next 1.0 version
                if ' ' in name:
                    old_name = name
                    name = '_'.join(name.split(' '))
                    log.warning(f"Deprecated Warning: The section {section_name} "
                                f"uses parameter {old_name} but there is no more "
                                "parameter with space in name. "
                                f"Will be substituted by {name}.\n")
                else:
                    old_name = name
                # end
                SYNONYMOUS_PROPERTIES = track_class.SYNONYMOUS_PROPERTIES
                # If the name is part of the synonymous we substitute by
                # the synonymous value
                if name in SYNONYMOUS_PROPERTIES and \
                   value in SYNONYMOUS_PROPERTIES[name]:
                    track_options[name] = SYNONYMOUS_PROPERTIES[name][value]
                elif name in track_class.STRING_PROPERTIES:
                    track_options[name] = value
                elif name in track_class.BOOLEAN_PROPERTIES:
                    try:
                        # I need to use old_name here else I get a KeyError:
                        track_options[name] = parser.getboolean(section_name,
                                                                old_name)
                        # In the next 1.0 should be:
                        # track_options[name] = parser.getboolean(section_name,
                        #                                         name)
                    except ValueError:
                        raise InputError(f"In section {section_name}, "
                                         f"{old_name} was set to {value}"
                                         " whereas we should have a boolean "
                                         "value. Please, use true or false.")
                        # In the next 1.0 should be:
                        #                f"{name} was set to {value}"
                    if value.lower() not in ['true', 'false']:
                        log.warning("Deprecation Warning: "
                                    f"In section {section_name}, {name} was "
                                    f"set to {value}"
                                    " whereas in the future only"
                                    " true and false value will be"
                                    " accepted.\n")
                elif name in track_class.FLOAT_PROPERTIES:
                    try:
                        track_options[name] = float(value)
                    except ValueError:
                        raise InputError(f"In section {section_name}, {name} "
                                         f"was set to {value}"
                                         " whereas we should have a float "
                                         "value.")
                    min_value, max_value = track_class.FLOAT_PROPERTIES[name]
                    if track_options[name] < min_value or \
                       track_options[name] > max_value:
                        raise InputError(f"In section {section_name}, {name} "
                                         f"was set to {value}"
                                         " whereas it should be between "
                                         f"{min_value} and {max_value}.")
                elif name in track_class.INTEGER_PROPERTIES:
                    try:
                        track_options[name] = int(value)
                    except ValueError:
                        raise InputError(f"In section {section_name}, {name} "
                                         f"was set to {value}"
                                         " whereas we should have an integer "
                                         "value.")
                    min_value, max_value = track_class.INTEGER_PROPERTIES[name]
                    if track_options[name] < min_value or \
                       track_options[name] > max_value:
                        raise InputError(f"In section {section_name}, {name} "
                                         f"was set to {value}"
                                         " whereas it should be between "
                                         f"{min_value} and {max_value}.")
                else:
                    unused_keys.append(name)
            # If there are unused keys they are printed in a warning.
            if len(unused_keys) > 0:
                log.warning(f"In section {section_name}, these parameters are "
                            f"unused:{unused_keys}.\n")
            # The track_options will be checked for the file paths:
            track_options = self.check_file_exists(track_options,
                                                   tracks_file_path,
                                                   track_options.get('file_type', 'no') in
                                                   ['hic_matrix', 'hic_matrix_square'])
            if 'file_type' in track_options:
                # The 'overlay_previous' is initialized:
                if 'overlay_previous' not in track_options:
                    track_options['overlay_previous'] = 'no'
                if track_options['overlay_previous'] not in ['no', 'yes', 'share-y']:
                    raise InputError(f"In section {section_name}, overlay_previous "
                                     f"was set to {track_options['overlay_previous']}."
                                     " Possible options are no, yes, share-y")
                # If there is no title:
                if 'title' not in track_options:
                    track_options['title'] = ''
                    if track_options['overlay_previous'] == 'no' and \
                        track_options['track_class'] not in [SpacerTrack,
                                                             XAxisTrack]:
                        log.warning("title not set for section "
                                    f"{track_options['section_name']}\n")
                # The track_options are added to the track_list
                track_list.append(track_options)
            else:
                type_list.append(track_options)
        # Now that they were all checked
        self.track_list = track_list
        self.type_list = type_list

    def close_files(self):
        """
        Close all opened files
        """
        for track in self.track_obj_list:
            track.__del__()

    @staticmethod
    def check_file_exists(track_dict, tracks_path, is_hic=False):
        """
        Checks if a file or list of files exists. If the file does not exists
        tries to check if the file may be relative to the track_file path,
        in such case the path is updated.
        :param track_dict: dictionary of track values. Should contain
                            a 'file' key containing the path of the file
                            or files to be checked separated by space
                            For example: file1 file2 file3
        :param tracks_path: path of the tracks file
        :return: None
        """
        for key in track_dict.keys():
            if key.endswith("file"):
                file_field_name = key
                # # THIS COULD BE REMOVED IN A NEXT 1.0 VERSION
                if file_field_name == 'boundaries_file':
                    log.warning("Deprecation Warning: "
                                "The boundaries_file is not used anymore."
                                " It will be ignored."
                                " Please use another track with the"
                                " `overlay_previous` option.\n")
                # # END
                file_names = [x for x in track_dict[file_field_name].split(" ") if x != '']
                full_path_file_names = []
                for file_name in file_names:
                    if is_hic and not file_name.endswith('.h5'):
                        file_name_to_check = file_name.split("::")[0]
                    else:
                        file_name_to_check = file_name

                    try:
                        open(file_name_to_check, 'r').close()
                        full_path_file_names.append(file_name)
                    except IOError:
                        try:
                            # try to find the file in the same path as the
                            # track file
                            name_with_tracks_path = tracks_path + "/" + file_name
                            name_with_tracks_path_to_check = tracks_path + "/" + file_name_to_check
                            open(name_with_tracks_path_to_check, 'r').close()
                            full_path_file_names.append(name_with_tracks_path)
                        except IOError:
                            raise InputError(f"File in section [{track_dict['section_name']}] "
                                             f"not found:\n{file_name}\n\n")

                track_dict[file_field_name] = " ".join(full_path_file_names)
        return track_dict

    @staticmethod
    def guess_filetype(track_dict, available_tracks):
        """

        :param track_dict: dictionary of track values with the 'file' key
                    containing a string path of the file or files.
                    Only the ending of the last file is used
                    in case when there are more files
        :param: available_tracks: list of available tracks

        :return: string file type detected
        """
        file_ = track_dict['file'].strip()
        file_type = None
        for track_type, track_class in available_tracks.items():
            for ending in track_class.SUPPORTED_ENDINGS:
                if file_.endswith(ending):
                    if file_type == track_class.TRACK_TYPE:
                        raise InputError("file_type already defined in other"
                                         " GenomeTrack")
                    else:
                        file_type = track_class.TRACK_TYPE

        if file_type is None:
            raise InputError(f"Section {track_dict['section_name']}: can not "
                             "identify file type. Please"
                             f" specify the file_type for '{file_}'")

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
        if start:
            log.info(time.time() - start)
        return time.time()


class SpacerTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = 'spacer'
    DEFAULTS_PROPERTIES = {}
    NECESSARY_PROPERTIES = []
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['overlay_previous',
                         'title', 'file_type']
    FLOAT_PROPERTIES = {'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def plot(self, ax, chrom_region, start_region, end_region):
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax)

    def plot_y_axis(self, ax, plot_ax):
        pass


class XAxisTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = 'x_axis'
    NECESSARY_PROPERTIES = []
    DEFAULTS_PROPERTIES = {'where': 'bottom',
                           'fontsize': 15}
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'where': ['top', 'bottom']}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['overlay_previous',
                         'title', 'where', 'file_type']
    FLOAT_PROPERTIES = {'fontsize': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwargs):
        super(XAxisTrack, self).__init__(*args, **kwargs)

    def plot(self, ax, chrom_region, region_start, region_end):

        start, end = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax)
        if self.properties['where'] == 'top':
            ax.axis["x"] = ax.new_floating_axis(0, 0.2, axis_direction="top")
            # In some old matplotlib version (<3.3.4) it was not set:
            ax.axis["x"].set_axis_direction("top")
            label_y_pos = 0.99
            vert_align = 'top'
        else:
            ax.axis["x"] = ax.new_floating_axis(0, 0.9)
            label_y_pos = 0.01
            vert_align = 'bottom'

        # First adjust the size of the label
        ax.axis["x"].major_ticklabels.set(size=self.properties['fontsize'])

        # Use the formatter
        ax.xaxis.set_major_formatter(MyBasePairFormatter())
        ax.text(0.5, label_y_pos, chrom_region, horizontalalignment='center',
                fontsize=self.properties['fontsize'],
                verticalalignment=vert_align, transform=ax.transAxes)

    def plot_y_axis(self, ax, plot_ax):
        pass
