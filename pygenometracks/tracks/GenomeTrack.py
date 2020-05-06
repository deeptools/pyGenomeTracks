# -*- coding: utf-8 -*-

from .. utilities import to_string, to_bytes
import logging
import numpy as np
from matplotlib import colors as mc
import matplotlib.pyplot as plt


class GenomeTrack(object):
    """
    The GenomeTrack object is a holder for all tracks that are to be plotted.
    For example, to plot a bedgraph file a new class that extends GenomeTrack
    should be created.

    It is expected that all GenomeTrack objects have a plot method.

    """
    SUPORTED_ENDINGS = []
    TRACK_TYPE = None
    OPTIONS_TXT = """
# title of track (plotted on the right side)
title =
# height of track in cm (ignored if the track is overlay on top the previous track)
height = 2
# if you want to plot the track upside-down:
# orientation = inverted
# if you want to plot the track on top of the previous track. Options are 'yes' or 'share-y'. For the 'share-y'
# option the y axis values is shared between this plot and the overlay plot. Otherwise, each plot use its own scale
#overlay_previous = yes
"""
    DEFAULTS_PROPERTIES = {}
    NECESSARY_PROPERTIES = []
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['file_type', 'orientation',  # For XAxisTrack and SpacerTrack these 2 are not used
                         'overlay_previous', 'title']
    FLOAT_PROPERTIES = {'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, properties_dict):
        FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
        logging.basicConfig(format=FORMAT)
        log = logging.getLogger(__name__)
        log.setLevel(logging.DEBUG)
        self.log = log
        self.properties = properties_dict
        self.set_properties_defaults()
        self.file_type = 'test'

    def set_properties_defaults(self):
        # put to default all properties which are not set:
        for prop in self.DEFAULTS_PROPERTIES:
            if prop not in self.properties:
                self.properties[prop] = self.DEFAULTS_PROPERTIES[prop]
        # use synonymous:
        for prop in self.SYNONYMOUS_PROPERTIES:
            synonymous = self.SYNONYMOUS_PROPERTIES[prop]
            if prop in self.properties and \
               self.properties[prop] in synonymous:
                self.properties[prop] = synonymous[self.properties[prop]]
        # check if properties are possible:
        for prop in self.POSSIBLE_PROPERTIES:
            possibles = self.POSSIBLE_PROPERTIES[prop]
            if self.properties[prop] not in possibles:
                default_value = self.DEFAULTS_PROPERTIES[prop]
                self.log.warning("*WARNING* {0}: '{1}' for section {2}"
                                 " is not valid. {0} has "
                                 "been set to "
                                 "{3}".format(prop,
                                              self.properties[prop],
                                              self.properties['section_name'],
                                              default_value))
                self.properties[prop] = default_value

    def plot_y_axis(self, ax, plot_axis, transform='no', log_pseudocount=0,
                    y_axis='tranformed'):
        """
        Plot the scale of the y axis with respect to the plot_axis
        Args:
            ax: axis to use to plot the scale
            plot_axis: the reference axis to get the max and min.
            transform: what was the transformation of the data
            log_pseudocount:
            y_axis: 'tranformed' or 'original'

        Returns:

        """
        if not self.properties.get('show_data_range', True):
            return

        def value_to_str(value):
            # given a numeric value, returns a
            # string that removes unneeded decimal places
            if value % 1 == 0:
                str_value = str(int(value))
            else:
                str_value = "{:.1f}".format(value)
            return str_value

        ymin, ymax = plot_axis.get_ylim()

        if transform == 'no' or y_axis == 'transformed':
            # This is a linear scale
            # plot something that looks like this:
            # ymax ┐
            #      │
            #      │
            # ymin ┘

            # the coordinate system used is the ax.transAxes (lower left corner (0,0), upper right corner (1,1)
            # this way is easier to adjust the positions such that the lines are plotted complete
            # and not only half of the width of the line.
            x_pos = [0, 0.5, 0.5, 0]
            y_pos = [0.01, 0.01, 0.99, 0.99]
        else:
            # ymid is the middle between ymin and ymax
            # if 0 is between ymin and ymax then ymid is 0
            if ymin * ymax < 0:
                ymid = 0
                ymid_pos = - ymin / (ymax - ymin)
            else:
                ymid = (ymin + ymax) / 2
                ymid_pos = 0.5

            if transform == 'log':
                ymin, ymid, ymax = np.exp([ymin, ymid, ymax]) - log_pseudocount
            elif transform == 'log2':
                ymin, ymid, ymax = np.exp2([ymin, ymid, ymax]) - \
                    log_pseudocount
            elif transform == 'log10':
                ymin, ymid, ymax = np.power(10, [ymin, ymid, ymax]) - \
                    log_pseudocount
            elif transform == 'log1p':
                ymin, ymid, ymax = np.expm1([ymin, ymid, ymax])
            elif transform == '-log':
                ymin, ymid, ymax = np.exp(- [ymin, ymid, ymax]) - \
                    log_pseudocount
            ymid_str = value_to_str(ymid)
            # plot something that looks like this:
            # ymax ┐
            #      │
            # ymid-|
            #      │
            # ymin ┘
            x_pos = [0, 0.5, 0.5, 0, 0.5, 0.5, 0]
            y_pos = [0.01, 0.01, ymid_pos, ymid_pos, ymid_pos, 0.99, 0.99]
        ymax_str = value_to_str(ymax)
        ymin_str = value_to_str(ymin)
        ax.plot(x_pos, y_pos, color='black', linewidth=1, transform=ax.transAxes)
        ax.text(-0.2, -0.01, ymin_str, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes)
        ax.text(-0.2, 1, ymax_str, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)
        if transform != 'no':
            if y_axis == 'original':
                ax.text(-0.2, ymid_pos, ymid_str, verticalalignment='center',
                        horizontalalignment='right', transform=ax.transAxes)
            else:
                if transform == 'log1p':
                    ymid_str = "log(1 + x)"
                else:
                    if log_pseudocount == 0:
                        ymid_str = "{}(x)".format(transform)
                    else:
                        ymid_str = "{}({} + x)".format(transform,
                                                       log_pseudocount)

                ax.text(0, 0.5, ymid_str, verticalalignment='center',
                        horizontalalignment='right', transform=ax.transAxes,
                        wrap=True)

        ax.patch.set_visible(False)

    def plot_label(self, label_ax, width_dpi, h_align='left'):
        if h_align == 'left':
            label_ax.text(0.05, 0.5, self.properties['title'],
                          horizontalalignment='left', size='large',
                          verticalalignment='center',
                          transform=label_ax.transAxes,
                          wrap=True)
        elif h_align == 'right':
            txt = label_ax.text(1, 0.5, self.properties['title'],
                                horizontalalignment='right', size='large',
                                verticalalignment='center',
                                transform=label_ax.transAxes,
                                wrap=True)
            # To be able to wrap to the left:
            txt._get_wrap_line_width = lambda: width_dpi
        else:
            txt = label_ax.text(0.5, 0.5, self.properties['title'],
                                horizontalalignment='center', size='large',
                                verticalalignment='center',
                                transform=label_ax.transAxes,
                                wrap=True)
            # To be able to wrap to the left:
            txt._get_wrap_line_width = lambda: width_dpi

    def process_type_for_coverage_track(self):
        default_plot_type = 'fill'
        self.plot_type = default_plot_type
        self.size = None

        if self.properties['type'].find(":") > 0:
            self.plot_type, size = self.properties['type'].split(":")
            try:
                self.size = float(size)
            except ValueError:
                self.log.warning("Invalid value: 'type = {}' in section: {}\n"
                                 "A number was expected after ':' and found "
                                 "'{}'. Will use default."
                                 "".format(self.properties['type'],
                                           self.properties['section_name'],
                                           size))
        else:
            self.plot_type = self.properties['type']

        if self.plot_type not in ['line', 'points', 'fill']:
            self.log.warning("Invalid: 'type = {}' in section: {}\n"
                             "Will use default."
                             "".format(self.properties['type'],
                                       self.properties['section_name']))
            self.plot_type = default_plot_type

    def process_color(self, param, colormap_possible=False,
                      bed_rgb_possible=False, colormap_only=False,
                      default_value_is_colormap=False):
        """
        Put a valid color/colormap in self.properties[param]
        Args:
            param: param to check/update the value
            colormap_possible: if the self.properties[param] can be a colormap
            bed_rgb_possible: if the self.properties[param] can be 'bed_rgb'
            colormap_only: if the self.properties[param] must be a colormap

        Returns:
            True if the self.properties[param] is a colormap
            False if the self.properties[param] is not a colormap

        """
        default_value = self.DEFAULTS_PROPERTIES[param]
        valid_color = None
        if bed_rgb_possible and self.properties[param] == 'bed_rgb':
            return False
        if mc.is_color_like(self.properties[param]):
            valid_color = self.properties[param]
        # It can be a tuple (for example (1, 0.88, 2./3) would be a valid color):
        elif self.properties[param][0] == '(':
            try:
                custom_color = eval(self.properties[param])
            except ValueError:
                self.log.warning("*WARNING* {}: '{}' for section {}"
                                 " is not valid. {} has "
                                 "been set to "
                                 "{}".format(param,
                                             self.properties[param],
                                             self.properties['section_name'],
                                             default_value))
                valid_color = default_value
            else:
                if mc.is_color_like(custom_color):
                    valid_color = custom_color
                else:
                    self.log.warning("*WARNING* {}: '{}' for section {}"
                                     " is not valid. {} has "
                                     "been set to "
                                     "{}".format(param,
                                                 self.properties[param],
                                                 self.properties['section_name'],
                                                 default_value))
                    valid_color = default_value
        if not colormap_possible:
            if valid_color is None:
                self.log.warning("*WARNING* {}: '{}' for section {}"
                                 " is not valid. It has "
                                 "been set to "
                                 "{}".format(param,
                                             self.properties[param],
                                             self.properties['section_name'],
                                             default_value))
                valid_color = default_value
            self.properties[param] = valid_color
            return False
        else:
            valid_colormap = None
            # We will try to process the color as a colormap
            if valid_color is None:
                # If someone what to use its own colormap,
                # he can specify the rgb values or color values:
                # For example:
                # colormap = ['white', (1, 0.88, 2./3), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
                if self.properties[param][0] == '[':
                    try:
                        custom_colors = eval(self.properties[param])
                    except (SyntaxError, NameError) as e:
                        self.log.warning("Warning: section {}, {} was set as {} but "
                                         "raises an error:\n{}\nIt will be ignored and"
                                         " default value will be used."
                                         "".format(self.properties['section_name'],
                                                   param, self.properties[param],
                                                   e))
                    else:
                        try:
                            valid_colormap = mc.LinearSegmentedColormap.from_list(
                                'custom', custom_colors, N=100)
                        except ValueError as e:
                            self.log.warning("Warning: section {}, {} was set as {} but "
                                             "raises an error:\n{}\nIt will be ignored and"
                                             " default value will be used."
                                             "".format(self.properties['section_name'],
                                                       param, self.properties[param],
                                                       e))
                else:
                    if self.properties[param] in dir(plt.cm):
                        valid_colormap = self.properties[param]
        # Here, colormap is possible
        # valid_color is None or a valid color or the default value
        # valid_colormap is None or a valid colormap
        if valid_color is None and valid_colormap is None:
            self.log.warning("*WARNING* {}: '{}' for section {}"
                             " is not valid. It has "
                             "been set to "
                             "{}".format(param,
                                         self.properties[param],
                                         self.properties['section_name'],
                                         default_value))
            self.properties[param] = default_value
            return default_value_is_colormap
        if colormap_only:
            if valid_colormap is None:
                # valid_color is not None
                self.log.warning("*WARNING* {}: '{}' for section {}"
                                 " is not valid: it is a color whereas"
                                 " a colormap was expected. It has "
                                 "been set to "
                                 "{}".format(param,
                                             self.properties[param],
                                             self.properties['section_name'],
                                             default_value))
                valid_colormap = default_value
            self.properties[param] = valid_colormap
            return True
        else:
            if valid_color is not None:
                self.properties[param] = valid_color
                return False
            else:
                self.properties[param] = valid_colormap
                return True

    @staticmethod
    def change_chrom_names(chrom):
        """
        Changes UCSC chromosome names to ensembl chromosome names
        and vice versa.
        """
        # TODO: mapping from chromosome names like mithocondria is missing
        if chrom.startswith('chr'):
            # remove the chr part from chromosome name
            chrom = chrom[3:]
        else:
            # prefix with 'chr' the chromosome name
            chrom = 'chr' + chrom

        return chrom

    @staticmethod
    def check_chrom_str_bytes(iteratable_obj, p_obj):
        # determine type
        if isinstance(p_obj, list) and len(p_obj) > 0:
            type_ = type(p_obj[0])
        else:
            type_ = type(p_obj)
        if not isinstance(type(next(iter(iteratable_obj))), type_):
            if type(next(iter(iteratable_obj))) is str:
                p_obj = to_string(p_obj)
            elif type(next(iter(iteratable_obj))) in [bytes, np.bytes_]:
                p_obj = to_bytes(p_obj)
        return p_obj

    def __del__(self):
        return
