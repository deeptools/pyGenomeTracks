from . GenomeTrack import GenomeTrack
import numpy as np
from .. utilities import get_length_w
from matplotlib import font_manager
from matplotlib.lines import Line2D

DEFAULT_SCALEBAR_COLOR = 'black'


class ScaleBarTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = 'scalebar'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# color of the scalebar
color = black
# To use transparency, you can use alpha
# default is 1
# alpha = 0.5
# line width:
# line_width = 0.5
# To set the position and the size of the scale bar
# 4 parameters can be set:
# x_center, scalebar_start_position, scalebar_end_position, size
# scalebar_start_position need to be smaller than x_center smaller than scalebar_end_position
# x_center: coordinate where the scale bar should be plotted (center)
# if not set and cannot be deduce from other parameters
# it will be in the middle of the plotted area
#x_center = 3100000
# size: in bp the length of the scale bar
# if not set and cannot be deduced from other parameters
# it will be like in UCSC:
# the higher number that begins with 1, 2 or 5 followed by 0s
# that is less than half the plotted area
#size = 100000
# Another example:
#scalebar_start_position = 2900000
#scalebar_end_position = 3100000
# where: where the size of the scale bar should
# appear among left, right, top, bottom
# default is left
#where = right
# fontsize: default is 12
#fontsize = 10
file_type = {TRACK_TYPE}
    """
    DEFAULTS_PROPERTIES = {'fontsize': 12,
                           'color': DEFAULT_SCALEBAR_COLOR,
                           'alpha': 1,
                           'line_width': 0.5,
                           'x_center': None,
                           'size': None,
                           'scalebar_start_position': None,
                           'scalebar_end_position': None,
                           'where': 'left'}
    NECESSARY_PROPERTIES = []
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'where': ['left', 'right',
                                     'top', 'bottom']}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['file_type', 'title', 'overlay_previous',
                         'color', 'where']
    FLOAT_PROPERTIES = {'fontsize': [0, np.inf],
                        'alpha': [0, 1],
                        'line_width': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {'x_center': [0, np.inf],
                          'scalebar_start_position': [0, np.inf],
                          'scalebar_end_position': [0, np.inf],
                          'size': [0, np.inf]}
    # The color can only be a color

    def set_properties_defaults(self):
        super(ScaleBarTrack, self).set_properties_defaults()
        self.process_color('color')
        # Check that all parameters for position are compatible and set x_center and size:
        param_set = [param for param in ['x_center', 'size',
                                         'scalebar_start_position',
                                         'scalebar_end_position']
                     if self.properties[param] is not None]
        # Check that the order is respected: start < center < end
        if 'scalebar_start_position' in param_set and 'scalebar_end_position' in param_set:
            assert self.properties['scalebar_start_position'] < self.properties['scalebar_end_position'], \
                f"scalebar_start_position ({self.properties['scalebar_start_position']})" \
                f"must be smaller than scalebar_end_position ({self.properties['scalebar_end_position']})."
        if 'scalebar_start_position' in param_set and 'x_center' in param_set:
            assert self.properties['scalebar_start_position'] < self.properties['x_center'], \
                f"scalebar_start_position ({self.properties['scalebar_start_position']})." \
                f"must be smaller than x_center ({self.properties['x_center']})"
        if 'x_center' in param_set and 'scalebar_end_position' in param_set:
            assert self.properties['x_center'] < self.properties['scalebar_end_position'], \
                f"x_center ({self.properties['x_center']})" \
                f"must be smaller than scalebar_end_position ({self.properties['scalebar_end_position']})."
        if len(param_set) == 4:
            # There is a risk of incompatibility
            computed_size = abs(self.properties['scalebar_end_position'] - self.properties['scalebar_start_position'])
            assert self.properties['size'] == computed_size, \
                f"`size` ({self.properties['size']}) is incompatible with "\
                f"`scalebar_start_position` and `scalebar_end_position`, should be {computed_size}"
            computed_center = (self.properties['scalebar_start_position'] + self.properties['scalebar_end_position']) / 2
            assert self.properties['x_center'] == computed_center, \
                f"`x_center` ({self.properties['x_center']}) is incompatible with "\
                f"`scalebar_start_position` and `scalebar_end_position`, should be {computed_center}"
        elif len(param_set) == 3:
            # There is a risk of incompatibility
            # I need to set x_center and size when needed
            if self.properties['x_center'] is None:
                computed_size = self.properties['scalebar_end_position'] - self.properties['scalebar_start_position']
                assert self.properties['size'] == computed_size, \
                    f"`size` ({self.properties['size']}) is incompatible with "\
                    f"`scalebar_start_position` and `scalebar_end_position`, should be {computed_size}"
                # I add x_center to properties:
                self.properties['x_center'] = (self.properties['scalebar_start_position'] + self.properties['scalebar_end_position']) / 2
            elif self.properties['size'] is None:
                computed_center = (self.properties['scalebar_start_position'] + self.properties['scalebar_end_position']) / 2
                assert self.properties['x_center'] == computed_center, \
                    f"`x_center` ({self.properties['x_center']}) is incompatible with "\
                    f"`scalebar_start_position` and `scalebar_end_position`, should be {computed_center}"
                # I add size to properties:
                self.properties['size'] = self.properties['scalebar_end_position'] - self.properties['scalebar_start_position']
            elif self.properties['scalebar_start_position'] is None:
                computed_size = 2 * (self.properties['scalebar_end_position'] - self.properties['x_center'])
                assert self.properties['size'] == computed_size, \
                    f"`size` ({self.properties['size']}) is incompatible with "\
                    f"`x_center` and `scalebar_end_position`, should be {computed_size}"
            else:  # self.properties['scalebar_start_position'] is None:
                computed_size = 2 * (self.properties['x_center'] - self.properties['scalebar_start_position'])
                assert self.properties['size'] == computed_size, \
                    f"`size` ({self.properties['size']}) is incompatible with "\
                    f"`scalebar_start_position` and `x_center`, should be {computed_size}"

    def plot(self, ax, chrom_region, start_region, end_region):
        # Get center and size from properties
        # Get the center from the properties
        x_center = self.properties['x_center']
        if x_center is None:
            # Try to guess it from other parameters:
            if self.properties['scalebar_start_position'] is not None and self.properties['scalebar_end_position'] is not None:
                # I deduce x_center:
                x_center = (self.properties['scalebar_start_position'] + self.properties['scalebar_end_position']) / 2
                # I add size to properties:
                self.properties['size'] = self.properties['scalebar_end_position'] - self.properties['scalebar_start_position']
            elif self.properties['scalebar_start_position'] is not None and self.properties['size'] is not None:
                x_center = self.properties['scalebar_start_position'] + self.properties['size'] / 2
            elif self.properties['scalebar_end_position'] is not None and self.properties['size'] is not None:
                x_center = self.properties['scalebar_end_position'] - self.properties['size'] / 2
            else:
                # Else put it in the middle
                x_center = (end_region + start_region) / 2
        # Get the size form the properties
        size = self.properties['size']
        if size is None:
            # Try to guess it from other parameters:
            if self.properties['scalebar_start_position'] is not None:
                size = 2 * (x_center - self.properties['scalebar_start_position'])
            elif self.properties['scalebar_end_position'] is not None:
                size = 2 * (self.properties['scalebar_end_position'] - x_center)
            else:
                # We put the size that is less than half the plotted region
                # Which begins with 1, 2 or 5 followed by 0s
                half_plotted_region = int((end_region - start_region) / 2)
                first_char_hpr = int(str(half_plotted_region)[0])
                if first_char_hpr >= 5:
                    first_char = 5
                elif first_char_hpr >= 2:
                    first_char = 2
                else:
                    first_char = 1
                size = first_char * 10**(len(str(half_plotted_region)) - 1)
        # We only plot if it will be visible
        if x_center - size / 2 > end_region or x_center + size / 2 < start_region:
            return

        # We adjust the unit to make it pretty
        if size < 1e3:
            size_label = f"{int(size) if float(size).is_integer() else size} bases"
        elif size < 1e6:
            new_size = size / 1e3
            size_label = f"{int(new_size) if new_size.is_integer() else new_size} kb"
        else:
            new_size = size / 1e6
            size_label = f"{int(new_size) if new_size.is_integer() else new_size} Mb"

        # We now draw the |-----|
        scalebar_start_position = x_center - size / 2
        scalebar_end_position = x_center + size / 2
        xdatas = [(scalebar_start_position, scalebar_start_position),  # vertical line start
                  (scalebar_start_position, scalebar_end_position),  # horizontal line
                  (scalebar_end_position, scalebar_end_position)]  # vertical line end
        ydatas = [(0, 1),  # vertical line start
                  (0.5, 0.5),  # horizontal line
                  (0, 1)]  # vertical line end
        for i in range(3):
            ax.add_line(Line2D(xdatas[i], ydatas[i],
                        color=self.properties['color'],
                        linewidth=self.properties['line_width'],
                        alpha=self.properties['alpha']))

        # Process the font size:
        fp = font_manager.FontProperties(size=self.properties['fontsize'])
        # We now write the size_label:
        where = self.properties['where']
        if where in ['top', 'bottom']:
            pt_per_cm = 28.3464567
            font_height_cm = self.properties['fontsize'] / pt_per_cm
            relative_height = font_height_cm / self.properties['height']
            if where == 'top':
                ax.text(x=x_center, y=1, s=size_label,
                        horizontalalignment='center',
                        verticalalignment='bottom', fontproperties=fp)
                ax.set_ylim(0, 1 / (1 - relative_height))
            else:
                ax.text(x=x_center, y=0, s=size_label,
                        horizontalalignment='center',
                        verticalalignment='top', fontproperties=fp)
                ax.set_ylim(1 - 1 / (1 - relative_height), 1)
        else:
            small_space = get_length_w(ax.get_figure().get_figwidth(),
                                       start_region, end_region,
                                       self.properties['fontsize'])
            xlim = ax.get_xlim()
            if where == 'left' and xlim[0] < xlim[1]:
                ax.text(x=scalebar_start_position - 2 * small_space, y=0.5, s=size_label,
                        horizontalalignment='right',
                        verticalalignment='center', fontproperties=fp)
            elif where == 'left' and xlim[0] >= xlim[1]:
                ax.text(x=scalebar_end_position + 2 * small_space, y=0.5, s=size_label,
                        horizontalalignment='right',
                        verticalalignment='center', fontproperties=fp)
            elif where == 'right' and xlim[0] >= xlim[1]:
                ax.text(x=scalebar_start_position - 2 * small_space, y=0.5, s=size_label,
                        horizontalalignment='left',
                        verticalalignment='center', fontproperties=fp)
            else:  # where == 'left' and xlim[0] < xlim[1]
                ax.text(x=scalebar_end_position + 2 * small_space, y=0.5, s=size_label,
                        horizontalalignment='left',
                        verticalalignment='center', fontproperties=fp)

        return ax

    def plot_y_axis(self, ax, plot_axis):
        return
