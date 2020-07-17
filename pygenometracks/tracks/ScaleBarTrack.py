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
# x_center: coordinate where the scale bar should be plotted (center)
# if not set it will be in the middle of the plotted area
#x_center = 3100000
# size: in bp the length of the scale bar
# if not set it will be like in UCSC:
# the higher number that begins with 1, 2 or 5 followed by 0s
# that is less than half the plotted area
#size = 100000
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
                           'where': 'left'}
    NECESSARY_PROPERTIES = []
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'where': ['left', 'right',
                                     'top', 'bottom']}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['file_type', 'title', 'color', 'where']
    FLOAT_PROPERTIES = {'fontsize': [0, np.inf],
                        'alpha': [0, 1],
                        'line_width': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {'x_center': [0, np.inf],
                          'size': [0, np.inf]}
    # The color can only be a color

    def plot(self, ax, chrom_region, start_region, end_region):
        # Get the center from the properties
        x_center = self.properties['x_center']
        if x_center is None:
            # Else put it in the middle
            x_center = (end_region + start_region) / 2
        # Get the size form the properties
        size = self.properties['size']
        if size is None:
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
            size_label = f"{size} bases"
        elif size < 1e6:
            new_size = size / 1e3
            size_label = f"{int(new_size) if new_size.is_integer() else new_size} kb"
        else:
            new_size = size / 1e6
            size_label = f"{int(new_size) if new_size.is_integer() else new_size} Mb"

        # We now draw the |-----|
        x_left = x_center - size / 2
        x_right = x_center + size / 2
        xdatas = [(x_left, x_left),  # left line
                  (x_left, x_right),  # horizontal line
                  (x_right, x_right)]  # right line
        ydatas = [(0, 1),  # left line
                  (0.5, 0.5),  # horizontal line
                  (0, 1)]  # right line
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
            if where == 'left':
                ax.text(x=x_left - 2 * small_space, y=0.5, s=size_label,
                        horizontalalignment='right',
                        verticalalignment='center', fontproperties=fp)
            else:
                ax.text(x=x_right + 2 * small_space, y=0.5, s=size_label,
                        horizontalalignment='left',
                        verticalalignment='center', fontproperties=fp)

        return ax

    def plot_y_axis(self, ax, plot_axis):
        return
