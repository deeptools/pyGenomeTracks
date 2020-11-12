from . GenomeTrack import GenomeTrack
import numpy as np
from .. utilities import InputError

DEFAULT_HLINES_COLOR = 'black'


class HLinesTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = 'hlines'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# color of the lines
color = black
# To use transparency, you can use alpha
# default is 1
# alpha = 0.5
# the default for min_value and max_value is 'auto' which means that the scale will go
# roughly from the minimum value found in the region plotted to the maximum value found.
min_value = 0
#max_value = auto
# line width:
# line_width = 0.5
# options for line_style are 'solid', 'dashed', 'dotted', and 'dashdot'
#line_style = solid
# y values where horizontal lines should be plotted separated by comma:
y_values = 10, 200
# set show_data_range to false to hide the text on the upper-left showing the data range
show_data_range = true
file_type = {TRACK_TYPE}
    """
    DEFAULTS_PROPERTIES = {'max_value': None,
                           'min_value': None,
                           'show_data_range': True,
                           'orientation': None,
                           'color': DEFAULT_HLINES_COLOR,
                           'alpha': 1,
                           'line_width': 0.5,
                           'line_style': 'solid'}
    NECESSARY_PROPERTIES = ['y_values']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           'line_style': ['solid', 'dashed',
                                          'dotted', 'dashdot']}
    BOOLEAN_PROPERTIES = ['show_data_range']
    STRING_PROPERTIES = ['y_values', 'file_type', 'overlay_previous',
                         'orientation', 'line_style',
                         'title', 'color']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'alpha': [0, 1],
                        'line_width': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}
    # The color can only be a color

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        y_array = self.properties['y_values'].split(',')
        try:
            self.y_values = [float(x) for x in y_array]
        except Exception as detail:
            raise InputError(f"y_values ({self.properties['y_values']}) not"
                             f" valid. \nError: {detail}")

    def set_properties_defaults(self):
        super(HLinesTrack, self).set_properties_defaults()
        self.process_color('color')

    def plot(self, ax, chrom_region, start_region, end_region):
        self.log.debug(f"y_values: {self.y_values}")
        for y_value in self.y_values:
            self.log.debug(f"y_value: {y_value}")
            ax.axhline(y=y_value,
                       linewidth=self.properties['line_width'],
                       color=self.properties['color'],
                       alpha=self.properties['alpha'],
                       linestyle=self.properties['line_style'])
        self.adjust_ylim(ax)

        return ax
