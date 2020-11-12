from . GenomeTrack import GenomeTrack
from .. readBed import ReadBed
# To remove next 1.0
from .. readGtf import ReadGtf
# End to remove
from .. utilities import opener, count_lines, temp_file_from_intersect
import matplotlib
from matplotlib.patches import Polygon
from intervaltree import IntervalTree, Interval
import numpy as np
from tqdm import tqdm

DEFAULT_BED_COLOR = '#1f78b4'
AROUND_REGION = 100000


class BedLikeTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = None
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# optional: line_width
#line_width = 0.5
# optional, default is black. To remove the border, simply set 'border_color' to none
#border_color = black
"""

    DEFAULTS_PROPERTIES = {'orientation': None,
                           'color': DEFAULT_BED_COLOR,
                           'border_color': 'black',
                           'line_width': 0.5,
                           # To remove in next 1.0
                           'prefered_name': 'transcript_name',
                           'merge_transcripts': False,
                           # end to remove
                           'region': None  # Cannot be set manually but is set by tracksClass
                           }
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted']}
    # To remove in next 1.0
    BOOLEAN_PROPERTIES = ['merge_transcripts']
    # And replace by:
    # BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['file', 'file_type',
                         'overlay_previous', 'orientation',
                         'title', 'color', 'border_color',
                         # To remove in next 1.0
                         'prefered_name']
    FLOAT_PROPERTIES = {'line_width': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        GenomeTrack.__init__(self, *args, **kwarg)
        self.bed_type = None  # once the bed file is processed,
        # this is bed3, bed4, bed5, bed6, bed8, bed9 or bed12
        self.interval_tree = {}  # interval tree of the bed regions
        self.interval_tree, min_score, max_score = self.process_bed(self.properties['region'])
        if self.colormap is not None:
            if self.properties.get('min_value', None) is not None:
                min_score = self.properties['min_value']
            if self.properties.get('max_value', None) is not None:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.colormap)
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def set_properties_defaults(self):
        super(BedLikeTrack, self).set_properties_defaults()
        self.colormap = None
        self.parametersUsingColormap = []
        # check if the color given is a color map
        is_colormap = self.process_color('color', colormap_possible=True,
                                         bed_rgb_possible=True,
                                         default_value_is_colormap=False)
        if is_colormap:
            self.colormap = self.properties['color']
            self.parametersUsingColormap.append('color')

        # check if border_color and color_utr are colors
        # if they are part of self.properties
        # (for example, TADsTracks do not have color_utr)
        for param in [p for p in ['border_color', 'color_utr']
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

    def get_bed_handler(self, plot_regions=None):
        if not self.properties.get('global_max_row', False):
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
                                 self.properties['merge_transcripts'])
            total_length = bed_file_h.length
        else:
            # end of remove
            total_length = count_lines(opener(file_to_open),
                                       asBed=True)
            bed_file_h = ReadBed(opener(file_to_open))

        return(bed_file_h, total_length)

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

        if self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, 0)
        else:
            ax.set_ylim(0, ymax)

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

    def plot_y_axis(self, ax, plot_axis):
        if self.colormap is not None:
            self.plot_custom_cobar(ax, fraction=1)
