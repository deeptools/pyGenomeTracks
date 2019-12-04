from . GenomeTrack import GenomeTrack
import numpy as np
from .. utilities import plot_coverage
import pyBigWig

DEFAULT_BIGWIG_COLOR = '#33a02c'


class BigWigTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.bw', '.bigwig']
    TRACK_TYPE = 'bigwig'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
color = #666666
# To use transparency, you can use alpha
# default is 1
# alpha = 0.5
# the default for min_value and max_value is 'auto' which means that the scale will go
# from the minimum value found in the region plotted to the maximum value found.
min_value = 0
#max_value = auto
# The number of bins takes the region to be plotted and divides it into the number of bins specified
# Then, at each bin the bigwig mean value is computed and plotted.
# A lower number of bins produces a coarser tracks
number_of_bins = 700
# to convert missing data (NaNs) into zeros. Otherwise, missing data is not plotted.
nans_to_zeros = true
# The possible summary methods are given by pyBigWig:
# mean/average/stdev/dev/max/min/cov/coverage/sum
# default is mean
summary_method = mean
# for type, the options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similarly points:ms sets the point size (markersize (ms) to the given float
# type = line:0.5
# type = points:0.5
# set show_data_range to false to hide the text on the upper-left showing the data range
show_data_range = true
file_type = {}
    """.format(TRACK_TYPE)

    DEFAULTS_PROPERTIES = {'max_value': None,
                           'min_value': None,
                           'show_data_range': True,
                           'orientation': None,
                           'color': DEFAULT_BIGWIG_COLOR,
                           'negative_color': None,
                           'alpha': 1,
                           'nans_to_zeros': False,
                           'summary_method': 'mean',
                           'number_of_bins': 700,
                           'type': 'fill'}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           'summary_method': ['mean', 'average', 'max', 'min',
                                              'stdev', 'dev', 'coverage',
                                              'cov', 'sum']}
    BOOLEAN_PROPERTIES = ['nans_to_zeros', 'show_data_range']
    STRING_PROPERTIES = ['file', 'file_type', 'overlay_previous',
                         'orientation', 'summary_method',
                         'title', 'color', 'negative_color',
                         'type']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'alpha': [0, 1],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {'number_of_bins': [1, np.inf]}
    # The color can only be a color
    # negative_color can only be a color or None

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.bw = pyBigWig.open(self.properties['file'])

    def set_properties_defaults(self):
        super(BigWigTrack, self).set_properties_defaults()
        super(BigWigTrack, self).process_type_for_coverage_track()
        if self.properties['negative_color'] is None:
            self.properties['negative_color'] = self.properties['color']

    def plot(self, ax, chrom_region, start_region, end_region):
        formated_region = "{}:{}-{}".format(chrom_region, start_region, end_region)

        if chrom_region not in self.bw.chroms().keys():
            chrom_region_before = chrom_region
            chrom_region = self.change_chrom_names(chrom_region)
            if chrom_region not in self.bw.chroms().keys():
                self.log.warning("*Warning*\nNeither " + chrom_region_before
                                 + " nor " + chrom_region + " existss as a "
                                 "chromosome name inside the bigwig file. "
                                 "This will generate an empty track!!\n")
                return

        chrom_region = self.check_chrom_str_bytes(self.bw.chroms().keys(), chrom_region)

        if chrom_region not in self.bw.chroms().keys():
            self.log.warning("Can not read region {} from bigwig file:\n\n"
                             "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                             "and that the region is valid".format(formated_region, self.properties['file']))

        # on rare occasions pyBigWig may throw an error, apparently caused by a corruption
        # of the memory. This only occurs when calling trackPlot from different
        # processors. Reloading the file solves the problem.
        num_tries = 0
        scores_per_bin = None
        while num_tries < 5:
            num_tries += 1
            try:
                scores_per_bin = np.array(self.bw.stats(chrom_region, start_region,
                                                        end_region, nBins=self.properties['number_of_bins'],
                                                        type=self.properties['summary_method'])).astype(float)
                if self.properties['nans_to_zeros'] and np.any(np.isnan(scores_per_bin)):
                    scores_per_bin[np.isnan(scores_per_bin)] = 0
            except Exception as e:
                self.bw = pyBigWig.open(self.properties['file'])

                self.log.warning("error found while reading bigwig scores ({}).\nTrying again. Iter num: {}".
                                 format(e, num_tries))
                pass
            else:
                if num_tries > 1:
                    self.log.warning("After {} the scores could be computed".format(num_tries))
                break

        x_values = np.linspace(start_region, end_region, self.properties['number_of_bins'])

        plot_coverage(ax, x_values, scores_per_bin, self.plot_type, self.size,
                      self.properties['color'],
                      self.properties['negative_color'],
                      self.properties['alpha'])

        ymin, ymax = ax.get_ylim()
        if self.properties['max_value'] is not None:
            ymax = self.properties['max_value']
        if self.properties['min_value'] is not None:
            ymin = self.properties['min_value']

        if self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, ymin)
        else:
            ax.set_ylim(ymin, ymax)

        return ax
