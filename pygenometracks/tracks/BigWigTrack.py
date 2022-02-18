from . GenomeTrack import GenomeTrack
import numpy as np
from .. utilities import plot_coverage, InputError, transform, change_chrom_names
import pyBigWig

DEFAULT_BIGWIG_COLOR = '#33a02c'


class BigWigTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.bw', '.bigwig']
    TRACK_TYPE = 'bigwig'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
color = #666666
# To use a different color for negative values
#negative_color = red
# To use transparency, you can use alpha
# default is 1
# alpha = 0.5
# the default for min_value and max_value is 'auto' which means that the scale will go
# roughly from the minimum value found in the region plotted to the maximum value found.
min_value = 0
#max_value = auto
# The number of bins takes the region to be plotted and divides it
# into the number of bins specified
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
# set show_data_range to false to hide the text on the left showing the data range
show_data_range = true
# to compute operations on the fly on the file
# or between 2 bigwig files
# operation will be evaluated, it should contains file or
# file and second_file,
# we advice to use nans_to_zeros = true to avoid unexpected nan values
#operation = 0.89 * file
#operation = - file
#operation = file - second_file
#operation = log2((1 + file) / (1 + second_file))
#operation = max(file, second_file)
#second_file = path for the second file
# To log transform your data you can also use transform and log_pseudocount:
# For the transform values:
# 'log1p': transformed_values = log(1 + initial_values)
# 'log': transformed_values = log(log_pseudocount + initial_values)
# 'log2': transformed_values = log2(log_pseudocount + initial_values)
# 'log10': transformed_values = log10(log_pseudocount + initial_values)
# '-log': transformed_values = - log(log_pseudocount + initial_values)
# For example:
#tranform = log
#log_pseudocount = 2
# When a transformation is applied, by default the y axis
# gives the transformed values, if you prefer to see
# the original values:
#y_axis_values = original
# If you want to have a grid on the y-axis
#grid = true
file_type = {TRACK_TYPE}
    """

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
                           'type': 'fill',
                           'transform': 'no',
                           'log_pseudocount': 0,
                           'y_axis_values': 'transformed',
                           'second_file': None,
                           'operation': 'file',
                           'grid': False}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None}}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           'summary_method': ['mean', 'average', 'max', 'min',
                                              'stdev', 'dev', 'coverage',
                                              'cov', 'sum'],
                           'transform': ['no', 'log', 'log1p', '-log', 'log2',
                                         'log10'],
                           'y_axis_values': ['original', 'transformed']}
    BOOLEAN_PROPERTIES = ['nans_to_zeros', 'show_data_range', 'grid']
    STRING_PROPERTIES = ['file', 'file_type', 'overlay_previous',
                         'orientation', 'summary_method',
                         'title', 'color', 'negative_color',
                         'transform', 'y_axis_values',
                         'type', 'second_file', 'operation']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'log_pseudocount': [- np.inf, np.inf],
                        'alpha': [0, 1],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {'number_of_bins': [1, np.inf]}
    # The color can only be a color
    # negative_color can only be a color or None

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.bw = pyBigWig.open(self.properties['file'])
        self.bw2 = None
        if 'second_file' in self.properties['operation']:
            if self.properties['second_file'] is None:
                raise InputError(f"operation: {self.properties['operation']}"
                                 " requires to set the parameter"
                                 " second_file.")
            else:
                self.bw2 = pyBigWig.open(self.properties['second_file'])

    def set_properties_defaults(self):
        super(BigWigTrack, self).set_properties_defaults()
        super(BigWigTrack, self).process_type_for_coverage_track()
        self.process_color('color')
        if self.properties['negative_color'] is None:
            self.properties['negative_color'] = self.properties['color']
        else:
            self.process_color('negative_color')
        if self.properties['operation'] != 'file':
            self.checkoperation()
            if self.properties['transform'] != 'no':
                raise InputError("'operation' and 'transform' cannot be set at"
                                 " the same time.")
            if self.properties['y_axis_values'] == 'original':
                self.log.warning("*Warning* 'operation' is used and "
                                 "'y_axis_values' was set to 'original'. "
                                 "'y_axis_values' can only be set to "
                                 "'original' when 'transform' is used.\n"
                                 " It will be set as 'transformed'.\n")
                self.properties['y_axis_values'] = 'transformed'

    def plot(self, ax, chrom_region, start_region, end_region):

        temp_end_region, temp_nbins, scores_per_bin = self.get_scores('self.bw', self.properties['file'],
                                                                      chrom_region, start_region, end_region)
        if scores_per_bin is None:
            self.log.warning("Scores could not be computed. This will generate an empty track\n")
            return

        if self.properties['nans_to_zeros'] and np.any(np.isnan(scores_per_bin)):
            scores_per_bin[np.isnan(scores_per_bin)] = 0

        x_values = np.linspace(start_region, temp_end_region, temp_nbins)
        # compute the operation
        operation = self.properties['operation']
        # Substitute log by np.log to make it evaluable:
        operation = operation.replace('log', 'np.log')
        if operation == 'file':
            pass
        elif 'second_file' not in operation:
            try:
                new_scores_per_bin = eval('[' + operation + ' for file in scores_per_bin]')
                new_scores_per_bin = np.array(new_scores_per_bin)
            except Exception as e:
                raise Exception("The operation in section "
                                f"{self.properties['section_name']} could not "
                                f"be computed: {e}")
            else:
                scores_per_bin = new_scores_per_bin
        else:
            temp_end_region2, temp_nbins2, scores_per_bin2 = self.get_scores('self.bw2', self.properties['second_file'],
                                                                             chrom_region, start_region, end_region)
            if scores_per_bin2 is None:
                self.log.warning("Scores for second_file could not be computed. This will generate an empty track\n")
                return

            if self.properties['nans_to_zeros'] and np.any(np.isnan(scores_per_bin2)):
                scores_per_bin2[np.isnan(scores_per_bin2)] = 0

            x_values2 = np.linspace(start_region, temp_end_region2, temp_nbins2)
            if not np.all(x_values == x_values2):
                raise Exception('The two bigwig files are not compatible on this region:'
                                f'{chrom_region}:{start_region}-{end_region}')
            # compute the operation
            try:
                new_scores_per_bin = eval('[' + operation
                                          + ' for file, second_file in'
                                          ' zip(scores_per_bin,'
                                          ' scores_per_bin2)]')
                new_scores_per_bin = np.array(new_scores_per_bin)
            except Exception as e:
                raise Exception("The operation in section "
                                f"{self.properties['section_name']} could not "
                                f"be computed: {e}")
            else:
                scores_per_bin = new_scores_per_bin

        transformed_scores = transform(scores_per_bin,
                                       self.properties['transform'],
                                       self.properties['log_pseudocount'],
                                       self.properties['file'])

        plot_coverage(ax, x_values, transformed_scores, self.plot_type,
                      self.size,
                      self.properties['color'],
                      self.properties['negative_color'],
                      self.properties['alpha'],
                      self.properties['grid'])

        self.adjust_ylim(ax)

        return ax

    def plot_y_axis(self, ax, plot_axis):
        super(BigWigTrack, self).plot_y_axis(ax, plot_axis,
                                             self.properties['transform'],
                                             self.properties['log_pseudocount'],
                                             self.properties['y_axis_values'],
                                             self.properties['grid'])

    def get_scores(self, bw_var, bw_file, chrom_region, start_region, end_region):
        bw = eval(bw_var)
        scores_per_bin = None
        if chrom_region not in bw.chroms().keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in bw.chroms().keys():
                self.log.warning("*Warning*\nNeither " + chrom_region_before
                                 + " nor " + chrom_region + " exists as a "
                                 "chromosome name inside the bigwig file. "
                                 "No score will be computed for"
                                 f" {bw_file}.\n")
                scores_per_bin = np.array([np.nan] * self.properties['number_of_bins'])

        if scores_per_bin is None and start_region > bw.chroms()[chrom_region]:
            self.log.warning("*Warning*\nThe region to plot starts beyond the"
                             " chromosome size. No score will be computed for"
                             f" {bw_file}.\n"
                             f"{chrom_region} size: {bw.chroms()[chrom_region]}"
                             f". Region to plot {start_region}-{end_region}\n")
            scores_per_bin = np.array([np.nan] * self.properties['number_of_bins'])

        if scores_per_bin is None and end_region > bw.chroms()[chrom_region]:
            self.log.warning("*Warning*\nThe region to plot extends beyond the"
                             " chromosome size. Please check.\n"
                             f"{chrom_region} size: {bw.chroms()[chrom_region]}"
                             f". Region to plot {start_region}-{end_region}\n")
            temp_end_region = bw.chroms()[chrom_region]
            temp_nbins = int(self.properties['number_of_bins'] * (temp_end_region - start_region) / (end_region - start_region))
        else:
            temp_end_region = end_region
            temp_nbins = self.properties['number_of_bins']
        # on rare occasions pyBigWig may throw an error, apparently caused by a corruption
        # of the memory. This only occurs when calling trackPlot from different
        # processors. Reloading the file solves the problem.
        if scores_per_bin is None:
            num_tries = 0
            while num_tries < 5:
                num_tries += 1
                try:
                    scores_per_bin = np.array(bw.stats(chrom_region, start_region,
                                                       temp_end_region, nBins=temp_nbins,
                                                       type=self.properties['summary_method'])).astype(float)
                except Exception as e:
                    bw = pyBigWig.open(self.properties['file'])

                    self.log.warning("error found while reading bigwig scores "
                                     f"({e}).\nTrying again."
                                     f" Iter num: {num_tries}.\n")
                    pass
                else:
                    if num_tries > 1:
                        self.log.warning(f"After {num_tries} the scores could be computed.\n")
                    break
        return temp_end_region, temp_nbins, scores_per_bin

    def __del__(self):
        try:
            self.bw.close()
        except AttributeError:
            pass
        if self.bw2 is not None:
            self.bw2.close()
