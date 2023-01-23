from . GenomeTrack import GenomeTrack
from .. utilities import file_to_intervaltree, plot_coverage, InputError, transform, change_chrom_names
import numpy as np
import pyBigWig
import tempfile
import os
import pysam

DEFAULT_BEDGRAPH_COLOR = '#a6cee3'


class BedGraphTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.bg', '.bg.gz', '.bg.bgz',
                         '.bedgraph', '.bedgraph.gz', '.bedgraph.bgz',
                         '.bedGraph', '.bedGraph.gz', '.bedGraph.bgz',
                         '.bdg', '.bdg.gz', '.bdg.bgz']
    TRACK_TYPE = 'bedgraph'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
color = green
# To use a different color for negative values
#negative_color = red
# To use transparency, you can use alpha
# default is 1
# alpha = 0.5
# the default for min_value and max_value is 'auto' which means that the scale will go
# roughly from the minimum value found in the region plotted to the maximum value found.
min_value = 0
#max_value = auto
# to convert missing data (NaNs) into zeros. Otherwise, missing data is not plotted.
nans_to_zeros = true
# for type, the options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similarly points:ms sets the point size (markersize (ms) to the given float
# type = line:0.5
# type = points:0.5
# If you want to plot a 4C track where you want to link
# the non-missing data (NaNs) together and only use the
# middle of the region instead of the region itself:
# Default is false.
# use_middle = true
# By default the bedgraph is plotted at the base pair
# resolution. This can lead to very large pdf/svg files
# If plotting large regions.
# If you want to decrase the size of your file.
# You can either rasterize the bedgraph profile by using:
# rasterize = true
# Or use a summary method on a given number of bin:
# The possible summary methods are given by pyBigWig:
# mean/average/stdev/dev/max/min/cov/coverage/sum
# summary_method = mean
# number_of_bins = 700
# set show_data_range to false to hide the text on the left showing the data range
show_data_range = true
# to compute operations on the fly on the file
# or between 2 bedgraph files
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
# '-log': transformed_values = log(log_pseudocount + initial_values)
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
                           'color': DEFAULT_BEDGRAPH_COLOR,
                           'negative_color': None,
                           'alpha': 1,
                           'nans_to_zeros': False,
                           'use_middle': False,
                           'summary_method': None,
                           'rasterize': False,
                           'number_of_bins': 700,
                           'type': 'fill',
                           'region': None,  # Cannot be set manually but is set by tracksClass
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
                                              'cov', 'sum', None],
                           'transform': ['no', 'log', 'log1p', '-log', 'log2',
                                         'log10'],
                           'y_axis_values': ['original', 'transformed']}
    BOOLEAN_PROPERTIES = ['show_data_range', 'nans_to_zeros',
                          'use_middle', 'rasterize', 'grid']
    STRING_PROPERTIES = ['file', 'file_type', 'overlay_previous',
                         'orientation', 'summary_method',
                         'title', 'color', 'negative_color',
                         'type', 'transform', 'y_axis_values',
                         'second_file', 'operation']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'log_pseudocount': [- np.inf, np.inf],
                        'alpha': [0, 1],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {'number_of_bins': [1, np.inf]}
    # The color can only be a color
    # negative_color can only be a color or None

    def __init__(self, properties_dict):
        super(BedGraphTrack, self).__init__(properties_dict)
        self.tbx2 = None
        self.load_file()

        self.interval_tree2 = None

        if 'second_file' in self.properties['operation']:
            if self.properties['second_file'] is None:
                raise InputError(f"operation: {self.properties['operation']}"
                                 " requires to set the parameter"
                                 " second_file.")
            else:
                # First try to open it as a Tabix file
                try:
                    # from the tabix file is not possible to know the
                    # global min and max
                    self.tbx2 = pysam.TabixFile(self.properties['second_file'])
                except IOError:
                    # load the file as an interval tree
                    self.interval_tree2, __, __ = file_to_intervaltree(self.properties['second_file'],
                                                                       self.properties['region'])

    def set_properties_defaults(self):
        super(BedGraphTrack, self).set_properties_defaults()
        super(BedGraphTrack, self).process_type_for_coverage_track()
        self.process_color('color')
        if self.properties['negative_color'] is None:
            self.properties['negative_color'] = self.properties['color']
        else:
            self.process_color('negative_color')

        if 'second_file' in self.properties['operation'] and \
           self.properties['second_file'] is not None and \
           self.properties['summary_method'] is None and \
           not self.properties['use_middle']:
            self.log.warning("When an operation is computed"
                             " between 2 files"
                             " a summary_method needs to be"
                             " used. Will use mean.\n")
            self.properties['summary_method'] = 'mean'

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

    def load_file(self):
        self.tbx = None
        # try to load a tabix file if available
        try:
            # from the tabix file is not possible to know the
            # global min and max
            self.tbx = pysam.TabixFile(self.properties['file'])
        except IOError:
            # load the file as an interval tree
            self.interval_tree, __, __ = file_to_intervaltree(self.properties['file'],
                                                              self.properties['region'])

        self.num_fields = None

    def _get_row_data(self, row, tbx_var='self.tbx'):
        """
        Returns the chrom, start, end and fields from either a tabix or a
        interval tree.
        Args:
            row: if tabix, the row comes from self.tbx.fetch otherwise
            comes from sorted(interval_tree[chrom] ...

        Returns:
            start, end, fields where values is a list

        """
        tbx = eval(tbx_var)
        if tbx is not None:
            fields = row.split("\t")
            values = fields[3:]
            start = int(fields[1])
            end = int(fields[2])

        else:
            values = row.data
            start = row.begin
            end = row.end

        # set the num_fields value
        # it is expected that the number of fields per row
        # is equal. This value is used for regions not covered
        # in the file and that should be represented as nans
        if self.num_fields is None:
            self.num_fields = len(values)
        return start, end, values

    def get_scores(self, chrom_region, start_region, end_region,
                   return_nans=True, tbx_var='self.tbx', inttree_var='self.interval_tree'):
        """
        Retrieves the score (or scores or whatever fields are in a bedgraph like file) and the positions
        for a given region. If return_nans is True the pos_list goes until at least end_region.
        In case there is no item in the region. It returns [], []
        Args:
            chrom_region:
            start_region:
            end_region:
        Returns:
            tuple:
                scores_list, pos_list
        """
        score_list = []
        pos_list = []
        tbx = eval(tbx_var)
        if tbx is not None:
            if chrom_region not in tbx.contigs:
                chrom_region_before = chrom_region
                chrom_region = change_chrom_names(chrom_region)
                if chrom_region not in tbx.contigs:
                    self.log.warning("*Warning*\nNeither "
                                     + chrom_region_before + " nor "
                                     + chrom_region + " exists as a "
                                     "chromosome name inside the bedgraph "
                                     "file. This will generate an empty "
                                     "track!!\n")
                    return score_list, pos_list

            iterator = tbx.fetch(chrom_region, start_region, end_region)

        else:
            inttree = eval(inttree_var)
            if chrom_region not in list(inttree):
                chrom_region_before = chrom_region
                chrom_region = change_chrom_names(chrom_region)
                if chrom_region not in list(inttree):
                    self.log.warning("*Warning*\nNo interval was found when "
                                     "overlapping with both "
                                     f"{chrom_region_before}:{start_region}-{end_region}"
                                     f" and {chrom_region}:{start_region}-{end_region}"
                                     " inside the bedgraph file. "
                                     "This will generate an empty "
                                     "track!!\n")
                    return score_list, pos_list
            iterator = iter(sorted(inttree[chrom_region][start_region - 10000:end_region + 10000]))

        prev_end = start_region
        for row in iterator:
            start, end, values = self._get_row_data(row, tbx_var)
            # if the region is not consecutive with respect to the previous
            # nan values are added.
            if return_nans and prev_end < start:
                score_list.append(np.repeat(np.nan, self.num_fields))
                pos_list.append((prev_end, start))
            prev_end = end
            score_list.append(values)
            pos_list.append((start, end))

        if self.num_fields is None:
            # This means there is no value:
            return score_list, pos_list

        # Add a last value if needed:
        if prev_end < end_region and return_nans:
            score_list.append(np.repeat(np.nan, self.num_fields))
            pos_list.append((prev_end, end_region))

        return score_list, pos_list

    def plot(self, ax, chrom_region, start_region, end_region):
        score_list, pos_list = self.get_scores(chrom_region, start_region, end_region)
        if pos_list == []:
            self.adjust_ylim(ax)
            return
        try:
            score_list = [float(x[0]) for x in score_list]
        except ValueError as ve:
            if "could not convert string to float: 'NA'" in str(ve):
                self.log.warning("*Warning*\nNA were found in the bedgraph"
                                 " will be replaced by nan.\n")
                score_list = [float(x[0]) if x[0] != 'NA' else float('nan') for x in score_list]
            else:
                raise ve
        if self.properties['use_middle']:
            x_values = np.asarray([(t[0] + t[1]) / 2
                                   for i, t in enumerate(pos_list)
                                   if not np.isnan(score_list[i])],
                                  dtype=float)
            score_list = np.asarray([x for x in score_list if not np.isnan(x)],
                                    dtype=float)
        elif self.properties['summary_method'] is not None:
            score_list, x_values = self.get_values_as_bigwig(score_list,
                                                             pos_list,
                                                             chrom_region,
                                                             start_region,
                                                             end_region)
        else:
            score_list, x_values = self.get_values_as_bdg(score_list,
                                                          pos_list)
        # compute the operation
        operation = self.properties['operation']
        # Substitute log by np.log to make it evaluable:
        operation = operation.replace('log', 'np.log')
        if operation == 'file':
            pass
        elif 'second_file' not in operation:
            try:
                new_score_list = eval('[' + operation + ' for file in score_list]')
                new_score_list = np.array(new_score_list)
            except Exception as e:
                raise Exception("The operation in section "
                                f"{self.properties['section_name']} could not "
                                f"be computed: {e}")
            else:
                score_list = new_score_list

        else:
            score_list2, pos_list2 = self.get_scores(chrom_region, start_region, end_region,
                                                     tbx_var='self.tbx2',
                                                     inttree_var='self.interval_tree2')
            if pos_list2 == []:
                self.adjust_ylim(ax)
                return
            try:
                score_list2 = [float(x[0]) for x in score_list2]
            except ValueError as ve:
                if "could not convert string to float: 'NA'" in str(ve):
                    self.log.warning("*Warning*\nNA were found in the bedgraph"
                                     " will be replaced by nan.\n")
                    score_list2 = [float(x[0]) if x[0] != 'NA' else float('nan') for x in score_list2]
                else:
                    raise ve
            if self.properties['use_middle']:
                x_values2 = np.asarray([(t[0] + t[1]) / 2
                                        for i, t in enumerate(pos_list2)
                                        if not np.isnan(score_list2[i])],
                                       dtype=float)
                score_list2 = np.asarray([x for x in score_list2 if not np.isnan(x)],
                                         dtype=float)
                if not all([x1 == x2 for x1, x2 in zip(x_values, x_values2)]):
                    # The x are not compatible we need to extrapolate:
                    new_x = sorted(np.unique(np.concatenate((x_values, x_values2), axis=0)))
                    new_score_list = []
                    i = 0
                    new_score_list2 = []
                    i2 = 0
                    for x in new_x:
                        if x_values[i] == x:
                            new_score_list.append(score_list[i])
                            i += 1
                        elif i == 0:
                            new_score_list.append(np.nan)
                        elif i == len(score_list):
                            new_score_list.append(np.nan)
                        else:
                            y1 = score_list[i - 1]
                            y2 = score_list[i]
                            x1 = x_values[i - 1]
                            x2 = x_values[i]
                            new_score_list.append(y1 + (y2 - y1) * (x - x1) / (x2 - x1))
                        if x_values2[i2] == x:
                            new_score_list2.append(score_list2[i2])
                            i2 += 1
                        elif i2 == 0:
                            new_score_list.append(np.nan)
                        elif i2 == len(score_list):
                            new_score_list.append(np.nan)
                        else:
                            y1 = score_list2[i2 - 1]
                            y2 = score_list2[i2]
                            x1 = x_values2[i2 - 1]
                            x2 = x_values2[i2]
                            new_score_list2.append(y1 + (y2 - y1) * (x - x1) / (x2 - x1))
                    x_values = np.array(new_x)
                    score_list = new_score_list
                    score_list2 = new_score_list2
            else:
                score_list2, x_values2 = self.get_values_as_bigwig(score_list2,
                                                                   pos_list2,
                                                                   chrom_region,
                                                                   start_region,
                                                                   end_region)
            # compute the operation
            try:
                new_score_list = eval('[' + operation + ' for file,second_file in zip(score_list, score_list2)]')
                new_score_list = np.array(new_score_list)
            except Exception as e:
                raise Exception("The operation in section "
                                f"{self.properties['section_name']} could not"
                                f" be computed: {e}")
            else:
                score_list = new_score_list

        transformed_scores = transform(score_list,
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

        if self.properties['rasterize']:
            ax.set_rasterized(True)

    def get_values_as_bigwig(self, score_list, pos_list, chrom_region,
                             start_region, end_region):
        # A temporary file is created
        id, temp_bigwig_file = tempfile.mkstemp(suffix='.bw')
        # We write into it
        bw = pyBigWig.open(temp_bigwig_file, 'w')
        # The list of chromosome size is the name of the current chromosome to
        # the last position of the pos_list
        bw.addHeader([(chrom_region, pos_list[-1][1])])
        # The starts, ends, score are stored
        bw.addEntries([chrom_region] * len(pos_list),
                      [p[0] for p in pos_list],
                      ends=[p[1] for p in pos_list],
                      values=score_list)
        bw.close()
        # The temporary file is opened
        bw = pyBigWig.open(temp_bigwig_file)
        # The scores are the summary:
        scores_per_bin = np.array(bw.stats(chrom_region,
                                           start_region,
                                           end_region,
                                           nBins=self.properties['number_of_bins'],
                                           type=self.properties['summary_method'])).astype(float)
        os.remove(temp_bigwig_file)
        if self.properties['nans_to_zeros'] and np.any(np.isnan(scores_per_bin)):
            scores_per_bin[np.isnan(scores_per_bin)] = 0
        x_values = np.linspace(start_region, end_region,
                               self.properties['number_of_bins'])

        return scores_per_bin, x_values

    def get_values_as_bdg(self, score_list, pos_list):
        # the following two lines will convert the score_list and the
        # tuples in pos list (where is item is tuple(start, end)
        # into an x value (pos_list) and a y value (score_list)
        # where x = start1, end1, star2, end2 ...
        # and y = score1, score1, score2, score2 ...

        # convert [1, 2, 3 ...] in [1, 1, 2, 2, 3, 3 ...]
        score_list = np.repeat(score_list, 2)
        # convert [(0, 10), (10, 20), (20, 30)] into [0, 10, 10, 20, 20, 30]
        x_values = np.asarray(sum(pos_list, tuple()), dtype=float)

        if self.properties['nans_to_zeros']:
            score_list[np.isnan(score_list)] = 0

        return score_list, x_values

    def plot_y_axis(self, ax, plot_axis):
        super(BedGraphTrack, self).plot_y_axis(ax, plot_axis,
                                               self.properties['transform'],
                                               self.properties['log_pseudocount'],
                                               self.properties['y_axis_values'],
                                               self.properties['grid'])

    def __del__(self):
        if self.tbx is not None:
            self.tbx.close()
        if self.tbx2 is not None:
            self.tbx2.close()
