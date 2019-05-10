from . GenomeTrack import GenomeTrack
import numpy as np

DEFAULT_BIGWIG_COLOR = '#33a02c'


class BigWigTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.bw', '.bigwig']
    TRACK_TYPE = 'bigwig'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
color = #666666
# the default for min_value and max_value is 'auto' which means that the scale will go
# from the minimum value found in the region plotted to the maximum value found.
min_value = 0
#max_value = auto
# The number of bins takes the region to be plotted and divides it into the number of bins specified
# Then, at each bin the bigwig mean value is computed and plotted.
# A lower number of bins produces a coarser tracks
number of bins = 500
# to convert missing data (NaNs) into zeros. Otherwise, missing data is not plotted.
nans to zeros = True
# the summary method by default is mean. Other
# methods are min and max
summary method = mean
# for type, the options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similarly points:ms sets the point size (markersize (ms) to the given float
# type = line:0.5
# type = points:0.5
# set show data range to no to hide the text on the upper-left showing the data range
show data range = yes
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        import pyBigWig
        self.bw = pyBigWig.open(self.properties['file'])
        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BIGWIG_COLOR

        if 'negative color' not in self.properties:
            self.properties['negative color'] = self.properties['color']

        if 'summary method' not in self.properties:
            self.properties['summary method'] = 'mean'

        if 'nans to zeros' not in self.properties:
            self.properties['nans to zeros'] = False

        self.plot_type = 'fill'
        self.size = None

        if 'type' in self.properties:
            if self.properties['type'].find(":") > 0:
                self.plot_type, size = self.properties['type'].split(":")
                try:
                    self.size = float(size)
                except ValueError:
                    exit("Invalid value: 'type = {}' in section: {}\n"
                         "A number was expected and found '{}'".format(self.properties['type'],
                                                                       self.properties['section_name'],
                                                                       size))
            else:
                self.plot_type = self.properties['type']

        if self.plot_type not in ['line', 'points', 'fill']:
            exit("Invalid: 'type = {}' in section: {}\n".format(self.properties['type'],
                                                                self.properties['section_name']))

    def plot(self, ax, chrom_region, start_region, end_region):
        formated_region = "{}:{}-{}".format(chrom_region, start_region, end_region)

        num_bins = 700
        if 'number of bins' in self.properties:
            try:
                num_bins = int(self.properties['number of bins'])
            except TypeError:
                num_bins = 700
                self.log.warn("'number of bins' value: {} for bigwig file {} "
                              "is not valid. Using default value (700)".format(self.properties['number of bins'],
                                                                               self.properties['file']))

        if chrom_region not in self.bw.chroms().keys():
            chrom_region_before = chrom_region
            chrom_region = self.change_chrom_names(chrom_region)
            if chrom_region not in self.bw.chroms().keys():
                self.log.error("*Error*\nNeither " + chrom_region_before + " "
                               "nor " + chrom_region + " exits as a chromosome"
                               " name inside the bigwig file.\n")
                return

        chrom_region = self.check_chrom_str_bytes(self.bw.chroms().keys(), chrom_region)

        if chrom_region not in self.bw.chroms().keys():
            self.log.warn("Can not read region {} from bigwig file:\n\n"
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
                                                        end_region, nBins=num_bins,
                                                        type=self.properties['summary method'])).astype(float)
                if self.properties['nans to zeros'] and np.any(np.isnan(scores_per_bin)):
                    scores_per_bin[np.isnan(scores_per_bin)] = 0
            except Exception as e:
                import pyBigWig
                self.bw = pyBigWig.open(self.properties['file'])

                self.log.warning("error found while reading bigwig scores ({}).\nTrying again. Iter num: {}".
                                 format(e, num_tries))
                pass
            else:
                if num_tries > 1:
                    self.log.warning("After {} the scores could be computed".format(num_tries))
                break

        x_values = np.linspace(start_region, end_region, num_bins)
        if self.plot_type == 'line':
            if self.properties['color'] == self.properties['negative color']:
                ax.plot(x_values, scores_per_bin, '-', linewidth=self.size, color=self.properties['color'])
            else:
                import warnings
                warnings.warn('Line plots with a different negative color might not look pretty')
                pos_x_values = x_values.copy()
                pos_x_values[scores_per_bin < 0] = np.nan
                ax.plot(pos_x_values, scores_per_bin, '-', linewidth=self.size, color=self.properties['color'])

                neg_x_values = x_values.copy()
                neg_x_values[scores_per_bin >= 0] = np.nan
                ax.plot(neg_x_values, scores_per_bin, '-', linewidth=self.size, color=self.properties['negative color'])

        elif self.plot_type == 'points':
            ax.plot(x_values[scores_per_bin >= 0], scores_per_bin[scores_per_bin >= 0], '.',
                    markersize=self.size, color=self.properties['color'])
            ax.plot(x_values[scores_per_bin < 0], scores_per_bin[scores_per_bin < 0], '.',
                    markersize=self.size, color=self.properties['negative color'])

        else:
            ax.fill_between(x_values, scores_per_bin, linewidth=0.1, color=self.properties['color'],
                            facecolor=self.properties['color'], where=scores_per_bin >= 0, interpolate=True)
            ax.fill_between(x_values, scores_per_bin, linewidth=0.1, color=self.properties['negative color'],
                            facecolor=self.properties['negative color'], where=scores_per_bin < 0, interpolate=True)

        ymin, ymax = ax.get_ylim()
        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            ymax = self.properties['max_value']
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            ymin = self.properties['min_value']

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, ymin)
        else:
            ax.set_ylim(ymin, ymax)

        return ax
