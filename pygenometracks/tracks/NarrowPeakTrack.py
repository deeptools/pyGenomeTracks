from . GenomeTrack import GenomeTrack
from . BedGraphTrack import BedGraphTrack
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


class NarrowPeakTrack(BedGraphTrack):
    SUPPORTED_ENDINGS = ['.narrowPeak']
    TRACK_TYPE = 'narrow_peak'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
#max_value = 0.70
show data range = yes
show labels = yes
# the narrowPeak format provides the information of the
# peak summit. By default this information is used
# although some peaks may look crooked.
use summit = yes
# type of plot: either box or peak
# box will plot a rectangle of the peak width
# peak will plot the shape of the peak, whose height is the
# narrowPeak file signal value (usually peak coverage)
type = peak
# if the peaks look too thin, the can be adjusted
width adjust = 1.5
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, properties_dict):
        super(NarrowPeakTrack, self).__init__(properties_dict)
        self.patches = []

    def set_properties_defaults(self):
        if 'color' not in self.properties:
            self.properties['color'] = '#FF000080'  # red, alpha=0.55
        if 'show data range' not in self.properties:
            self.properties['show data range'] = 'yes'
        if 'show labels' not in self.properties:
            self.properties['show labels'] = 'yes'
        if 'use summits' not in self.properties:
            self.properties['use summits'] = 'yes'
        if 'width adjust' not in self.properties:
            self.properties['width adjust'] = 1.5
        else:
            self.properties['width adjust'] = float(self.properties['width adjust'])
        if 'type' not in self.properties:
            self.properties['type'] = 'peak'

    def peak_plot(self, start, end, height, center=None, width_adjust=1.5):
        # uses bezier curves to plot a shape that
        # looks like a peak
        from matplotlib.path import Path
        import matplotlib.patches as patches
        peak_width = float(end - start)
        if center is None:
            center = peak_width / 2 + start
        if width_adjust != 1:
            start -= width_adjust * peak_width / 2
            end += width_adjust * peak_width / 2
            peak_width *= width_adjust

        path_data = [
            (Path.MOVETO, (start, 0)),
            (Path.CURVE4, (start + peak_width / 2, 0)),
            (Path.CURVE4, (start + peak_width * 0.4, height)),
            (Path.CURVE4, (center, height)),
            (Path.CURVE4, (end - peak_width * 0.4, height)),
            (Path.CURVE4, (end - peak_width / 2, 0)),
            (Path.CURVE4, (end, 0))]

        codes, verts = zip(*path_data)
        path = Path(verts, codes)
        return patches.PathPatch(path)

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        score_list, pos_list = self.get_scores(chrom_region, start_region, end_region, return_nans=False)

        self.patches = []

        max_signal = -1
        for idx, peak in enumerate(score_list):
            name, score, strand, signal_value, p_value, q_value, summit = peak

            signal_value = float(signal_value)
            p_value = float(p_value)
            q_value = float(q_value)
            summit = int(summit)
            start, end = pos_list[idx]
            if summit > 0:
                summit = start + summit
            else:
                summit = None
            if self.properties['type'] == 'box':
                self.patches.append(Rectangle((start, 0), end - start, 100, edgecolor='black',))
                max_signal = 110
            else:
                if signal_value > max_signal:
                    max_signal = signal_value
                self.patches.append(self.peak_plot(start, end, signal_value, center=summit,
                                                   width_adjust=self.properties['width adjust']))

            x_pos = start + float(end - start) / 2
            y_pos = 0 - max_signal * 0.05
            if self.properties['show labels'] != 'no':
                ax.text(x_pos, y_pos, "{}\np-val:{:.1f}\nq-val:{:.1f}".format(name, p_value, q_value),
                        horizontalalignment='center', size='smaller', verticalalignment='top')

        collection = PatchCollection(self.patches, facecolor=self.properties['color'])
        ax.add_collection(collection)

        if 'max_value' not in self.properties or self.properties['max_value'] == 'auto':
            self.properties['max_value'] = max_signal

        ymax = self.properties['max_value']
        if self.properties['show labels'] != 'no':
            if self.properties['type'] == 'box':
                ymin = ymax * -3
            else:
                ymin = ymax * -0.8
        else:
            ymin = 0
        ax.set_ylim(ymin, ymax)
        if self.properties['show data range'] != 'no' and self.properties['type'] != 'box':
            if float(ymax) % 1 == 0:
                ymax_print = int(ymax)
            else:
                ymax_print = "{:.1f}".format(ymax)
            ymin_print = 0
            small_x = 0.01 * (end_region - start_region)
            # by default show the data range
            ax.text(start_region - small_x, ymax - ymax * 0.01,
                    "[{}-{}]".format(ymin_print, ymax_print), horizontalalignment='left', verticalalignment='top')

        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')
