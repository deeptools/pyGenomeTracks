# -*- coding: utf-8 -*-
from __future__ import division
from . BedGraphTrack import BedGraphTrack
from . GenomeTrack import GenomeTrack
import json
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


class EpilogosTrack(BedGraphTrack):
    """
    The data format for this type of track can be found
    at http://wiki.wubrowse.org/QuantitativeCategorySeries.
    In summary, first three columns are chrom, start,end,
    Next columns are a jason formatted line which contains
    two attributes, id and qcat array. E.g.:

     id:8,qcat:[ [-0.0079,6], [-0.0056,17], [-0.0035,13], [-0.0023,5], ..

    """
    SUPPORTED_ENDINGS = ['.qcat', '.qcat.bgz']
    TRACK_TYPE = 'epilogos'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
    """.format(TRACK_TYPE)

    def set_properties_defaults(self):
        if 'max_value' not in self.properties or self.properties['max_value'] == 'auto':
            self.properties['max_value'] = None

        if 'min_value' not in self.properties or self.properties['min_value'] == 'auto':
            self.properties['min_value'] = None

        if 'type' not in self.properties:
            self.properties['type'] = 'matrix'

        if 'pos score in bin' not in self.properties:
            self.properties['pos score in bin'] = 'center'

        if 'show data range' not in self.properties:
            self.properties['show data range'] = 'yes'

        if 'plot horizontal lines' not in self.properties:
            self.properties['plot horizontal lines'] = 'no'

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        """
        Plots a bedgraph matrix file, that instead of having
        a single value per bin, it has several values.
        """
        from matplotlib import cm
        cmap = cm.get_cmap('tab20b')

        values_list, pos_list = self.get_scores(chrom_region, start_region, end_region)
        ymin = float('Inf')
        ymax = -float('Inf')

        edgecolor = 'black'
        if len(values_list) > 1000:
            edgecolor = 'none'
            linewidth = 0
        elif 1000 > len(values_list) > 500:
            linewidth = 0.01
        else:
            linewidth = 0.5

        rects = []
        for idx, qcat_json in enumerate(values_list):
            # the qcat_json is a pseudo json line, that misses
            # the { } and the quotes. The following lines fix that.
            qcat_json = qcat_json[0]
            qcat_json = '{' + qcat_json.replace('id', '"id"').replace('qcat', '"qcat"') + '}'
            qcat = json.loads(qcat_json)

            # find the minimum negative value
            neg_values = [x[0] for x in qcat['qcat'] if x[0] < 0]
            if len(neg_values) > 0:
                min_neg_sum = sum([x[0] for x in qcat['qcat'] if x[0] < 0])
            else:
                min_neg_sum = 0

            # Draw a rectangle for each value
            y_low = min_neg_sum
            for qcat_value, qcat_id in qcat['qcat']:
                if qcat_value == 0:
                    continue
                height = abs(qcat_value)
                if height + y_low > ymax:
                    ymax = height + y_low
                if y_low < ymin:
                    ymin = y_low
                start, end = pos_list[idx]

                # Rectangle(xy, width, height, angle=0.0, **kwargs)
                # Draw a rectangle with lower left at xy = (x, y) with specified width, height and rotation angle.
                rects.append(Rectangle((start, y_low), end - start, height,
                                       edgecolor=edgecolor, facecolor=cmap(qcat_id / 15), linewidth=linewidth))
                y_low += height
        collection = PatchCollection(rects, match_original=True)
        ax.add_collection(collection)
        ax.set_ylim(ymin - ymin * 0.01, ymax + ymax * 0.01)
        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center', transform=label_ax.transAxes)
