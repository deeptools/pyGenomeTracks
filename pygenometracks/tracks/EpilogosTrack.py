# -*- coding: utf-8 -*-
from __future__ import division
from . BedGraphTrack import BedGraphTrack
from . GenomeTrack import GenomeTrack
import json
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
import numpy as np


class EpilogosTrack(BedGraphTrack):
    """
    The data format for this type of track can be found
    at http://wiki.wubrowse.org/QuantitativeCategorySeries.
    In summary, first three columns are chrom, start,end,
    Next columns are a json formatted line which contains
    two attributes, id and qcat array. E.g.:

     id:8,qcat:[ [-0.0079,6], [-0.0056,17], [-0.0035,13], [-0.0023,5], ..

    """
    SUPPORTED_ENDINGS = ['.qcat', '.qcat.bgz']
    TRACK_TYPE = 'epilogos'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# The categories file should contain the color information for each category id
# A categories file should look like:
# {{
# "categories":{{
#           "1":["Active TSS","#ff0000"],
#           "2":["Flanking Active TSS","#ff4500"],
#           "3":["Transcr at gene 5\" and 3\"","#32cd32"],
#           "4":["Strong transcription","#008000"],
#           "5":["Weak transcription","#006400"]
# 	}}
#}}
categories_file = <path to json categories file>
# optional. If not given, it is guessed from the file ending.
file_type = {TRACK_TYPE}
    """
    DEFAULTS_PROPERTIES = {'categories_file': None,
                           'orientation': None,
                           'region': None}  # Cannot be set manually but is set by tracksClass
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted']}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['categories_file', 'file', 'file_type',
                         'overlay_previous', 'orientation',
                         'title']
    FLOAT_PROPERTIES = {'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, properties_dict):
        GenomeTrack.__init__(self, properties_dict)

        self.load_file()

    def set_properties_defaults(self):
        GenomeTrack.set_properties_defaults(self)
        # load categories file
        if self.properties['categories_file'] is not None:
            with open(self.properties['categories_file']) as f:
                self.categories = json.load(f)['categories']
        else:
            self.categories = None

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plots a bedgraph matrix file, that instead of having
        a single value per bin, it has several values.
        """
        cmap = cm.get_cmap('tab20b')

        values_list, pos_list = self.get_scores(chrom_region, start_region, end_region)
        if pos_list == []:
            return

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
            if not isinstance(qcat_json, str):
                # This would happen if the qcat file has a missing value.
                # The missing value is filled with np.repeat(np.nan, ..) and should be skipped here.
                continue
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
                # use color from categories file is given
                if self.categories is not None:
                    qcat_color = self.categories[str(qcat_id)][1]
                else:
                    qcat_color = cmap(qcat_id / 15)
                height = abs(qcat_value)
                if height + y_low > ymax:
                    ymax = height + y_low
                if y_low < ymin:
                    ymin = y_low
                start, end = pos_list[idx]

                # Rectangle(xy, width, height, angle=0.0, **kwargs)
                # Draw a rectangle with lower left at xy = (x, y) with specified width, height and rotation angle.
                rects.append(Rectangle((start, y_low), end - start, height,
                                       edgecolor=edgecolor, facecolor=qcat_color, linewidth=linewidth))
                y_low += height
        collection = PatchCollection(rects, match_original=True)
        ax.add_collection(collection)
        ax.set_ylim(ymin - ymin * 0.01, ymax + ymax * 0.01)

        if self.properties['orientation'] == 'inverted':
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(ymax, ymin)

    def plot_y_axis(self, ax, plot_axis):
        GenomeTrack.plot_y_axis(self, ax, plot_axis)

    def __del__(self):
        if self.tbx is not None:
            self.tbx.close()
