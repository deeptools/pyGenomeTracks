# -*- coding: utf-8 -*-
from . GenomeTrack import GenomeTrack


class TextTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.txt']  # this is used by make_tracks_file to guess the type of track based on file name
    TRACK_TYPE = 'text'
    OPTIONS_TXT = """
height = 3
title =
text =
# x position of text in the plot (in bp)
x position =
"""
    DEFAULTS_PROPERTIES = {'text': 'hello world'}
    NECESSARY_PROPERTIES = ['x_position']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = []
    STRING_PROPERTIES = ['text', 'title', 'file_type']
    FLOAT_PROPERTIES = {'height': [0, np.inf],
                        'x_position': [0, np.inf]}
    INTEGER_PROPERTIES = {}
    def plot(self, ax, chrom, region_start, region_end):
        """
        This example simply plots the given title at a fixed
        location in the axis. The chrom, region_start and region_end
        variables are not used.
        Args:
            ax: matplotlib axis to plot
            chrom_region: chromosome name
            start_region: start coordinate of genomic position
            end_region: end coordinate
        """
        # print text at position x = self.properties['x position'] and y = 0.5 (center of the plot)
        ax.text(self.properties['x_position'], 0.5, self.properties['text'])
