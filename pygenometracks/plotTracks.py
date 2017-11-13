

"""
[spacer]
# This will add a space between two tracks. No options are required, but
# for finer control, the amount of space required can be given
# height value in centimeters
height = 1.5

[bigwig]
file = /data/manke/group/ramirez/HiC-ThomasLing/data/external/Graveley_mRNA-seq/GSM390060_Kc167-4_spa.bw
title = Kc RNA-seq (Cherbas et al.)
color = black
min_value = 0
#max_value = auto
height = 1.5
number of bins = 500
nans to zeros = True
# options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similary points:ms sets the point size (markersize (ms) to the given float
type = line
# type = line:0.5
# type = points:0.5
#optional in case it can not be guessed by the file ending
file_type = bigwig
# the orientation of the y axis can be inverted.
# useful for plotting
orientation = inverted


[simple bed]
file = file.bed
title = peaks
color = red
# optional border color. Set to none for no border color
border_color = black
height = 0.5
# optional. If not given it is guessed from the file ending (file has to end in .bed)
file_type = bed
# optional display. The default display is to accommodate as many rows as needed to avoid overlapping of the regions
# The 'collapsed' option is useful for chromatin segmentation bed files
# that should be displayed in one line. Usually, this bed files encode the color of the region and thus
# they can be distinguished.
# display = collapsed
# the interlaced option uses two rows to plot the regions.
# display = interlaced

[bed genes]
# example of a genes track
# has the same options as a simple
# bed but if the type=genes is given
# the the file is interpreted as gene
# file. If the bed file contains the exon
# structure then this is plotted. Otherwise
# a region **with direction** is plotted.
file = genes.bed
title = genes
color = darkblue
#if color is a valid colormap, then the score is mapped
# to the colormap
#color = RdBlGn
height = 5
# to turn off/on printing of labels
labels = off
# optional. If not given is guessed from the file ending
file_type = bed
# optional: font size can be given if default are not good
fontsize = 10
# optional
# display = collapsed
# display = interleaved
#optional, default is black
#border_color = black
# style to plot the genes when they have exon information
#style = UCSC
#style = flybase
# maximum number of gene rows to be plotted. This
# field is useful to limit large number of close genes
# to be printed over many rows. When several images want
# to be combined this must be set to get equal size
# genes in all images
#gene rows = 10
# by default the ymax is the number of
# rows occupied by the genes in the region plotted. However,
# by setting this option, the global maximum is used instead.
# This is useful to combine images that are all consistent and
# have the same number of rows.
#global max row = yes


[chrom states]
# this is a case of a bed file that is plotted 'collapsed'
# useful to plot chromatin states if the bed file contains
# the color to plot
file = chromatinStates.bed
title = chromatin states
# color is replaced by the color in the bed file
# in this case
color = black
# optional boder color. Set to none for no border color
border_color = black
# default behaviour when plotting intervals from a
# bed file is to 'expand' them such that they
# do not overlap. The display = collapsed
# directive overlaps the intervals.
display = collapsed
height = 0.3

[bedgraph]
file = file.bg
title = bedgraph track
color = green
height = 0.2
# optional, otherwise guseed from file ending
file_type = bedgraph

[arcs]
# an arc connecting two regions can be plotted
# the file format is
#   chr1 start1 end1 chr2 start2 end2 score
# for example:
#   chr1 100 200 chr1 250 300 0.5
title =  arcs
color = red
# orientation = inverted
# if line width is not given, the score is used to set the line width
#line width = 0.5
file = arcs.txt

[vlines]
# vertical dotted lines from the top to the bottom of the figure
# can be drawn. For this a bed file is required
# but only the first two columns (chromosome name and start
# are used.
# vlines can also be given at the command line as a list
# of genomic positions. However, sometimes to give a file
# is more convenient in case many lines want to be plotted.
file = regions.bed
type = vlines


"""

from __future__ import division
import sys
import argparse
import matplotlib
matplotlib.use('Agg')

import pygenometracks.tracksClass
from pygenometracks._version import __version__

DEFAULT_BED_COLOR = '#1f78b4'
DEFAULT_BIGWIG_COLOR = '#33a02c'
DEFAULT_BEDGRAPH_COLOR = '#a6cee3'
DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
DEFAULT_TRACK_HEIGHT = 3  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
# proportion of width dedicated to (figure, legends)
# DEFAULT_WIDTH_RATIOS = (0.95, 0.05)
DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0.12, 'top': 0.9}


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description='Plots genomic tracks. The tracks.ini file can be genarated using the `make_tracks_file` '
                    'program.',
        usage="%(prog)s --tracks tracks.ini --region chr1:1000000-4000000 -o image.png")

    parser.add_argument('--tracks',
                        help='File containing the instructions to plot the tracks ',
                        type=argparse.FileType('r'),
                        required=True,
                        )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--region',
                       help='Region to plot, the format is chr:start-end')

    group.add_argument('--BED',
                       help='Instead of a region, a file containing the regions to plot, in BED format, '
                       'can be given. If this is the case, multiple files will be created using a prefix '
                       'the value of --outFileName',
                       type=argparse.FileType('r')
                       )

    parser.add_argument('--width',
                        help='figure width in centimeters',
                        type=float,
                        default=DEFAULT_FIGURE_WIDTH)

    parser.add_argument('--height',
                        help='Figure height in centimeters. If not given, the figure height is computed '
                             'based on the heights of the tracks. If given, the track height are proportionally '
                             'scaled to match the desired figure height.',
                        type=float)

    parser.add_argument('--title', '-t',
                        help='Plot title',
                        required=False)

    parser.add_argument('--outFileName', '-out',
                        help='File name to save the image, file prefix in case multiple images '
                             'are stored',
                        required=True)

    parser.add_argument('--vlines',
                        help='Genomic cooordindates separated by space. E.g. '
                        '--vlines 150000 3000000 124838433 ',
                        type=int,
                        nargs='+'
                        )

    parser.add_argument('--fontSize',
                        help='Font size for the labels of the plot',
                        type=float,
                        )

    parser.add_argument('--dpi',
                        help='Resolution for the image in case the'
                             'ouput is a raster graphics image (e.g png, jpg)',
                        type=int,
                        default=72
                        )

    parser.add_argument('--trackLabelFraction',
                        help='By default the space dedicated to the track labels is 0.05 of the'
                             'plot width. This fraction can be changed with this parameter if needed.',
                        default=0.05,
                        type=float)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def get_region(region_string):
    """
    splits a region string into
    a chrom, start_region, end_region tuple
    The region_string format is chr:start-end
    """
    if region_string:
        # separate the chromosome name and the location using the ':' character
        chrom, position = region_string.strip().split(":")

        # clean up the position
        for char in ",.;|!{}()":
            position = position.replace(char, '')

        position_list = position.split("-")
        try:
            region_start = int(position_list[0])
        except IndexError:
            region_start = 0
        try:
            region_end = int(position_list[1])
        except IndexError:
            region_end = 1e15  # a huge number
        if region_start < 0:
            region_start = 0
        if region_end <= region_start:
            exit("Please check that the region end is larger than the region start.\n"
                 "Values given:\nstart: {}\nend: {}\n".format(region_start, region_end))

        return chrom, region_start, region_end


def main(args=None):

    args = parse_arguments().parse_args(args)
    trp = pygenometracks.tracksClass.PlotTracks(args.tracks.name, args.width, fig_height=args.height,
                                                fontsize=args.fontSize, dpi=args.dpi,
                                                track_label_width=args.trackLabelFraction)

    if args.BED:
        count = 0
        for line in args.BED.readlines():
            count += 1
            try:
                chrom, start, end = line.strip().split('\t')[0:3]
            except ValueError:
                continue
            try:
                start, end = map(int, [start, end])
            except ValueError as detail:
                sys.stderr.write("Invalid value found at line\t{}\t. {}\n".format(line, detail))
            file_name = "{}_{}:{}-{}".format(args.outFileName, chrom, start, end)
            if end - start < 200000:
                start -= 100000
                end += 100000
            sys.stderr.write("saving {}'\n".format(file_name))
            trp.plot(file_name, chrom, start, end, title=args.title)
    else:
        region = get_region(args.region)
        trp.plot(args.outFileName, *region, title=args.title)
