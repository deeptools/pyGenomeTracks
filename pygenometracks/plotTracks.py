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
number_of_bins = 500
nans_to_zeros = true
# options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similarly points:ms sets the point size (markersize (ms) to the given float
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
# optional border_color. Set to none for no border_color
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
# whether printing the labels
labels = false
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
#gene_rows = 10
# by default the ymax is the number of
# rows occupied by the genes in the region plotted. However,
# by setting this option, the global maximum is used instead.
# This is useful to combine images that are all consistent and
# have the same number of rows.
#global_max_row = true


[chrom states]
# this is a case of a bed file that is plotted 'collapsed'
# useful to plot chromatin states if the bed file contains
# the color to plot
file = chromatinStates.bed
title = chromatin states
# color is replaced by the color in the bed file
# in this case
color = black
# optional border_color. Set to none for no border_color
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
# if line_width is not given, the score is used to set the line width
#line_width = 0.5
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

import sys
import os
import argparse
import warnings

from pygenometracks.tracksClass import PlotTracks
from pygenometracks._version import __version__
from .utilities import InputError, get_region
import matplotlib.pyplot as plt

DEFAULT_FIGURE_WIDTH = 40  # in centimeters


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description='Plots genomic tracks on specified region(s). '
                    'Citations :\nRamirez et al.  High-resolution TADs reveal DNA '
                    'sequences underlying genome organization in flies. '
                    'Nature Communications (2018) doi:10.1038/s41467-017-02525-w\n'
                    'Lopez-Delisle et al.  pyGenomeTracks: reproducible'
                    ' plots for multivariate genomic datasets. '
                    'Bioinformatics (2020) doi:10.1093/bioinformatics/btaa692',
        usage="%(prog)s --tracks tracks.ini --region chr1:1000000-4000000 -o image.png")

    parser.add_argument('--tracks',
                        help='File containing the instructions to plot the tracks. '
                        'The tracks.ini file can be genarated using the `make_tracks_file` program.',
                        type=argparse.FileType('r'),
                        required=True,
                        )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--region',
                       help='Region to plot, the format is chr:start-end')

    group.add_argument('--BED',
                       help='Instead of a region, a file containing the regions to plot, in BED format, '
                       'can be given. If this is the case, multiple files will be created. '
                       'It will use the value of --outFileName as a template'
                       ' and put the coordinates between the file name and the extension.',
                       type=argparse.FileType('r')
                       )

    width_group = parser.add_mutually_exclusive_group()
    width_group.add_argument('--width',
                             help=f'figure width in centimeters (default is {DEFAULT_FIGURE_WIDTH})',
                             type=float,
                             default=DEFAULT_FIGURE_WIDTH)
    width_group.add_argument('--plotWidth',
                             help='width in centimeters of the plotting (central) part',
                             type=float,
                             default=None)

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

    parser.add_argument('--fontSize',
                        help='Font size for the labels of the plot (default is 0.3 * figure width)',
                        type=float)

    parser.add_argument('--dpi',
                        help='Resolution for the image in case the'
                             ' ouput is a raster graphics image (e.g png, jpg) (default is 72)',
                        type=int,
                        default=72
                        )

    parser.add_argument('--trackLabelFraction',
                        help='By default the space dedicated to the track labels is 0.05 of the'
                             ' plot width. This fraction can be changed with this parameter if needed.',
                        default=0.05,
                        type=float)

    parser.add_argument('--trackLabelHAlign',
                        help='By default, the horizontal alignment of the track '
                             'labels is left. This alignemnt can be changed to '
                             'right or center.',
                        default='left',
                        choices=['left', 'right', 'center'])

    parser.add_argument('--decreasingXAxis',
                        help='By default, the x-axis is increasing. '
                             'Use this option if you want to see all tracks'
                             ' with a decreasing x-axis.',
                        action='store_true')

    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    # Identify the regions to plot:
    if args.BED:
        regions = []
        for line in args.BED.readlines():
            try:
                chrom, start, end = line.strip().split('\t')[0:3]
            except ValueError:
                continue
            try:
                start, end = map(int, [start, end])
            except ValueError as detail:
                warnings.warn(f"Invalid value found at line\t{line}\t. {detail}\n")
                continue
            regions.append((chrom, start, end))
    else:
        regions = [get_region(args.region)]

    if len(regions) == 0:
        raise InputError("There is no valid regions to plot.")

    # Create all the tracks
    trp = PlotTracks(args.tracks.name, args.width, fig_height=args.height,
                     fontsize=args.fontSize, dpi=args.dpi,
                     track_label_width=args.trackLabelFraction,
                     plot_regions=regions, plot_width=args.plotWidth)

    # Create dir if dir does not exists:
    # Modified from https://stackoverflow.com/questions/12517451/automatically-creating-directories-with-file-output
    os.makedirs(os.path.dirname(os.path.abspath(args.outFileName)), exist_ok=True)

    # Plot them
    if args.BED:
        name = args.outFileName.split(".")
        file_suffix = name[-1]
        file_prefix = ".".join(name[:-1])
        for chrom, start, end in regions:
            file_name = f"{file_prefix}_{chrom}-{start}-{end}.{file_suffix}"
            if end - start < 200000:
                warnings.warn("A region shorter than 200kb has been "
                              "detected! This can be too small to return "
                              "a proper TAD plot!\n")
            sys.stderr.write(f"saving {file_name}\n")
            current_fig = trp.plot(file_name, chrom, start, end, title=args.title,
                                   h_align_titles=args.trackLabelHAlign,
                                   decreasing_x_axis=args.decreasingXAxis)
            plt.close(current_fig)
    else:
        current_fig = trp.plot(args.outFileName, *regions[0], title=args.title,
                               h_align_titles=args.trackLabelHAlign,
                               decreasing_x_axis=args.decreasingXAxis)
        plt.close(current_fig)
    trp.close_files()
