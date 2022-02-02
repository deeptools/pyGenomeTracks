import argparse
import os
import sys
from pygenometracks.tracksClass import PlotTracks

from pygenometracks._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Facilitates the creation of a configuration file for pyGenomeTracks. The program takes a list '
                                                 'of files and does the boilerplate for the configuration file.',
                                     usage="%(prog)s --trackFiles <bigwig file> <bed file> etc. -o tracks.ini")

    # define the arguments
    parser.add_argument('--trackFiles', '-f',
                        help='Files to use in for the tracks. The ending of the file is used to define the type of '
                             'track. E.g. `.bw` for bigwig, `.bed` for bed etc. For a arcs or links file, the file '
                             'ending recognized is `.arcs` or `.links`',
                        nargs='+',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('--out', '-o',
                        help='File to save the tracks',
                        metavar='output file',
                        type=argparse.FileType('w'),
                        required=True)

    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    available_tracks = PlotTracks.get_available_tracks()

    args.out.write("""
[x-axis]
#optional
#fontsize = 20
# default is bottom meaning below the axis line
# where = top

[spacer]
# height of space in cm (optional)
height = 0.5

""")
    for file_h in args.trackFiles:
        track_added = False
        label = ".".join(os.path.basename(file_h.name).split(".")[0:-1])
        for track_type, track_class in available_tracks.items():
            for ending in track_class.SUPPORTED_ENDINGS:
                if file_h.name.endswith(ending):
                    default_values = track_class.OPTIONS_TXT
                    default_values = default_values.replace("title =", f"title = {label}")
                    args.out.write(f"\n[{label}]\nfile = {file_h.name}\n{default_values}")

                    sys.stdout.write(f"Adding {track_type} file: {file_h.name}\n")
                    track_added = True

        if track_added is False:
            sys.stdout.write(f"WARNING: file format not recognized for: {file_h.name}\n")
