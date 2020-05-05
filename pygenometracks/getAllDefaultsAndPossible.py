"""
This python script will generate two files:
- docs/content/all_default_properties_rst.txt
This file is a rst table with all the defaults values for each parameter
for each track class. This table is included in the readthedocs
- docs/content/all_possible_properties.txt
This file is a markdown list with possible values.
This can also be used in the readthedocs
"""
from pygenometracks.tracksClass import PlotTracks, XAxisTrack
import numpy as np
import os.path

not_used_string = ''
not_set_string = 'not set'
track_separator = ', '

# Here is the keyword that was used in version 3.1.2
# for boolean values
# GOOD_PRACTICES = {'labels': {True: 'on', False: 'off'},
#                   'show_data_range': {True: 'yes', False: 'no'},
#                   'plot_horizontal_lines': {True: 'yes', False: 'no'},
#                   'use_middle': {True: 'yes', False: 'no'},
#                   'rasterize': {True: 'yes', False: 'no'},
#                   'global_max_row': {True: 'yes', False: 'no'},
#                   'show_masked_bins': {True: 'yes', False: 'no'},
#                   'show_labels': {True: 'yes', False: 'no'},
#                   'use_summit': {True: 'yes', False: 'no'},
#                   # 'skip': {True: 'yes', False: 'no'},
#                   'merge_transcripts': {True: 'on', False: 'off'},
#                   'nans_to_zeros': {True: 'True', False: 'False'}}

putStarsAfter = ['second_file', 'operation']
starPut = False
starText = ('\n\n* While pyGenomeTracks can convert coverage tracks on the fly,'
            ' this might be a time-consuming step, especially on large files and'
            ' if you want to replot many times. In this situation, we recommend'
            ' using the deepTools suite to convert your files in advance. For'
            ' example [bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)')

def main():
    all_tracks = PlotTracks.get_available_tracks()
    my_prefered_order_tracks_names = [None, 'epilogos', 'links',
                                      'domains', 'bed', 'gtf', 'narrow_peak',
                                      'bigwig', 'bedgraph', 'bedgraph_matrix',
                                      'hlines', 'hic_matrix']
    my_prefered_order_tracks_names = [k for k in my_prefered_order_tracks_names
                                      if k in all_tracks]
    other_tracks = list(set(all_tracks.keys())
                        - set(my_prefered_order_tracks_names))
    # Get all possible and default parameters
    all_default_parameters = {}
    all_tracks_with_default = []
    all_possible_parameters = {}
    for track_type in my_prefered_order_tracks_names + other_tracks:
        track_class = all_tracks[track_type]
        if track_class == XAxisTrack:
            track_type = "x-axis"
        has_default = False
        for p, value in track_class.DEFAULTS_PROPERTIES.items():
            all_default_parameters[p] = all_default_parameters.get(p, {})
            all_default_parameters[p][track_type] = value
            has_default = True
        if has_default:
            all_tracks_with_default += [track_type]

        for p, value in track_class.POSSIBLE_PROPERTIES.items():
            all_possible_parameters[p] = all_possible_parameters.get(p, {})
            all_possible_parameters[p][track_type] = value

        for p in track_class.BOOLEAN_PROPERTIES:
            all_possible_parameters[p] = all_possible_parameters.get(p, {})
            all_possible_parameters[p][track_type] = ["true", "false"]

    # For the default they are summarized in a matrix
    mat = np.empty((len(all_default_parameters) + 2, len(all_tracks) + 1),
                   dtype='U25')
    mat[0, 0] = 'parameter'
    mat[1, 0] = '--'
    for j, track_type in enumerate(all_tracks_with_default, start=1):
        mat[0, j] = track_type
        mat[1, j] = '-'
        for i, p in enumerate(all_default_parameters):
            if j == 1:
                if p in putStarsAfter:
                    mat[i + 2, 0] = p + "*"
                    starPut = True
                else:
                    mat[i + 2, 0] = p

            default = all_default_parameters[p].get(track_type,
                                                    not_used_string)

            if isinstance(default, bool):
                default = str(default).lower()

            if default is None:
                default = not_set_string

            mat[i + 2, j] = default
    # The matrix is written in a file to be able to use it in the readthedocs
    max_char = max([len(mat[i, j]) for i in range(mat.shape[0]) for j in range(mat.shape[1])])
    mat[1, :] = ['=' * max_char] * mat.shape[1]
    np.savetxt(os.path.join("docs", "content", "all_default_properties_rst.txt"),
               mat, fmt='%-{}s'.format(max_char), delimiter="  ",
               header='  '.join(mat[1, :]), footer='  '.join(mat[1, :]),
               comments='')
    with open(os.path.join("docs", "content", "all_default_properties_rst.txt"), 'a') as f:
        f.write(starText)
    # For the possible:
    with open(os.path.join("docs", "content", "all_possible_properties.txt"),
              'w') as fo:
        for p, possible_dic in all_possible_parameters.items():
            possible_values = {}
            for track_type, pv in possible_dic.items():
                if len(possible_values) == 0:
                    possible_values[track_type] = pv
                else:
                    added = False
                    for k in possible_values:
                        if possible_values[k] == pv:
                            possible_values[k + track_separator + track_type] = pv
                            del possible_values[k]
                            added = True
                            break
                    if not added:
                        possible_values[track_type] = pv
            fo.write("- **" + p + "**:\n\n")
            for name in [k for k in possible_values]:
                reformated_possible = ", ".join([v for v in possible_values[name]
                                                if v is not None])
                if None in possible_values[name]:
                    reformated_possible += ", " + not_set_string
                fo.write("  - for *" + name + "*: " + reformated_possible + "\n\n")


if __name__ == "__main__":
    main()
