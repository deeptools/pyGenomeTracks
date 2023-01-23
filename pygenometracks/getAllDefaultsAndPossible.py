"""
This python script will generate two files:
- docs/content/all_default_properties_rst.txt
This file is a rst table with all the defaults values for each parameter
for each track class or vtype. This table is included in the readthedocs
- docs/content/all_possible_properties.txt
This file is a markdown list with possible values.
This can also be used in the readthedocs
"""
from pygenometracks.tracksClass import PlotTracks, DEFAULT_TRACK_HEIGHT
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
starText = ('\n\n'
            r'\* While pyGenomeTracks can convert coverage tracks on the fly,'
            ' this might be a time-consuming step, especially on large files and'
            ' if you want to replot many times. In this situation, we recommend'
            ' using the deepTools suite to convert your files in advance. For'
            ' example `bamCoverage <https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html>`_'
            ' or `bamCompare <https://deeptools.readthedocs.io/en/develop/content/tools/bamCompare.html>`_')

defaut_for_all = {'overlay_previous': 'no'}
possible_for_all = {'overlay_previous': ['no', 'yes', 'share-y']}


def main():
    all_tracks = PlotTracks.get_available_tracks()
    my_prefered_order_tracks_names = ['x_axis', 'epilogos', 'links',
                                      'domains', 'bed', 'gtf', 'narrow_peak',
                                      'bigwig', 'bedgraph', 'bedgraph_matrix',
                                      'hlines', 'hic_matrix', 'hic_matrix_square',
                                      'maf', 'scalebar']
    my_prefered_order_tracks_names = [k for k in my_prefered_order_tracks_names
                                      if k in all_tracks]
    other_tracks = list(set([k for k in all_tracks.keys() if k is not None])
                        - set(my_prefered_order_tracks_names))

    all_types = PlotTracks.get_available_types()
    # Get all possible and default parameters
    all_default_parameters = {}
    all_tracks_with_default = []
    all_possible_parameters = {}
    for track_type in my_prefered_order_tracks_names + other_tracks:
        track_class = all_tracks[track_type]
        for p, value in defaut_for_all.items():
            all_default_parameters[p] = all_default_parameters.get(p, {})
            all_default_parameters[p][track_type] = value
        has_default = False
        for p, value in track_class.DEFAULTS_PROPERTIES.items():
            if p != 'region':
                all_default_parameters[p] = all_default_parameters.get(p, {})
                all_default_parameters[p][track_type] = value
                has_default = True
        if has_default:
            all_tracks_with_default += [track_type]

        for p, value in possible_for_all.items():
            all_possible_parameters[p] = all_possible_parameters.get(p, {})
            all_possible_parameters[p][track_type] = value
        for p, value in track_class.POSSIBLE_PROPERTIES.items():
            all_possible_parameters[p] = all_possible_parameters.get(p, {})
            all_possible_parameters[p][track_type] = value

        for p in track_class.BOOLEAN_PROPERTIES:
            all_possible_parameters[p] = all_possible_parameters.get(p, {})
            all_possible_parameters[p][track_type] = ["true", "false"]

    for my_type in all_types:
        type_class = all_types[my_type]
        has_default = False
        for p, value in type_class.DEFAULTS_PROPERTIES.items():
            if p != 'region':
                all_default_parameters[p] = all_default_parameters.get(p, {})
                all_default_parameters[p][my_type] = value
                has_default = True
        if has_default:
            all_tracks_with_default += [my_type]

        for p, value in type_class.POSSIBLE_PROPERTIES.items():
            all_possible_parameters[p] = all_possible_parameters.get(p, {})
            all_possible_parameters[p][my_type] = value

        for p in type_class.BOOLEAN_PROPERTIES:
            all_possible_parameters[p] = all_possible_parameters.get(p, {})
            all_possible_parameters[p][my_type] = ["true", "false"]

    # For the default they are summarized in a matrix
    mat = np.empty((len(all_default_parameters) + 2, len(all_tracks_with_default) + 1),
                   dtype='U100')
    mat[0, 0] = 'parameter'
    mat[1, 0] = '--'
    for j, track_type in enumerate(all_tracks_with_default, start=1):
        mat[0, j] = f":doc:`tracks/{track_type}`"
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
               mat, fmt=f'%-{max_char}s', delimiter="  ",
               header='  '.join(mat[1, :]), footer='  '.join(mat[1, :]),
               comments='')
    # For export as csv remove the row with '=='
    mat_csv = np.delete(mat, 1, 0)
    # update the header:
    for j, track_type in enumerate(all_tracks_with_default, start=1):
        mat_csv[0, j] = track_type
    np.savetxt(os.path.join("docs", "content", "all_default_properties.csv"),
               mat_csv, fmt='%s', delimiter=",")

    if starPut:
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
            fo.write(f"- **{p}**:\n\n")
            for name in [k for k in possible_values]:
                reformated_possible = ", ".join([v for v in possible_values[name]
                                                if v is not None])
                if None in possible_values[name]:
                    reformated_possible += ", " + not_set_string
                fo.write(f"  - for *{name}*: {reformated_possible}\n\n")

    # For the track description:
    for track_type in my_prefered_order_tracks_names + other_tracks:
        starPut = False
        track_class = all_tracks[track_type]
        with open(os.path.join("docs", "content", "tracks", "auto", f"{track_type}_deduced_from_code.txt"),
                  'w') as fo:
            fo.write("Necessary:\n")
            fo.write("^^^^^^^^^^\n")
            for n in track_class.NECESSARY_PROPERTIES:
                fo.write("- **" + n + "**\n\n")
            fo.write("Optional:\n")
            fo.write("^^^^^^^^^\n")
            fo.write("- **title**: Put here a title which will apprear on the right.\n\n")
            if track_type == "hic_matrix":
                fo.write("- **height**: If you do not set it, the height will "
                         "be adjusted so that each bin is a square, else you "
                         "can choose any float above 0.\n\n")
            elif track_type == "x_axis":
                fo.write("- **height**: If you do not set it, the height will "
                         "be ``fontsize / 8``, else you "
                         "can choose any float above 0.\n\n")
            else:
                default = DEFAULT_TRACK_HEIGHT
                fo.write(f"- **height**: `{default}` (default) or float above 0.\n\n")
            for p in all_default_parameters:
                if track_type in all_default_parameters[p]:
                    if p in putStarsAfter:
                        line_start = (f"- **{p}"
                                      r"\***: ")
                        starPut = True
                    else:
                        line_start = f"- **{p}**: "
                    default = all_default_parameters[p][track_type]
                    if isinstance(default, bool):
                        if default:
                            fo.write(f"{line_start}`true` (default) or false.\n\n")
                        else:
                            fo.write(f"{line_start}`false` (default) or true.\n\n")
                    else:
                        if default is None:
                            default_str_no_ind = "by default this option is not set"
                            default_str = f"{default_str_no_ind} but you can also put:"
                        else:
                            default_str_no_ind = f"`{default}` (default)"
                            default_str = f"{default_str_no_ind} or"

                        if track_type in all_possible_parameters.get(p, {}):
                            others = [str(v) for v in all_possible_parameters[p][track_type] if v != default]
                            if len(others) > 1:
                                others_str = ", ".join(others[:-1]) + " or " + others[-1]
                            else:
                                others_str = others[0]
                            fo.write(f"{line_start}{default_str} {others_str}.\n\n")
                        elif p in track_class.FLOAT_PROPERTIES or p in track_class.INTEGER_PROPERTIES:
                            if p in track_class.FLOAT_PROPERTIES:
                                range_constrains = track_class.FLOAT_PROPERTIES[p]
                                type_p = "float"
                            else:
                                range_constrains = track_class.INTEGER_PROPERTIES[p]
                                type_p = "integer"
                            range_str = ""
                            if range_constrains[0] != - np.inf:
                                range_str = f"{range_str} above {range_constrains[0]}"
                            if range_constrains[1] != np.inf:
                                range_str = f"{range_str} below {range_constrains[1]}"
                            fo.write(f"{line_start}{default_str} any {type_p}{range_str}\n\n")
                        else:
                            fo.write(f"{line_start}{default_str_no_ind}\n\n")
            if starPut:
                fo.write(starText)
        with open(os.path.join("docs", "content", "tracks", "auto", f"{track_type}_options_text.txt"),
                  'w') as fo:
            fo.write(track_class.OPTIONS_TXT)

    # For the type description:
    for my_type in all_types:
        type_class = all_types[my_type]
        starPut = False
        with open(os.path.join("docs", "content", "tracks", "auto", f"{my_type}_deduced_from_code.txt"),
                  'w') as fo:
            fo.write("Necessary:\n")
            fo.write("^^^^^^^^^^\n")
            fo.write(f"- **type**: {my_type}\n")
            for n in type_class.NECESSARY_PROPERTIES:
                fo.write("- **" + n + "**\n\n")
            fo.write("Optional:\n")
            fo.write("^^^^^^^^^\n")
            for p in all_default_parameters:
                if my_type in all_default_parameters[p]:
                    if p in putStarsAfter:
                        line_start = (f"- **{p}"
                                      r"\***: ")
                        starPut = True
                    else:
                        line_start = f"- **{p}**: "
                    default = all_default_parameters[p][my_type]
                    if isinstance(default, bool):
                        if default:
                            fo.write(f"{line_start}`true` (default) or false.\n\n")
                        else:
                            fo.write(f"{line_start}`false` (default) or true.\n\n")
                    else:
                        if default is None:
                            default_str_no_ind = "by default this option is not set"
                            default_str = f"{default_str_no_ind} but you can also put:"
                        else:
                            default_str_no_ind = f"`{default}` (default)"
                            default_str = f"{default_str_no_ind} or"

                        if my_type in all_possible_parameters.get(p, {}):
                            others = [str(v) for v in all_possible_parameters[p][my_type] if v != default]
                            if len(others) > 1:
                                others_str = ", ".join(others[:-1]) + " or " + others[-1]
                            else:
                                others_str = others[0]
                            fo.write(f"{line_start}{default_str} {others_str}.\n\n")
                        elif p in type_class.FLOAT_PROPERTIES or p in type_class.INTEGER_PROPERTIES:
                            if p in type_class.FLOAT_PROPERTIES:
                                range_constrains = type_class.FLOAT_PROPERTIES[p]
                                type_p = "float"
                            else:
                                range_constrains = type_class.INTEGER_PROPERTIES[p]
                                type_p = "integer"
                            range_str = ""
                            if range_constrains[0] != - np.inf:
                                range_str = f"{range_str} above {range_constrains[0]}"
                            if range_constrains[1] != np.inf:
                                range_str = f"{range_str} below {range_constrains[1]}"
                            fo.write(f"{line_start}{default_str} any {type_p}{range_str}\n\n")
                        else:
                            fo.write(f"{line_start}{default_str_no_ind}\n\n")
            if starPut:
                fo.write(starText)
        with open(os.path.join("docs", "content", "tracks", "auto", f"{my_type}_options_text.txt"),
                  'w') as fo:
            fo.write(type_class.OPTIONS_TXT)


if __name__ == "__main__":
    main()
