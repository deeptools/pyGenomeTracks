from pygenometracks.tracksClass import PlotTracks
import numpy as np

not_used_string = ''
not_set_string = 'not set'
track_separator = ','

# Here is the keyword we want people to use
# for boolean values
# The first one is True, the second is False
GOOD_PRACTICES = {'labels': {True: 'on', False: 'off'},
                  'show data range': {True: 'yes', False: 'no'},
                  'plot horizontal lines': {True: 'yes', False: 'no'},
                  'use middle': {True: 'yes', False: 'no'},
                  'rasterize': {True: 'yes', False: 'no'},
                  'global max row': {True: 'yes', False: 'no'},
                  'show_masked_bins': {True: 'yes', False: 'no'},
                  'show labels': {True: 'yes', False: 'no'},
                  'use summit': {True: 'yes', False: 'no'},
                  # 'skip': {True: 'yes', False: 'no'},
                  'merge transcripts': {True: 'on', False: 'off'}}


def main():
    all_tracks = PlotTracks.get_available_tracks()

    # Get all possible and default parameters
    all_default_parameters = {}
    all_tracks_with_default = []
    all_possible_parameters = {}
    for track_type, track_class in all_tracks.items():
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

    # For the default they are summarized in a matrix
    mat = np.empty((len(all_default_parameters) + 2, len(all_tracks) + 1),
                   dtype='U25')
    mat[0, 0] = 'parameter'
    mat[1, 0] = '--'
    j = 1
    for track_type in all_tracks_with_default:
        mat[0, j] = track_type
        mat[1, j] = '-'
        for i, p in enumerate(all_default_parameters):
            if j == 1:
                mat[i + 2, 0] = p

            default = all_default_parameters[p].get(track_type,
                                                    not_used_string)

            if p in GOOD_PRACTICES and default is not not_used_string:
                default = GOOD_PRACTICES[p][default]

            if default is None:
                default = not_set_string

            mat[i + 2, j] = default
        j += 1
    # The matrix is written in a file to be able to use it in the README.md
    np.savetxt("all_default_properties.txt", mat, fmt='%s', delimiter=" | ")

    # For the possible:
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
        print("- **" + p + "**:")
        for name in [k for k in possible_values]:
            reformated_possible = ", ".join([v for v in possible_values[name] if v is not None])
            if None in possible_values[name]:
                reformated_possible += ", " + not_set_string
            print("  - for *" + name + "*: " + reformated_possible)
    for p, pv in GOOD_PRACTICES.items():
        names = []
        for track_type in all_tracks_with_default:
            if track_type in all_default_parameters[p]:
                names += [track_type]
        print("- **" + p + "**:")
        name = ", ".join(names)
        reformated_possible = ", ".join([v for k, v in pv.items()])
        print("  - for *" + name + "*: " + reformated_possible)


if __name__ == "__main__":
    main()
