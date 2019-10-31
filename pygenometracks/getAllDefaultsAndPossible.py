from pygenometracks.tracks import *
from pygenometracks.tracksClass import PlotTracks
import numpy as np

not_used_string = ''
not_set_string = 'not set'
track_separator = ','

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
    mat[1, 0] = '-'
    j = 1
    for track_type in all_tracks_with_default:
        mat[0, j] = track_type
        mat[1, j] = '-'
        for i, p in enumerate(all_default_parameters):
            if j == 1:
                mat[i + 2, 0] = p
            default = all_default_parameters[p].get(track_type,
                                                    not_used_string)
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
        print(p + ":")
        for name in [k for k in possible_values]:
            reformated_possible = ", ".join([v for v in possible_values[name] if v is not None])
            if None in possible_values[name]:
                reformated_possible += ", " + not_set_string
            print("\tfor " + name + ": " + reformated_possible)

if __name__ == "__main__":
    main()
