from tempfile import NamedTemporaryFile
import os.path
import filecmp
import pygenometracks.makeTracksFile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")
# the relative path is needed to compare the two .ini files correctly.
# Otherwise the paths will differ.
relative_path = os.path.relpath(ROOT)


def test_make_tracks():
    outfile = NamedTemporaryFile(suffix='.ini', prefix='pyGenomeTracks_test_', delete=False)
    args = "--trackFiles {0} {1} {2} {3} --out {4}" \
           "".format(os.path.join(relative_path, 'Li_et_al_2015.h5'),
                     os.path.join(relative_path, 'bigwig_chrx_2e6_5e6.bw'),
                     os.path.join(relative_path, 'tad_classification.bed'),
                     os.path.join(relative_path, 'epilog.qcat.bgz'),
                     outfile.name).split()
    print("using args: {}".format(" ".join(args)))
    pygenometracks.makeTracksFile.main(args)

    if filecmp.cmp(outfile.name,
                   os.path.join(ROOT, 'master_tracks.ini')) is False:
        import difflib

        diff = difflib.unified_diff(open(outfile.name).readlines(),
                                    os.path.join(ROOT, 'master_tracks.ini').readlines(), lineterm='')
        print(''.join(list(diff)))
    assert(filecmp.cmp(outfile.name,
                       os.path.join(ROOT, 'master_tracks.ini')) is True)

    os.remove(outfile.name)
