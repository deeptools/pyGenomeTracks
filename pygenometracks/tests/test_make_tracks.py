from tempfile import NamedTemporaryFile
import os.path
import filecmp
import pygenometracks.makeTracksFile

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
# the relative path is needed to compare the two .ini files correctly. Otherwise
# the paths will differ.
relative_path = os.path.relpath(ROOT)


def test_make_tracks():
    outfile = NamedTemporaryFile(suffix='.ini', prefix='pyGenomeTracks_test_', delete=False)
    args = "--trackFiles {0}/Li_et_al_2015.h5 {0}/bigwig_chrx_2e6_5e6.bw {0}/tad_classification.bed " \
           "{0}/epilog.qcat.bgz " \
           "--out {1}".format(relative_path, outfile.name).split()
    print("using args: {}".format(" ".join(args)))
    pygenometracks.makeTracksFile.main(args)

    if filecmp.cmp(outfile.name, ROOT + 'master_tracks.ini') is False:
        import difflib

        diff = difflib.unified_diff(open(outfile.name).readlines(),
                                    open(ROOT + 'master_tracks.ini').readlines(), lineterm='')
        print(''.join(list(diff)))
    assert(filecmp.cmp(outfile.name, ROOT + 'master_tracks.ini') is True)

    os.remove(outfile.name)
