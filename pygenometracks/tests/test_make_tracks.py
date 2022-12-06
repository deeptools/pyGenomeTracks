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
    outfile = NamedTemporaryFile(suffix='.ini', prefix='pyGenomeTracks_test_',
                                 delete=False)
    all_paths = [os.path.join(relative_path, 'Li_et_al_2015.h5'),
                 os.path.join(relative_path, 'bigwig_chrx_2e6_5e6.bw'),
                 os.path.join(relative_path, 'tad_classification.bed'),
                 os.path.join(relative_path, 'epilog.qcat.bgz')]
    args = f"--trackFiles {all_paths[0]} {all_paths[1]} {all_paths[2]}"\
           f" {all_paths[3]} --out {outfile.name}".split()
    pygenometracks.makeTracksFile.main(args)

    if filecmp.cmp(outfile.name,
                   os.path.join(ROOT, 'master_tracks.ini')) is False:
        import difflib

        diff = difflib.unified_diff(open(outfile.name).readlines(),
                                    open(os.path.join(ROOT, 'master_tracks.ini')).readlines(), lineterm='')
        print(''.join(list(diff)))

    assert filecmp.cmp(outfile.name,
                       os.path.join(ROOT, 'master_tracks.ini')) is True

    os.remove(outfile.name)
