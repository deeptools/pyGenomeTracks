from tempfile import NamedTemporaryFile
import os.path
import filecmp
import pygenometracks.makeTracksFile

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def test_hicPlotTads():
    outfile = NamedTemporaryFile(suffix='.ini', prefix='pyGenomeTracks_test_', delete=False)
    args = "--trackFiles {0}/bigwig_chrx_2e6_5e6.bw {0}/tad_classification.bed  " \
           "--out {1}".format(ROOT, outfile.name).split()
    pygenometracks.makeTracksFile.main(args)

    assert(filecmp.cmp(outfile.name, ROOT + '/master_tracks.ini') is True)

    os.remove(outfile.name)
