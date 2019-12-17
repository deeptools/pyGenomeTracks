from tempfile import NamedTemporaryFile
import os.path
import filecmp
import pygenometracks.fromGtfToBed12

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")


def test_convert_gtf():
    outfile = NamedTemporaryFile(suffix='.bed',
                                 prefix='pyGenomeTracks_test_gtf_',
                                 delete=False)
    args = " --mergeTranscripts --useGene --out {} {}" \
           "".format(outfile.name,
                     os.path.join(ROOT, 'dm3_subset_BDGP5.78.gtf.gz')).split()
    print("using args: {}".format(" ".join(args)))
    pygenometracks.fromGtfToBed12.main(args)

    same_files = filecmp.cmp(outfile.name,
                             os.path.join(ROOT,
                                          'dm3_subset_BDGP5.78.asbed.bed'))
    if not same_files:
        import difflib

        diff = difflib.unified_diff(open(outfile.name).readlines(),
                                    open(os.path.join(ROOT,
                                                      'dm3_subset_BDGP5.78.asbed.bed')).readlines(),
                                    lineterm='')
        print(''.join(list(diff)))
    assert(same_files)

    os.remove(outfile.name)
