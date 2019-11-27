# This bash script can be used to regenerate all test output when a global change happens.
# It needs to be launched from the root folder.
make_tracks_file --trackFiles pygenometracks/tests/test_data/Li_et_al_2015.h5 pygenometracks/tests/test_data/bigwig_chrx_2e6_5e6.bw pygenometracks/tests/test_data/tad_classification.bed pygenometracks/tests/test_data/epilog.qcat.bgz -o pygenometracks/tests/test_data/master_tracks.ini
pgt --tracks ./pygenometracks/tests/test_data/bed_and_gtf_tracks.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_and_gtf.png
pgt --tracks ./pygenometracks/tests/test_data/bedgraph.ini --region X:2850000-3150000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_bedgraph.png
pgt --tracks ./pygenometracks/tests/test_data/bigwig.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_bigwig.png
pgt --tracks ./pygenometracks/tests/test_data/alpha.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_alpha.png
pgt --tracks ./pygenometracks/tests/test_data/epilogos.ini --region X:3100000-3150000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_epilogos.png
pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_plot_hic.png
pgt --tracks ./pygenometracks/tests/test_data/narrow_peak.ini --region X:2760000-2802000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_narrowPeak.png
pgt --tracks ./pygenometracks/tests/test_data/narrow_peak2.ini --region X:2760000-2802000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_narrowPeak2.png
pgt --tracks ./pygenometracks/tests/test_data/browser_tracks.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_plot.png
pgt --tracks ./pygenometracks/tests/test_data/bedgraph_useMid.ini --region chr2:73,800,000-75,744,000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_bedgraph_useMid.png
pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_rasterize_height.ini --region X:2500000-2600000 --trackLabelFraction 0.23 --width 38 --dpi 10 -o ./pygenometracks/tests/test_data/master_plot_hic_rasterize_height.pdf
pgt --tracks ./pygenometracks/tests/test_data/bedgraph_useMid.ini --region chr2:73,800,000-75,744,000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_bedgraph_useMid.pdf
pgt --tracks ./pygenometracks/tests/test_data/bed_maxLab_tracks.ini --region X:2000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_maxLab.png
pgt --tracks ./pygenometracks/tests/test_data/browser_tracks.ini --region Y:0-1000000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_plot_3.png
pgt --tracks ./pygenometracks/tests/test_data/browser_tracks.ini --region X:0-1000000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_plot_2.png
