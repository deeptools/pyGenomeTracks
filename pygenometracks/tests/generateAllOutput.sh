# This bash script can be used to regenerate all test output when a global change happens.
# It needs to be launched from the root folder.
# test_bed_and_gtf_tracks:
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_and_gtf_tracks.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_and_gtf.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_arrow_tracks.ini --region X:3120000-3150000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_arrow.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_arrow_tracks.ini --region X:3130000-3140000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_arrow_zoom.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_maxLab_tracks.ini --region X:2000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_maxLab.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_flybase_tracks.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_flybase.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_tssarrow_tracks.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_tssarrow.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_tssarrow_tracks.ini --region X:3020000-3070000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_tssarrow_zoom.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_tssarrow_tracks.ini --region X:3130000-3150000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_tssarrow_zoom2.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_all_labels_inside.ini --region X:3100000-3200000 --trackLabelFraction 0.2 --width 38 --dpi 130 --trackLabelHAlign right -o ./pygenometracks/tests/test_data/master_bed_all_label_inside.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_colormap_genes.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_colormap_genes.png

# test_bedGraphMatrixTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph.ini --region X:2850000-3150000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_bedgraph.png

# test bedGraphTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph_useMid.ini --region chr2:73,800,000-75,744,000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_bedgraph_useMid.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph_useMid.ini --region chr2:73,800,000-75,744,000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_bedgraph_useMid.pdf

# test bigWigTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/bigwig.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_bigwig.png
bin/pgt --tracks ./pygenometracks/tests/test_data/alpha.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_alpha.png
bin/pgt --tracks ./pygenometracks/tests/test_data/hlines.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_hlines.png

# test_epilogosTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/epilogos.ini --region X:3100000-3150000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_epilogos.png

# test_hiCMatrixTracks:
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_plot_hic.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_rasterize_height.ini --region X:2500000-2600000 --trackLabelFraction 0.23 --width 38 --dpi 10 -o ./pygenometracks/tests/test_data/master_plot_hic_rasterize_height.pdf

# test_make_tracks:
bin/make_tracks_file --trackFiles pygenometracks/tests/test_data/Li_et_al_2015.h5 pygenometracks/tests/test_data/bigwig_chrx_2e6_5e6.bw pygenometracks/tests/test_data/tad_classification.bed pygenometracks/tests/test_data/epilog.qcat.bgz -o pygenometracks/tests/test_data/master_tracks.ini

# test_narrowPeakTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/narrow_peak.ini --region X:2760000-2802000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_narrowPeak.png
bin/pgt --tracks ./pygenometracks/tests/test_data/narrow_peak2.ini --region X:2760000-2802000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_narrowPeak2.png

# test_plot_tracks:
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_plot.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks.ini --region Y:0-1000000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_plot_3.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks.ini --region X:0-1000000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_plot_2.png
bin/pgt --tracks ./pygenometracks/tests/test_data/empty.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_empty.png

# test_tracks_label:
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.2.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.5.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.5_fs20.png --fontSize 20
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --trackLabelHAlign right --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.5_ral.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --trackLabelHAlign center --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.5_cal.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --trackLabelHAlign center --width 38 --dpi 250 -o ./pygenometracks/tests/test_data/master_title_0.5_cal_d250.png

# tests arcs
bin/pgt --tracks ./pygenometracks/tests/test_data/short_long_arcs.ini --region chr11:40000000-46000000  --trackLabelFraction 0.2 --width 38 --dpi 130  -o pygenometracks/tests/test_data/master_short_long_arcs.png
