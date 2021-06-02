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
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_all_labels_inside.ini --region X:3215000-3240000 --trackLabelFraction 0.2 --width 38 --dpi 130 --trackLabelHAlign right --decreasingXAxis -o ./pygenometracks/tests/test_data/master_bed_all_label_inside_dec.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_colormap_genes.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_colormap_genes.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_shuffle.ini --BED ./pygenometracks/tests/test_data/regions_chr1XY.bed --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_bed_shuffle.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bed_vlines.ini --BED ./pygenometracks/tests/test_data/regionsXfakeChr.bed --trackLabelFraction 0.5 --width 38 --dpi 130 --trackLabelHAlign center -o ./pygenometracks/tests/test_data/master_bed_vlines.png

# non_classical_bed
bin/pgt --tracks pygenometracks/tests/test_data/bed_unusual_formats.ini --region X:20000-40000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_bed_unusual_formats.png
bin/pgt --tracks pygenometracks/tests/test_data/strange_strand.ini --region chr1:0-500 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_strange_strand.png
bin/pgt --tracks pygenometracks/tests/test_data/invalid_rgb.ini --region chr1:0-500 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_invalid_rgb.png
bin/pgt --tracks pygenometracks/tests/test_data/invalid_blockCount.ini --region chrX:15000-24000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_invalid_blockCount.png
bin/pgt --tracks pygenometracks/tests/test_data/invalid_CDScoo.ini --region chrX:15000-24000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_invalid_CDScoo.png
bin/pgt --tracks pygenometracks/tests/test_data/invalid_blocks.ini --region chrX:15000-24000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_invalid_blocks.png
bin/pgt --tracks pygenometracks/tests/test_data/invalid_score.ini --region chrX:15000-24000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_invalid_score.png
bin/pgt --tracks pygenometracks/tests/test_data/bed_different_UTR.ini --region chr1:0-500 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_different_UTR.png

# test_bedGraphMatrixTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph.ini --region X:2850000-3150000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_bedgraph.png

# test bedGraphTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph_useMid.ini --BED ./pygenometracks/tests/test_data/regions_imbricated_chr2.bed --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_bedgraph_useMid.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph_useMid.ini --region chr2:73,800,000-75,744,000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_bedgraph_useMid.pdf
bin/pgt --tracks ./pygenometracks/tests/test_data/operation_bdg.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_operation_bdg.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph_withNA.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_bedgraph_withNA.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph_negative.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_negative.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph_end_not_covered.ini --region chr7:100-400 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_bedgraph_end_not_covered.png
bin/pgt --tracks ./pygenometracks/tests/test_data/bedgraph_useMid_op.ini --region chr2:74000000-74800000 --trackLabelFraction 0.2 --dpi 130  -o ./pygenometracks/tests/test_data/master_bedgraph_useMid_op.png

# test bigWigTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/bigwig.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_bigwig.png
bin/pgt --tracks ./pygenometracks/tests/test_data/alpha.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_alpha.png
bin/pgt --tracks ./pygenometracks/tests/test_data/hlines.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_hlines.png
bin/pgt --tracks ./pygenometracks/tests/test_data/operation.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_operation.png
bin/pgt --tracks ./pygenometracks/tests/test_data/operation.ini --region fakeChr:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_operation_fakeChr.png
bin/pgt --tracks ./pygenometracks/tests/test_data/grid.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_grid.png
bin/pgt --tracks ./pygenometracks/tests/test_data/example_bigwig.ini --region X:2,500,000-3,000,000 -o ./pygenometracks/tests/test_data/master_example_bigwig.png
bin/pgt --tracks ./pygenometracks/tests/test_data/example_op.ini --region 2L:0-1000 -o ./pygenometracks/tests/test_data/master_operation_2L.png
bin/pgt --tracks ./pygenometracks/tests/test_data/example_bigwig.ini --region X:2,500,000-3,000,000 --width 12 -o ./pygenometracks/tests/test_data/master_example_bigwig_width12.png
bin/pgt --tracks ./pygenometracks/tests/test_data/example_bigwig.ini --region X:2,500,000-3,000,000 --plotWidth 12 -o ./pygenometracks/tests/test_data/master_example_bigwig_plotwidth12.png
bin/pgt --tracks ./pygenometracks/tests/test_data/example_bigwig.ini --region X:2,500,000-3,000,000 --plotWidth 12 --trackLabelFraction 0.5 -o ./pygenometracks/tests/test_data/master_example_bigwig_plotwidth12Lab0.5.png
bin/pgt --tracks ./pygenometracks/tests/test_data/example_op2.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_operation_scien.png

# test_epilogosTrack:
bin/pgt --tracks ./pygenometracks/tests/test_data/epilogos.ini --region X:3100000-3150000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_epilogos.png

# test_hiCMatrixTracks:
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_plot_hic.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_rasterize_height.ini --BED ./pygenometracks/tests/test_data/regions_XY.bed --trackLabelFraction 0.23 --width 38 --dpi 10 -o ./pygenometracks/tests/test_data/master_plot_hic_rasterize_height.pdf
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_log-log.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_plot_hic_log-log.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 --dpi 130 --decreasingXAxis -o ./pygenometracks/tests/test_data/master_plot_hic_dec.png
bin/pgt --tracks ./pygenometracks/tests/test_data/mcool.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_mcool.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_small_test.ini --BED ./pygenometracks/tests/test_data/regions_chr1XY.bed --trackLabelFraction 0.23 --width 38 -o ./pygenometracks/tests/test_data/master_plot_hic_small_test.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_small_test.ini --region 1:0-200000 --trackLabelFraction 0.23 --width 38 -o ./pygenometracks/tests/test_data/master_plot_hic_small_test.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_small_test.ini --region chrM:0-20000 --trackLabelFraction 0.23 --width 38 -o ./pygenometracks/tests/test_data/master_plot_hic_small_test_chrM.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_small_test.ini --region chr1:0-5000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_hic_small_test_small_region.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_small_test.ini --region Y:90000000-100000000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_hic_small_test_above_chrY.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_small_test_2.ini --region chr1:0-100000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_hic_small_test_2.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_one_interaction_cool.ini --BED ./pygenometracks/tests/test_data/regions_chr1XY.bed --trackLabelFraction 0.23 --width 38 -o ./pygenometracks/tests/test_data/master_plot_hic_one_interaction_withBED.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_one_interaction_cool.ini --region chrY:0-1000000 --trackLabelFraction 0.23 --width 38 -o ./pygenometracks/tests/test_data/master_plot_hic_one_interaction_withRegion_chrY-0-1000000.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_one_interaction_diag_h5.ini --region chrY:0-1000000 --trackLabelFraction 0.23 --width 38 -o ./pygenometracks/tests/test_data/master_plot_hic_one_interaction_diag_chrY-0-1000000.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_small_test_3.ini --region chr1:0-200000 --trackLabelFraction 0.23 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_hic_small_test_3.png
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks_hic_force_scale.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 --dpi 130 --fontSize 16 -o ./pygenometracks/tests/test_data/master_plot_hic_force_scale.png

# test_logScale:
bin/pgt --tracks ./pygenometracks/tests/test_data/log1p.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_log1p.png
bin/pgt --tracks ./pygenometracks/tests/test_data/log.ini  --region chr2:73,800,000-75,744,000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_log.png
bin/pgt --tracks ./pygenometracks/tests/test_data/log_grid.ini  --region chr2:73,800,000-75,744,000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_log_grid.png
bin/pgt --tracks ./pygenometracks/tests/test_data/log1p_grid.ini --region X:2700000-3100000 --trackLabelFraction 0.2 --dpi 130 -o ./pygenometracks/tests/test_data/master_log1p_grid.png
bin/pgt --tracks ./pygenometracks/tests/test_data/log_more.ini  --region chr2:73,800,000-75,744,000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_log_more.png

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
bin/pgt --tracks ./pygenometracks/tests/test_data/browser_tracks.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130  --decreasingXAxis -o ./pygenometracks/tests/test_data/master_plot_dec.png
bin/pgt --tracks ./pygenometracks/tests/test_data/firstTrackOverlay.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130  -o ./pygenometracks/tests/test_data/master_empty2.png


# test_tracks_label:
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.2.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.5.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.5_fs20.png --fontSize 20
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --trackLabelHAlign right --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.5_ral.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --trackLabelHAlign center --width 38 --dpi 130 -o ./pygenometracks/tests/test_data/master_title_0.5_cal.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --trackLabelFraction 0.5 --trackLabelHAlign center --width 38 --dpi 250 -o ./pygenometracks/tests/test_data/master_title_0.5_cal_d250.png
bin/pgt --tracks ./pygenometracks/tests/test_data/title.ini --region X:3000000-3500000 --height 10 --title force_height -o ./pygenometracks/tests/test_data/master_title_force_height.png

# tests arcs
bin/pgt --tracks ./pygenometracks/tests/test_data/short_long_arcs.ini --region chr11:40000000-46000000  --trackLabelFraction 0.2 --width 38 --dpi 130  -o pygenometracks/tests/test_data/master_short_long_arcs.png
bin/pgt --tracks ./pygenometracks/tests/test_data/arcs_use_middle.ini --region X:3000000-3300000  --trackLabelFraction 0.2 --width 38 --dpi 130  -o pygenometracks/tests/test_data/master_arcs_use_middle.png
bin/pgt --tracks ./pygenometracks/tests/test_data/arcs_no_score.ini --region X:3000000-3300000  --trackLabelFraction 0.2 --width 38 --dpi 130  -o pygenometracks/tests/test_data/master_arcs_no_score.png

# tests scaleBar
bin/pgt --tracks pygenometracks/tests/test_data/scale_bar.ini --region X:3200000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_scale_bar_zoom.png
bin/pgt --tracks pygenometracks/tests/test_data/scale_bar.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_scale_bar.png
bin/pgt --tracks pygenometracks/tests/test_data/scale_bar_startend.ini --region X:3000000-3600000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_scale_bar_startend.png
bin/pgt --tracks pygenometracks/tests/test_data/scale_bar_startend.ini --region X:3000000-3600000 --trackLabelFraction 0.2 --width 38 --dpi 130 --decreasingXAxis -o pygenometracks/tests/test_data/master_scale_bar_startend_dec.png
bin/pgt --tracks pygenometracks/tests/test_data/scale_bar_startend.ini --region X:2000000000-2500000000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_scale_bar_startend_outside.png
bin/pgt --tracks pygenometracks/tests/test_data/scale_bar_startend.ini --region X:3199500-3201000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_scale_bar_startend_superzoom.png

# tests maf
bin/pgt --tracks pygenometracks/tests/test_data/first_maf.ini --BED pygenometracks/tests/test_data/regions_maf.bed --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_first_maf.png
bin/pgt --tracks pygenometracks/tests/test_data/first_maf.ini --region chr1:0-1000 --trackLabelFraction 0.2 --width 38 --dpi 130 -o pygenometracks/tests/test_data/master_first_maf_empty_chr.png
