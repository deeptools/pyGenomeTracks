Possible parameters
===================
Here is a table to summarize which are the parameters that can be use for each of the `file_type` and which is the default value:
Empty means this parameter is not used.
not set means that by default the parameter is commented.

.. <!--- Start of default table -->
parameter | x-axis | epilogos | links | domains | bed | narrow_peak | bigwig | bedgraph | bedgraph_matrix | hlines | hic_matrix
-- | - | - | - | - | - | - | - | - | - | - | -
where | bottom |  |  |  |  |  |  |  |  |  | 
fontsize | 15 |  |  | 12 | 12 |  |  |  |  |  | 
categories_file |  | not set |  |  |  |  |  |  |  |  | 
orientation |  | not set | not set | not set | not set | not set | not set | not set | not set | not set | not set
links_type |  |  | arcs |  |  |  |  |  |  |  | 
line_width |  |  | not set | 0.5 | 0.5 |  |  |  |  | 0.5 | 
line_style |  |  | solid |  |  |  |  |  |  | solid | 
color |  |  | blue | #1f78b4 | #1f78b4 | #FF000080 | #33a02c | #a6cee3 |  | black | 
alpha |  |  | 0.8 |  |  |  | 1 | 1 |  | 1 | 
max_value |  |  | not set | not set | not set | not set | not set | not set | not set | not set | not set
min_value |  |  | not set | not set | not set |  | not set | not set | not set | not set | not set
border_color |  |  |  | black | black |  |  |  |  |  | 
interval_height |  |  |  | 100 | 100 |  |  |  |  |  | 
prefered_name |  |  |  | transcript_name | transcript_name |  |  |  |  |  | 
merge_transcripts |  |  |  | false | false |  |  |  |  |  | 
labels |  |  |  |  | true |  |  |  |  |  | 
style |  |  |  |  | flybase |  |  |  |  |  | 
display |  |  |  |  | stacked |  |  |  |  |  | 
max_labels |  |  |  |  | 60 |  |  |  |  |  | 
global_max_row |  |  |  |  | false |  |  |  |  |  | 
gene_rows |  |  |  |  | not set |  |  |  |  |  | 
arrow_interval |  |  |  |  | 2 |  |  |  |  |  | 
arrowhead_included |  |  |  |  | false |  |  |  |  |  | 
color_utr |  |  |  |  | grey |  |  |  |  |  | 
height_utr |  |  |  |  | 1 |  |  |  |  |  | 
show_data_range |  |  |  |  |  | true | true | true | true | true | 
show_labels |  |  |  |  |  | true |  |  |  |  | 
use_summit |  |  |  |  |  | true |  |  |  |  | 
width_adjust |  |  |  |  |  | 1.5 |  |  |  |  | 
type |  |  |  |  |  | peak | fill | fill | matrix |  | 
negative_color |  |  |  |  |  |  | not set | not set |  |  | 
nans_to_zeros |  |  |  |  |  |  | false | false |  |  | 
summary_method |  |  |  |  |  |  | mean | not set |  |  | 
number_of_bins |  |  |  |  |  |  | 700 | 700 |  |  | 
use_middle |  |  |  |  |  |  |  | false |  |  | 
rasterize |  |  |  |  |  |  |  | false | true |  | true
pos_score_in_bin |  |  |  |  |  |  |  |  | center |  | 
plot_horizontal_lines |  |  |  |  |  |  |  |  | false |  | 
region |  |  |  |  |  |  |  |  |  |  | not set
depth |  |  |  |  |  |  |  |  |  |  | 100000
show_masked_bins |  |  |  |  |  |  |  |  |  |  | false
scale_factor |  |  |  |  |  |  |  |  |  |  | 1
transform |  |  |  |  |  |  |  |  |  |  | no
colormap |  |  |  |  |  |  |  |  |  |  | RdYlBu_r
.. <!--- End of default table -->

Some parameters can take only discrete values.

They are summarized here:
.. <!--- Start of possible table -->
- **where**:
  - for *x-axis*: top, bottom
- **orientation**:
  - for *epilogos, links, domains, bed, narrow_peak, bigwig, bedgraph, bedgraph_matrix, hlines, hic_matrix*: inverted, not set
- **links_type**:
  - for *links*: arcs, triangles, loops
- **line_style**:
  - for *links, hlines*: solid, dashed, dotted, dashdot
- **style**:
  - for *bed*: flybase, UCSC
- **display**:
  - for *bed*: collapsed, triangles, interleaved, stacked
- **type**:
  - for *narrow_peak*: peak, box
  - for *bedgraph_matrix*: matrix, lines
- **summary_method**:
  - for *bigwig*: mean, average, max, min, stdev, dev, coverage, cov, sum
  - for *bedgraph*: mean, average, max, min, stdev, dev, coverage, cov, sum, not set
- **pos_score_in_bin**:
  - for *bedgraph_matrix*: center, block
- **transform**:
  - for *hic_matrix*: no, log, log1p, -log
- **labels**:
  - for *bed*: true, false
- **show_data_range**:
  - for *narrow_peak, bigwig, bedgraph, bedgraph_matrix, hlines*: true, false
- **plot_horizontal_lines**:
  - for *bedgraph_matrix*: true, false
- **use_middle**:
  - for *bedgraph*: true, false
- **rasterize**:
  - for *bedgraph, bedgraph_matrix, hic_matrix*: true, false
- **global_max_row**:
  - for *bed*: true, false
- **show_masked_bins**:
  - for *hic_matrix*: true, false
- **show_labels**:
  - for *narrow_peak*: true, false
- **use_summit**:
  - for *narrow_peak*: true, false
- **merge_transcripts**:
  - for *domains, bed*: true, false
- **nans_to_zeros**:
  - for *bigwig, bedgraph*: true, false
- **arrowhead_included**:
  - for *bed*: true, false
.. <!--- End of possible table -->
