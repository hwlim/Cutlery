Plan:
1. Footprint analysis using MEME suite
- Motivation: When sorted by motif score, the contrast of CnR footprint become less clear, which may help us to decide optimal threshold of fimo motif score
- meme motif search for seed motif (primarily because homer often gives limited length motifs:
using random 1000 sequence or
multiple meme runs for multiple disjoint subsets
- Footprint contrast sorted by meme motif score
- motif threshold + contrast threshold

2. Footprint visualization for two pairs of bigWig files for single anchor bed file
- Primary motif is to compare before/after bias correction side by side anchoring on the set of regions
- 2 x 3 panels: add 1) average plot and 2) heatmap to the existing visualization routine for an additional bigwig file

