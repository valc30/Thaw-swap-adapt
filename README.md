# Thaw-swap-adapt
We combine metagenomics and metatranscriptomics across a permafrost-to-water continuum to resolve the metaplasmidome. We reconstructed a vast, structured, environmentally responsive plasmid reservoir and tracked its host-community. This release gathers permutation statistical test python and R scripts associated with the publication "Thaw-swap-adapt, XX", Valentine Cyriaque, Martial Leroy, Anne-Charlotte Vranckx, Ruddy Wattiez, Jérôme Comte.

Python 3.10
- module_count_prob.py
Calculate probability to find more or less modules than expected by chance after plasmid node shuffeling  using infomap
Usage : \
python infomap_perm_nodes.py \
  --input 	network.txt \
  --n_perm 	1000 \
  --outdir 	results/ \
  --infomap 	/path/to/infomap \
  --seed        3 \
Input: An Infomap multilayer network .txt file in explicit format. The file contains a *Vertices section listing node IDs and names, followed by a *Multilayer section listing edges as layer_u node_u layer_v node_v weight. Each row encodes one interaction between two node-layer pairs.

R version 4.5.1
- modules_sample_taxo_statistics.R \
Calculate probability of finding the same taxonomic class or sample type within modules ; Identify how many modules are pure with each sample_type and shared modules \
Packages: ggplot2 4.0.2
- XX
  Usage \
  Packages \


