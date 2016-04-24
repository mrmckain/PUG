[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19132.svg)](http://dx.doi.org/10.5281/zenodo.19132)

PUG--Phylogenetic Placement of Polyploidy Using Genomes
=============
<b>Author</b>: Michael R. McKain<br>
</br>
Version 2.0, June 30, 2015<br>
</br>
Version 2.1, February 1, 2016
<br></br>
<b>Contact</b>: https://github.com/mrmckain
<h3>Description</h3>

PUG is designed to test the phylogenetic placment of a polyploid event over mulitple gene trees using a known species tree and a set of putative paralogs.
Despite the name, PUG works with genomic or transcriptomic data. 

Paralogs are queried in gene trees. Their placement, and the relative placement of taxa in the gene tree, are compared to the given species tree and if they are found to match a node of the species tree, the pair is counted towards supporting the polyploid event origin at that node. A single node may be represented in a gene tree by multiple paralog pairs, but that node is counted only once for the gene tree, if the "unique" option is used.  The "all pairs" option will all multiple mappings to the same node.

<h4>Requirements</h4>

TreeIO from BioPerl is required to use PUG.

<h4>Input</h4>

<b>Preparation of Input</b>:
	Each tree leaf needs to have a unique identifier for the species/accession associated with it. For example, a sequence from the species Andropogon virginicus might be--Androvirg.c1_g1_i1--while another from Andropogon gerardii would be--Androgera.c1_g1_i1. 

<b>Paralogs File</b>:
	A tab-delimited file of putative paralogs dervied from Ks analyses or synteny analysis. For each pair, a third column can be given that gives the source or the putative WGD event the use is testing.  If this is not supplied, then a default of "unknown" is used in downstream analysis. 

<b>Trees Directory</b>:
	A directory of tree files with bootstrap values. 

<b>Outgroups</b>:
	A comma-delmited list of outgroups where at least one must be in the gene tree for it to be used. These are needed to aid in accurate rooting of the phylogeny for the PUG search.

<b>Species Tree</b>:
	A species tree where the leaves are named so that they can be found in the gene trees. Following the example in "Preparation of Input", the species tree leaves for Andropogon virginicus and Andropogon gerardi would be Androvirg and Androgera, respectively.

<b>Name</b>:
	An optional input that allows you to have a run name prepended to all the outfiles.

<h4>Output</h4>

<b>Labeled_Species_Tree</b>:
        Newick formatted file containing user submitted species tree with internal nodes labeled to match those of the summary file.

<b>Labeled_Species_Tree.eps</b>:
	Postscript file containing user submitted species tree with internal nodes labeled to match those of the summary file.

<b>Gene_Tree_Results.txt</b>:
	Summary of paralog pairs and their relative placment in gene tree. Will finish this description to define each column in file.

<b>Gene_Trees_Pairs_Bad_Results.txt</b>:
	Summary of paralog pairs that did not fit any placement based on the species tree.

<b>Paralog_Pairs_Per_Tree.txt</b>:
	Gene tree files names, the paralog pairs identified in the gene tree, and the bootstrap value of the LCA of the paralogs.

<b>Summary_Results.txt</b>:
	Summary of all positive results per putative polyploid event, per labeled species tree node, and at each possible bootstrap value.  Further filtering of these results can 	  easily be done by the user.

<b>Paralog_Pairs_Nodes_Bootstraps.txt</b>:
	Filtered list of paralog pairs that map to a specific node, as marked in the labeled species tree. Orthogroup, node, bootstrap value, and pair names are given.

<b>SUMMARIZED_5080.txt</b>:
	Summary of mapped gene duplication by node.  Counts are given for thresholds of bootstrap values 80 and 50. These results are used to make the figures in the PUG_Figure_Maker.R script.

<b>Estimated_Putative_Paralogs.txt</b>:
	If the "estimate_paralogs" option is selected, all possible unique gene pairs from gene trees are identified and printed to this file.  Format is that of the typical input paralog file for PUG.

<h4>Usage</h4>

perl PUG.pl --paralogs file --trees directory --outgroups list --species species_tree [options]

<b>Options:</b>
         
	-paralogs    		Note: If this file is not provided then PUG will default with the "estimate_paralogs" option. Tab delmited file of putative in-paralogs used to query gene trees for placement of WGD event. The file is expected to have paralogs in the first two columns, and a third that has the putative event name. This allows the user to use paralog sets derived from multiple sources in the same analyses. It also also the user to have multiple events in the same analysis.
	-trees    			Directory of gene trees used to identify WGD placement. Trees should be in Newick format and have bootstrap values. The "bipartitions" file of RAxML works well.
	-outgroups     		Comma delimited list of outgroups for trees. At least one of these taxa must be in the tree for it to be considered.
	-species_tree     	Newick format tree with species relationships of all taxa present in gene trees. The names of these taxa must be represented in the names of the genes in the gene trees and in the paralogs file. Make sure they are unique and do not overlap with any other portion of the gene names.
	-name    			Identifier for this run. This is be a prefix to all output files.  [Default = "PUG"]
	-all_pairs			Flag that allows for counting of all paralog pairs, not just unique LCAs.
	-estimate_paralogs	Estimates all possible unique gene pairs from gene trees to test with PUG algorithm. This is now the default ifa paralog file is not provided.
	-help    			Brief help message
	-man    			Full documentation


<h4>Figure Maker</h4>

A script for making publication quality figures of a PUG analysis is now provided. These figures are similar to those found in <a href="http://gbe.oxfordjournals.org/content/8/4/1150.full?sid=9825eb8d-7200-4945-8671-a3224f4e4f29">McKain et. al. 2016</a> and an example is provided below.

<h4>Usage</h4>

Rscript PUG_Figure_Maker.R Labeled_Species_Tree.tre SUMMARIZED_5080.txt New_File_Name.pdf

<b>Example of PUG_Figure_Maker output.</b><br>
<img src="Final_Orchid_WGD.pdf" />





