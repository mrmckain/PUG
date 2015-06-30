PUG--Phylogenetic Placement of Polyploidy Using Genomes
=============
<b>Author</b>: Michael R. McKain<br>
</br>
Version 2.0, June 30, 2015
<br></br>
<b>Contact</b>: https://github.com/mrmckain
<h3>Description</h3>

PUG is designed to test the phylogenetic placment of a polyploid event over mulitple gene trees using a known species tree and a set of putative paralogs.
Despite the name, PUG works with genomic or transcriptomic data. 

Paralogs are queried in gene trees. Their placement, and the relative placement of taxa in the gene tree, are compared to the given species tree and if they are found to match a node of the species tree, the pair is counted towards supporting the polyploid event origin at that node. A single node may be represented in a gene tree by multiple paralog pairs, but that node is counted only once for the gene tree.

<h4>Requirements</h4>

TreeIO from BioPerl is required to use PUG.

<h4>Output</h4>

<b>Labeled_Species_Tree</b>:<br></br>
        Newick formatted file containing user submitted species tree with internal nodes labeled to match those of the summary file.

<b>Labeled_Species_Tree.eps</b>:<br></br>
	Postscript file containing user submitted species tree with internal nodes labeled to match those of the summary file.

<b>Gene_Tree_Results.txt</b>:<br></br>
	Summary of paralog pairs and their relative placment in gene tree. Will finish this description to define each column in file.

<b>Gene_Trees_Pairs_Bad_Results.txt</b>:<br></br>
	Summary of paralog pairs that did not fit any placement based on the species tree.

<b>Paralog_Pairs_Per_Tree.txt</b>:<br></br>
	Gene tree files names, the paralog pairs identified in the gene tree, and the bootstrap value of the LCA of the paralogs.

<b>Summary_Results.txt</b>:<br></br>
	Summary of all positive results per putative polyploid event, per labeled species tree node, and at each possible bootstrap value.  Further filtering of these results can 	  easily be done by the user.

<h4>Usage</h4>

perl PUG.pl --paralogs file --trees file --outgroups list --species species_tree [options]

<b>Options:</b>
         
	-paralogs    		Tab delmited file of putative in-paralogs used to query gene trees for placement of WGD event. The file is expected to have paralogs in the first two columns, and a third that has the putative event name. This allows the user to use paralog sets derived from multiple sources in the same analyses. It also also the user to have multiple events in the same analysis.
	-trees    			Directory of gene trees used to identify WGD placement. Trees should be in Newick format and have bootstrap values. The "bipartitions" file of RAxML works well.
	-outgroups     		Comma delimited list of outgroups for trees. At least one of these taxa must be in the tree for it to be considered.
	-species_tree     	Newick format tree with species relationships of all taxa present in gene trees. The names of these taxa must be represented in the names of the genes in the gene trees and in the paralogs file. Make sure they are unique and do not overlap with any other portion of the gene names.
	-name    			Identifier for this run. This is be a prefix to all output files.  [Default = "PUG"]
	-help    			Brief help message
	-man    			Full documentation



