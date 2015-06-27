#!/usr/bin/perl -w
#use strict;
use Bio::TreeIO;
use Bio::Tree::Draw::Cladogram;
use Getopt::Long;
use Pod::Usage;

my $man = 0;
my $help = 0;
my $paralogs;
my $trees;
my $outgroups;
my $species_tree;
my $prefix = "GAPId";


GetOptions('help|?' => \$help, man => \$man, "paralogs=s" => \$paralogs, "trees=s" => \$trees, "species_tree=s" => \$species_tree, "name=s" => \$prefix) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


=head1 NAME
GAPId -- Genome Assisted Polyploidy Identification
=head1 SYNOPSIS
perl gapid.pl -paralogs file -trees file -outgroups list [options] 
=head1 OPTIONS
    -paralogs          Tab delmited file of putative in-paralogs used to query gene trees for placement of WGD event.
    -trees             Directory of gene trees used to identify WGD placement. Trees should be in Newick format and have bootstrap values. 
    -outgroups         Comma delimited list of outgroups for trees. At least one of these taxa must be in the tree for it to be considered.
    -species_tree      Newick format tree with species relationships of all taxa present in gene trees.
    -name              Identifier for this run.  [Default = "GAPId"]
    -help              Brief help message
    -man               Full documentation
=cut

my %hypotheses;
my %species_taxa;
my %species_index;
my %node_index;
my $node_count=0;
my $treeio = new Bio::TreeIO(-format => "newick", -file => "$species");
while( $tree = $treeio->next_tree ) {
        for my $node ( $tree->get_nodes){
                if($node->is_Leaf){
                    $species_taxa{$node->id}=1;
                    $species_index{$node->id}=$node;
                }
                else{
                        $node_count++;
                        $node->id = $node_count;
                        $node_index{$node}=$node_count;
                        $hypotheses{$node}=0;
                }
        }
    my $obj1 = Bio::Tree::Draw::Cladogram->new(-bootstrap => 0,
                                           -tree    => $tree,
                                           -compact => 1);
        $obj1->print(-file => "$prefix\_Tree.eps");
}



my $tree;
open my $OUT2, ">", "$prefix\_Paralog_Pairs_Per_Tree.txt";
open my $OUT, ">", "$prefix\_Gene\_Tree\_Initial\_Results.txt";
open my $BADOUT, ">", "$prefix\_Gene_Trees_Bad_Results.txt";
my %results; 
my %dups;
my %dups2;
my %seqpairs;
my %seqpairs2;
my %lcas;



###Load putative in-paralogs into hash for searching.###
open my $DUPS, "<", $paralogs; 
while (<$DUPS>) {
    chomp;
    my ($pair1, $pair2) = split /\s+/;
    $dups{$pair1}{$pair2} = "event";
    push(@{$dups2{$pair1}}, $pair2);
}
close $DUPS;

###Create outgroups data set.###
my @outs = split(/,/,$outgroups);

###Go through gene family trees to identify paralogs.###
my @file = <$trees/*>;
for my $file (@file){
    my $treeio = new Bio::TreeIO(-format => "newick", -file => "$file", -internal_node_id => 'bootstrap');
        while( $tree = $treeio->next_tree ) {
            my %pairs;
            my $what;
            my $tick=0;
            my %taxa;
            my $taxon;
			my $outexists=0;
            my @taxa = $tree->get_leaf_nodes;
            foreach my $taxa (@taxa){
                $taxon = $taxa->id;
				for my $outspecies (@outs){
                    if($taxon =~ /$outspecies/){
                        $outexists++;
                    }
                }
                $taxa{$taxon}=$taxa;   
            }
			if($outexists == 0){
				next;
			}
            for my $taxa (sort keys %taxa){
                if (exists $dups{$taxa}){
                    for my $possible_dup (@{$dups2{$taxa}}){
                        if (exists $taxa{$possible_dup} && $possible_dup ne $taxa){
                            my $which = $possible_dup;
                            $what = $dups{$taxa}{$which};
                            #$what is the event that is stored for the duplication pair
                            $what .= $tick unless $tick == 0;
                            $tick++;
                            push(@{$pairs{$what}}, $taxa{$taxa});
                            my $nodes = $tree->find_node(-id => "$possible_dup");
                            next unless (defined $nodes); 
                            my $check = $nodes->id;
                            push(@{$pairs{$what}}, $nodes);
                            $seqpairs{$what}=$taxa;
                            $seqpairs2{$what}=$which;
                        }
                    }
                }
            }
            for my $un (sort keys %pairs){
                if (scalar @{$pairs{$un}} < 2){
                    delete $pairs{$un};
                }
            }
            for my $dupeves (sort keys %pairs){
                my @sisters;
                my $lca = $tree->get_lca(-nodes => \@{$pairs{$dupeves}});
                my ($tax1_spec, $tax2_spec, $tax1_bs, $tax2_bs);
                my %tax_spec;
                my $counter=1;
                my $branch_length;
                for my $pair_tax(@{$pairs{$dupeves}}){
                    my ($anc, $old_anc);
                    $anc=$pair_tax;
                    $old_anc = $pair_tax;
                    while ($anc ne $lca){
                        $old_anc = $anc;
                        $anc = $old_anc->ancestor;
                    }
                    if($old_anc eq $pair_tax){
                        my $spec_leaf = $old_anc->id;
                        $tax_spec{$counter}{$pair_tax}{$spec_leaf}="single";
                        $counter++;
                    }
                    else{
                        my $temp_bs = $old_anc->bootstrap;
                        for my $child ($old_anc->get_all_Descendents){
                            if($child->is_Leaf){
                                my $spec_leaf = $child->id;
                                $tax_spec{$counter}{$pair_tax}{$spec_leaf}=$temp_bs;
                            }
                        }
                        $counter++;
                    }
                }
                for my $pair_taxid (@{$pairs{$dupeves}}){
                    if (exists $tax_spec{1}{$pair_taxid}){
                        for my $species_para_clade (sort keys %{$tax_spec{1}{$pair_taxid}}){
                            $tax1_spec .= $species_para_clade . ",";
                            $tax1_bs = $tax_spec{1}{$pair_taxid}{$species_para_clade};
                        }
                    }
                    elsif (exists $tax_spec{2}{$pair_taxid}){
                        for my $species_para_clade (sort keys %{$tax_spec{2}{$pair_taxid}}){
                            $tax2_spec .= $species_para_clade . ",";
                            $tax2_bs = $tax_spec{2}{$pair_taxid}{$species_para_clade};
                        }
                    }
                }
                my $tree_size = scalar keys %taxa;
                my $bs = $lca->bootstrap;
                my $above_taxa;
                my $outside_taxa;
                my $lca_used=0;
                my %counters;
                            
                my $children_concat;
                for my $child ( $lca->get_all_Descendents ) {
                    push(@sisters, $child);
                    $children_concat .= $child;
                }

				if (exists $lcas{$children_concat}){
                    $lca_used = "1";
                    next;
                }

			    my $file_tree = substr($file, index($file, ".")+1);
				$file_tree=substr($file_tree, 0, index($file_tree, "_"));
				print $OUT2 "$file_tree\t$bs\t";
				my ($pair_tax1, $pair_tax2);
				my $i=1;
				foreach my $pair_tax (@{$pairs{$dupeves}}){
					$pair_tax = $pair_tax->id;
					if($i==1){
						$pair_tax1=$pair_tax;
					}
					else{
						$pair_tax2=$pair_tax;
					}
					$i++;
					#print $OUT2 "\t$pair_tax";
				}
				print $OUT2 "$dups{$pair_tax1}{$pair_tax2}\t";
				$pair_tax1 =~ s/XXX/:/g;
				$pair_tax1 = substr($pair_tax1, 0,-1);
				$pair_tax2 =~ s/XXX/:/g;
				$pair_tax2 = substr($pair_tax2, 0,-1);
				print $OUT2 "$pair_tax1\t$pair_tax2\n";
	
                                	
				$lcas{$children_concat}=1;
                my $counting = -1;
                my %lca_species;
                for my $sister (@sisters){
                    if ($sister->is_Leaf){
					    $sister = $sister->id;
                        $above_taxa .= $sister . ",";
                        $lca_species{$sister}=0;
                    }
				}
                my $species_id;
        		$lca = $lca->ancestor;
        		$counting++;
        		for my $species ($lca->get_all_Descendents){
                	if ($species->is_Leaf){
                        $species_id = $species->id;
                    }
                	else {
                        next;
                    }
                	if (exists $lca_species{$species_id}){
                        next;
                    }
                    else{
					   $species = $species->id;
                        $outside_taxa .= $species . ",";
                    }                        
                }
                my $result = &hypothesis_test($outside_taxa,$above_taxa);
                if($result == 0){
                    print $BADOUT "$file\t$dupeves\t$bs\t$lca_used\t$pair_tax1\t$tax1_bs\t$tax1_spec\t$pair_tax2\t$tax2_bs\t$tax2_spec\t$outside_taxa\t$above_taxa\n";
                }
                else{
                    $hypotheses{$node_index{$result}}{$bs}++;
                    print $OUT "$file\t$dupeves\t$bs\t$lca_used\t$pair_tax1\t$tax1_bs\t$tax1_spec\t$pair_tax2\t$tax2_bs\t$tax2_spec\t$outside_taxa\t$above_taxa\n";
                }
            }
        }
    }
}

my %summary;
for my $node_index_id (keys %node_index){
    for (my $i = 0; $i <= 100; $i++){
        $summary{$node_index_id}{$i}=0;
    }
}

for my $hyp_result (keys %hypotheses){
    for my $bs_value (keys %{$hypotheses{$hyp_result}}){
        $summary{$hyp_result}{$bs_value}++
    }
}


open my $SUMMARY, ">", $prefix . "_Summary_Results.txt";
print $SUMMARY "Species_Tree_Node";
for (my $i = 100; $i >= 0; $i--){
    print $SUMMARY "\t$i";
}
print $SUMMARY "\n";

for my $hyp_result (sort {$a <=> $b} keys %summary){
    print $SUMMARY "$hyp_result";
    for my $bs_value (sort {$b <=> $a} keys %{$summary{$hyp_result}}){
        print $SUMMARY "\t$bs_value";
    }
    print $SUMMARY "\n";
}
###Build hypotheses to test from species tree.###
sub hypothesis_test {
    my %badin;
    my %badout;
    my %goodout;
    my %goodin;
    my @in_taxa = split(/,/,$_[1]);
    my @out_taxa = split(/,/,$_[0]);
    my %in_taxa;
    my %out_taxa;
    ###Get species names for genes found above LCA of paralogs and in the sister clade to the LCA.###
    for my $intax (@in_taxa){
            for my $spectax (keys %species_taxa){
                    if($intax =~ /$spectax/){
                            $in_taxa{$spectax}=1;
                    }
            }
    }

    for my $outtax (@out_taxa){
            for my $spectax (keys %species_taxa){
                    if($outtax =~ /$spectax/){
                            $out_taxa{$spectax}=1;
                    }
            }
    }    
    
    ###Check if out and in share taxa. If they do, then quit. These cannot fulfill any hypothesis.
    for my $outtax (keys %out_taxa){
            if(exists $in_taxa{$outtax}){
                my $value=0;
                return $value;
            }
    }

    ###Read in species tree to query nodes for WGD.###
    my $treeio = new Bio::TreeIO(-format => "newick", -file => "$species");
    while( $tree = $treeio->next_tree ) {
        my @taxa = $tree->get_leaf_nodes;
        my $root = $tree->get_root_node;
        my @current_in_taxa;
        for my $current_taxon (keys %in_taxa){
            push(@current_in_taxa, $species_index{$current_taxon});
        }
        ###Identifies the node of the species tree with all of the taxa sharing the putative WGD.###
        my $lca = get_lca(-nodes => \@current_in_taxa);
        ###Get all taxa that should be above LCA.###
        my @good_taxa = $lca->get_all_Descendents;
        my %good_taxa;
        for my $gtax (@good_taxa){
            if($gtax->is_Leaf){ 
                $good_taxa{$gtax->id}=1;
            }
        }

        ###Get all taxa sister to LCA.###
        my $ancestor_lca = $lca->ancestor;
        my @bad_taxa = $ancestor->get_all_Descendents;
        my %bad_taxa;
        my %sister_taxa;
        for my $btax (@bad_taxa){
            if($btax->is_Leaf){ 
                if(!exists $good_taxa{$btax->id}){
                    $sister_taxa{$btax->id};
                }
            }
        }
        ###Put all non-in taxa into a hash for query.###
        for my $all_taxa (@taxa){
            $all_taxa = $all_taxa->id;
            if(!exists $good_taxa{$all_taxa}){
                $bad_taxa{$all_taxa}=1;
            }
        }

        ###Check sister group taxa.  If there are none, then we are at the base of the tree and cannot accurately test the polyploid event position.###
        if(scalar keys %bad_taxa == 0){
            my $value = 0;
            return $value;
        }

        for my $temp_spec (keys %in_taxa){
                if(!exists $good_taxa{$temp_spec}){
                    my $value = 0;
                    return $value;
                }
        }
        ###Checking the taxa in the gene tree sister to the LCA of the paralogs, we first see if any of those taxa should be in the species tree LCA. If true, we skip this tree because it does't
        ###fit a hypothesis.  We then ask if there is one of the sister clade members from the species tree in the gene tree sister clade, if true we return the LCA node.
        for my $temp_spec (keys %out_taxa){
            if(exists $good_taxa{$temp_spec}){
                my $value = 0;
                return $value;
            }
        }
        for my $temp_spec (keys %out_taxa){
            if(exists $sister_taxa{$temp_spec}){
                my $value = $lca;
                return $value;
            }
        }

        my $value = 0;
        return $value

    }
} 
