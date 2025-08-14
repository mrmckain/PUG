#!/user/env perl
use strict;
use warnings;


my %summary;
open my $file, "<", $ARGV[0];
while(<$file>){
	chomp;

	my @tarray = split/\s+/;
	#unless($tarray[5] eq "single" && $tarray[8] eq "single"){
	#	next;
	#}

	my %species_count;
	my @set1 = split(/,/, $tarray[6]);
	my @set2 = split(/,/, $tarray[9]);

	for my $tax1 (@set1){
		$tax1 =~ /^(.*?_.*?)_/;
		$species_count{$1}++;
	}
	for my $tax1 (@set2){
		$tax1 =~ /^(.*?_.*?)_/;
		$species_count{$1}++;
	}
	unless (scalar keys %species_count == 1){
		next;
	}
	$tarray[4] =~ /^(.*?_.*?)_/;
	my $id1 = $1;

	$tarray[7] =~ /^(.*?_.*?)_/;
	my $id2 = $1;

	my $bsv = $tarray[2];
	if($bsv >= 50){
		$summary{$id1}{BSV50}++;
	}
	if($bsv >= 80){
		$summary{$id1}{BSV80}++;
	}
}


open my $out, ">", "Tip_dupes_summary.txt";
print $out "ID\tBSV_50\tBSV_80\n";

for my $id (keys %summary){
	print $out "$id\t$summary{$id}{BSV50}\t$summary{$id}{BSV80}\n";
}