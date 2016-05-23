#!/usr/bin/perl

use strict;
use warnings;

die "Specify t or c, threshold and version number!" unless (@ARGV==3);

my $type = $ARGV[0];
my $threshold = $ARGV[1];
my $version = $ARGV[2]; 

my %genes=();
for (my $i=1;$i<=100;$i++) {
	open IN, "<../data/nscumcomp_$type/genes/nscumcomp_genes.txt.$i" or die $!;
	while (<IN>) {
		chomp;
		$genes{$_}++;
	}
	close IN;
}
print "Total number of genes is ". scalar(keys %genes) . "\n";

my $count=0;
open OUT, ">../data/nscumcomp_$type/${threshold}TimesGenes.txt.$type.$version" or die $!;
foreach my $key(keys %genes) {
	if ($genes{$key} >= $threshold) {
		$count++;
		print OUT "$key\n";
	}
}
print "There are $count genes which appear $threshold times or more\n";

close OUT;
