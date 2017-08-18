#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
die "segment4bam.pl <file> <segment length>\n" if($#ARGV !=1);

my $file=$ARGV[0];
my $segsize=$ARGV[1];

open(FILE, "samtools view -F 0x04 -b $file | samtools view |");
while (<FILE>) {
    chomp;
    my @c = split(/\t/, $_);
    my $posi = int($c[3]/$segsize);
    if($posi%2) {
	print "$_\n";
    }
}
close FILE;
