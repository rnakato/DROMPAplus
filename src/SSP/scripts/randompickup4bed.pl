#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
die "random_pickup4bam.pl <file> <num>\n" if($#ARGV !=1);

my $file=$ARGV[0];
my $num=$ARGV[1];
my $wc = `wc -l $file`;
my @c = split(/ /, $wc);
my $nread = $c[0];
my $p = $num/$nread;

open(FILE, $file) ||die "error: can't open $file.\n";
while (<FILE>) {
    my $x = rand();
    if($x < $p){ print $_;}
}
close FILE;
