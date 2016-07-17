#!/usr/bin/perl -w

open(File, $ARGV[0]) ||die "error: can't open $ARGV[0].\n";
while(<File>){
    next if($_ eq "\n");
    last if($_ =~ /Strand shift/);
    chomp;
    my @clm= split(/\t/, $_);
    print "$clm[1]\t";
}
close (File);

print "\n";
