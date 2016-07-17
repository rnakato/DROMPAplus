#!/usr/bin/perl -w

if($#ARGV==-1) {
    print "\tNSC\tRLSC\tfragment length\tBackground enrichment\tBackground uniformity\n";
    exit;
}

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
