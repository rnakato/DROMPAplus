#!/usr/bin/perl -w

if($#ARGV==-1) {
    print "\tNSC\tRLSC\tfragment length\tBackground enrichment\tBackground uniformity";
    print "\t# of reads in peak\t%\t# of reads in background\t%\t# of reads in repeats\t%\n";
    exit;
}

open(File, $ARGV[0]) ||die "error: can't open $ARGV[0].\n";

while(<File>){
    next if($_ eq "\n");
    last if($_ =~ /Strand shift/);
    chomp;
    my @clm= split(/\t/, $_);
    if($clm[0] =~ /nread/) {
	print "$clm[1]\t$clm[2]\t";
    }
    else {
	print "$clm[1]\t";
    }
}
close (File);

print "\n";
