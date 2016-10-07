#!/usr/bin/perl -w

if($#ARGV==-1) {
    print "\tFragment score\tRead score\tlen50\tlen100\tlen150\tlen500\tlen1000\tlen2000\tlen3000\tlen10000\tlen100000\tlen1000000";
    print "\tNSC\tRLSC\tfragment length\tBackground enrichment\tBackground uniformity\n";
    exit;
}

my $on=0;
open(File, $ARGV[0]) ||die "error: can't open $ARGV[0].\n";
while(<File>){
    next if($_ eq "\n");
    if($_ =~ /Strand shift/){
	$on=1;
	next;
    }
    chomp;
    my @clm= split(/\t/, $_);
    if(!$on){
	print "$clm[1]\t";
    } else {
	print "$clm[3]\t" if($clm[0] == -200 || $clm[0] == 0 || $clm[0] == 100 || $clm[0] == 150 || $clm[0] == 500 || $clm[0] == 1000 || $clm[0] == 1499);
    }
}
close (File);

print "\n";
