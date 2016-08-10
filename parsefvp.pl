#!/usr/bin/perl -w

if($#ARGV==-1) {
    print "\tFragment score\tRead score\tlen50\tlen100\tlen150\tlen500\tlen1000\tlen2000\tlen3000\tlen10000\tlen100000\tlen1000000";
    exit;
}

open(File, $ARGV[0]) ||die "error: can't open $ARGV[0].\n";

my $n=0;
while(<File>){
    next if($_ eq "\n");
    chomp;
    
    my @clm= split(/\t/, $_);
    $n++;   

    next if($n==3);
    last if($n>12);
    if($clm[0] =~ /insufficient/) {
	print "($clm[1])\t";
    } else {
	print "$clm[1]\t";
    }
}
close (File);
