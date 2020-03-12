#!/usr/bin/perl -w

use IO::Uncompress::Gunzip;

if ($#ARGV != 1) {
    print " makemappabilitytable.pl <genometable> <gzipped binary mappability file prefix>\n";
    exit;
}

$gtfile = $ARGV[0];
$dir = $ARGV[1];

$num=0;
open(InputFile, $gtfile) || die "Error: can't open $gtfile.\n";
while(<InputFile>) {
    next if($_ eq "\n");
    chomp;
    @clm= split(/\t/, $_);
    $name[$num] = $clm[0];
    $len[$num] = $clm[1];
    $num++;
}
close (InputFile);

for ($i=0; $i<$num; $i++) {
    my $filename = "${dir}_${name[$i]}_binary.txt.gz";

    my $fh_in = IO::Uncompress::Gunzip->new($filename);
    open(InputFile, $filename) || die "Error: can't open $filename.\n";
    while(my $line = readline $fh_in) {
	$count1 = ($line =~ tr/1/1/);
    }
    if (!$count1) { die "Error: no mappable seq data in $filename.\n"; }
    my $r = $count1 / $len[$i];
    print "$name[$i]\t$count1\n";
    close (InputFile);
}
