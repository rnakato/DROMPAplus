#!/usr/bin/perl -w

$file="";
if($ARGV[0] =~ /(.+)\.([0-9]+)\.tsv/){
    $file=$1;
}
if($file =~ /(.+)\/(.+)/){
    $file=$2;
}

$pcrbiasthre="";
$FRiP="";
$gcsummit="";
open(File, $ARGV[0]) || die "error: can't open $ARGV[0].\n";
while(<File>){
    next if($_ eq "\n");
    chomp;
    if($_ =~ /Redundancy threshold: (.+)/){
	$pcrbiasthre = $1;
    }elsif($_ =~ /Library complexity: (.+) \((.+)\)/){
	$tested_complexity = $1;
	$tested_reads = $2;
    }elsif($_ =~ /GC summit: (.+)/){
	$gcsummit = $1;
    }elsif($_ =~ /Genome/){
	chomp;
	my @clm= split(/\t/, $_);
	$total_reads = $clm[4];
	$plus = $clm[5];
	$minus = $clm[6];
	$total_remained = $clm[8];
	$total_filtered = $clm[11];
	if($gcsummit eq ""){
	    $depth = $clm[14];
	    $total_gc_base = $clm[17];
	    $nread_in_peak = $clm[19];
	    $FRiP = $clm[20];
	}else{
	    $depth = $clm[17];
	    $total_gc_base = $clm[20];
	    $nread_in_peak = $clm[22];
	    $FRiP = $clm[23];
	}
    }
}
close (File);

print "Sample\tMapped reads\t + strand\t - strand\tRedundancy threshold\tNonredundant\tRedundant\tComplexity for10M\tTested_reads\tRead depth\tGenome coverage";
print "\treads in peaks\tFRiP";
print "\tGC summit" if($gcsummit ne "");
print "\n";

print "$file\t$total_reads\t$plus\t$minus\t$pcrbiasthre\t$total_remained\t$total_filtered\t$tested_complexity\t$tested_reads\t$depth\t$total_gc_base";
print "\t$nread_in_peak\t$FRiP";
print "\t$gcsummit" if($gcsummit ne "");
print "\n";
