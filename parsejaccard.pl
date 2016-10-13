#!/usr/bin/perl -w

use Getopt::Long qw/:config posix_default no_ignore_case bundling auto_help/;

my $nline=3;
GetOptions('nline|n=s' => \$nline);

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
    if(!$on && $nline==3){
	print "$clm[1]\t";
    } elsif($on) {
	print "$clm[$nline]\t" if($clm[0] == -200 || $clm[0] == 0 || $clm[0] == 100 || $clm[0] == 150 || $clm[0] == 500 || $clm[0] == 1000 || $clm[0] == 1499 || $clm[0] == 5000 || !($clm[0]%100000));
    }
}
close (File);

print "\n";
