#!/usr/bin/env perl -w

$filename = $ARGV[0];
$len=0;
$lenall=0;

open(InputFile,$filename) ||die "error: can't open file.\n";
while($line = <InputFile>){
    next if($line eq "\n");
    chomp $line;
    if($line =~ ">"){
	print "$len\n" if($len);
	$lenall += $len;
	$len=0;
	if($' =~ /([A-Za-z0-9_]+)\s(.+)/){ $name = $1;}
	else{ $name = $';}
	print "$name\t";
	next;
    }else{
	chomp($line);
	$len += length($line);
    }
} 
close (InputFile);

print "$len\n";
