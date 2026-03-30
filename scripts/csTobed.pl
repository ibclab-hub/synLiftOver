#!/usr/bin/perl
#
#
use strict;
use warnings; 

my $f_in = shift;
my $ref  = shift;
my $tar = shift;
my $synid = 0;
open(F,"$f_in");
while(<F>){
	chomp;
	if ($_ =~ />(\d+)$/){
		$synid = $1;
		print $synid."\t";
	}else{
		if ($_ =~ /$ref\.(\S+):(\d+)-(\d+) (.)/){
			print $1."\t".$2."\t".$3."\t".$4."\t";
		}
		if ($_ =~ /$tar\.(\S+):(\d+)-(\d+) (.)/){
			print $1."\t".$2."\t".$3."\t".$4."\n";
		}
	}
}
close(F);
