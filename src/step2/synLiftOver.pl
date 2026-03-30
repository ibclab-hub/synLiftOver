#!/usr/bin/perl
#
#
use strict;
use warnings; 
use Parallel::ForkManager;
use FindBin qw($Bin);

my $f_in_bed = shift;
my $f_size = shift;
my $outdir = shift;
my $axtdir = shift;
my $minmatch = shift;
my $outputname = shift;
my $dup_allow = shift // 0; 
my $relax = shift // 0; 


my $split = $Bin."/splitBed.sh";
my $liftOver =  $Bin."/synLiftOver.out";

`mkdir -p $outdir/tmp/`;
#split
my $chrs = `bash $split $f_in_bed $outdir/tmp/`;
my @ar_chrs = split(/\n/,$chrs);


my $thread = scalar(@ar_chrs);
my $pm = new Parallel::ForkManager($thread);
foreach my $thischr (@ar_chrs){
	my $input_bed = "$outdir/tmp/$thischr.bed";
	my $output_bed = "$outdir/tmp/$thischr.out.bed";
	my $output_log = "$outdir/tmp/$thischr.log";
	my $axt = "$axtdir/$thischr.axt";
	$pm -> start and next;
	`$liftOver $input_bed $thischr $f_size $minmatch $axt $dup_allow $relax  1>$output_bed 2>$output_log  `;
	$pm -> finish;
}

$pm -> wait_all_children;

#merge all
`cat $outdir/tmp/*.out.bed > $outdir/$outputname`;


