#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Temp qw(tempfile);
use POSIX qw(_exit);
@ARGV >= 6 && @ARGV <= 8 or die "Usage: synLiftOver.pl source.bed target.size outdir axtdir minmatch outputname [dup_allow] [relax]\n";
my ($bed,$size,$outdir,$axtdir,$minmatch,$outputname,$dup,$relax)=@ARGV;
$dup //= 0; $relax //= 0;
my $binary="$Bin/synLiftOver.out";
-r $bed or die "Cannot read source BED: $bed\n";
-r $size or die "Cannot read target sizes: $size\n";
-x $binary or die "Missing executable: $binary\n";
$outputname =~ m{^[^/]+$} && $outputname ne '.' && $outputname ne '..' or die "outputname must be a filename\n";
make_path("$outdir/tmp");
open my $src,'<',$bed or die "Cannot open $bed: $!\n";
my (%handles,%seen,@chrs);
while (my $line=<$src>) {
 next if $line =~ /^\s*$/;
 my ($chr)=split /\s+/,$line;
 $chr =~ /^[A-Za-z0-9_.-]+$/ && $chr ne '.' && $chr ne '..' or die "Unsupported chromosome name: $chr\n";
 if (!$seen{$chr}++) {
  push @chrs,$chr;
  open my $h,'>',"$outdir/tmp/$chr.bed" or die "Cannot create split BED: $!\n";
  $handles{$chr}=$h;
 }
 print {$handles{$chr}} $line or die "Cannot write split BED: $!\n";
}
close $src or die "Cannot close source BED: $!\n";
for my $h (values %handles) {close $h or die "Cannot close split BED: $!\n";}
@chrs=sort @chrs;
my @available;
for my $chr (@chrs) {
 my $axt="$axtdir/$chr.axt";
 if (!-e $axt) {
  my $message="WARNING: no AXT for source chromosome $chr ($axt); its query intervals are unmapped\n";
  warn $message;
  open my $log,'>',"$outdir/tmp/$chr.log" or die "Cannot create missing-AXT log: $!\n";
  print {$log} $message; close $log or die "Cannot close log: $!\n";
  open my $empty,'>',"$outdir/tmp/$chr.out.bed" or die "Cannot reset missing-AXT output: $!\n";
  close $empty or die "Cannot close empty output: $!\n";
  next;
 }
 -f $axt && -r $axt or die "AXT exists but is not a readable file: $axt\n";
 push @available,$chr;
}
@chrs=@available;
my %children;
for my $chr (@chrs) {
 my $pid=fork(); defined $pid or die "Cannot fork: $!\n";
 if ($pid==0) {
  open STDOUT,'>',"$outdir/tmp/$chr.out.bed" or die "Cannot create output: $!\n";
  open STDERR,'>',"$outdir/tmp/$chr.log" or die "Cannot create log: $!\n";
  my $status=system {$binary} $binary,"$outdir/tmp/$chr.bed",$chr,$size,$minmatch,"$axtdir/$chr.axt",$dup,$relax;
  if ($status!=0) {print STDERR "ERROR: liftover failed (system status=$status)\n"; _exit(1);}
  _exit(0);
 }
 $children{$pid}=$chr;
}
my @failed;
for my $pid (sort {$a<=>$b} keys %children) {
 my $r=waitpid($pid,0);
 push @failed,$children{$pid} if $r!=$pid || $?!=0;
}
@failed and die "Liftover failed for: @failed; final output was not replaced. Check tmp/*.log\n";
my ($dst,$temp)=tempfile('.synLiftOver-XXXXXX',DIR=>$outdir,UNLINK=>0);
for my $chr (@chrs) {
 open my $h,'<',"$outdir/tmp/$chr.out.bed" or die "Cannot read output: $!\n";
 while (my $line=<$h>) {print {$dst} $line or die "Cannot merge output: $!\n";}
 close $h or die "Cannot close output: $!\n";
}
close $dst or die "Cannot close merged output: $!\n";
rename $temp,"$outdir/$outputname" or die "Cannot publish output: $!\n";
