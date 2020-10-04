#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;
my $usage = q/Usage:
   get_mcounts.pl <run_folder>
   e.g. get_mcounts.pl R19209
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
# --
my $smpdir=shift(@ARGV) || die("Error: no input folder given!\n");
{ local $/='/';chomp($smpdir); }
my ($smp)=($smpdir=~m/([^\/]+)$/);
my $fbw="$smpdir/quant_r/reads_$smp.collapsed.fa_mapped2mature.bwt";
my $fr="$smpdir/reads_$smp.collapsed.fa";
die("Error: cannot find $fr\n") unless -f $fr;
die("Error: cannot find $fbw\n") unless -f $fbw;
my $numreads=0; #after expanding collapsed
my $ncr=0; #number of collapsed reads
open(FA, $fr) || die("Error opening file $fr\n");
while (<FA>) {
 if (m/^>(\S+)/) {
   my $n=$1;
   $ncr++;
   my ($x)=($n=~m/_x(\d+)$/);
   if ($x==0) { die("Error parsing multiplicity for read $x\n"); }
   $numreads+=$x;
 }
}
close(FA);

print STDERR 'Total reads : '.sprintf('%10d',$numreads)." ($ncr when collapsed)\n";
my $rmapped=0;
my (@rl, @ml, %rh, %mh);
open(BW, $fbw) || die("Error opening file $fbw\n");
while (<BW>) {
 chomp;
 my ($r, $strand, $m, $pos, $seq, $q, $mism)=split(/\t/);
 my $rhd=$rh{$r};
 if (++$rh{$r} == 1) { #first time seeing this read
   my ($x)=($r=~m/_x(\d+)$/);
   if ($x==0) { die("Error parsing multiplicity for read $x\n"); }
   $rmapped+=$x;
   #$rhd=[];
   push(@rl, $r); # @rl will be list of mapped reads
   #$rh{$r}=$rhd; 
 }
 #push(@$rhd, $m); # $rh{$r}= [$m1, $m2, ...]
 my $mhd=$mh{$m};
 unless ($mhd) {
   $mhd=[];
   push(@ml, $m); # @ml will be a list of potentially expressed mRNAs 
   $mh{$m}=$mhd;
 }
 push(@$mhd, $r); # $mh{$m}= [$r1, $r2, ...]
}
close(BW);
my $pm=sprintf('%.2f%%', ($rmapped*100.00)/$numreads);
print STDERR "     Mapped : ".sprintf('%10d',$rmapped)." ($pm)\n";
print STDERR ' ('.scalar(@rl). " collapsed reads mapped)\n";
open(MC, ">$smpdir/${smp}.mcounts.tsv") || die ("Error creating ${smp}.mcounts.tsv\n");
print MC "mature_miRNA\trcount\n";
foreach my $m (@ml) {
   my $mhd=$mh{$m}; #set of reads mapped to $m
   my $mrcount=0.0; #mapped reads for $m
   foreach my $r (@$mhd) {
     # is read multimapped?
     my $df=$rh{$r};
     my ($x)=($r=~m/_x(\d+)$/);
     if ($x==0) { die("Error parsing multiplicity for read $x\n"); }
     $mrcount += $x/$df;
   }
   print MC "$m\t".sprintf('%.2f', $mrcount)."\n";
}
close(MC);
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************

