#!/usr/bin/perl
use strict;
my $USAGE = q/Usage:
   get_mcounts.pl <subdir>..
   e.g. get_mcounts.pl R*
/;
umask 0002;

die("$USAGE\n") unless @ARGV>0;
my @mirs;
my $fxcsv='all_samples.tsv';
my $fxcsvnorm='all_samples.norm.tsv';
my $fstats_all='all_samples.stats.tsv';

open(XCSV, ">$fxcsv") || die("Error creating file $fxcsv\n");
open(XCSVNORM, ">$fxcsvnorm") || die("Error creating file $fxcsvnorm\n");
open(STATS, ">$fstats_all") || die("Error creating file $fstats_all\n");
print STATS "sample\ttotal_reads\tmapped\n";

# --
#my $smpdir=shift(@ARGV) || die("Error: no input folder given!\n");
#{ local $/='/';chomp($smpdir); }
foreach my $d (@ARGV) {
  { local $/='/';chomp($d); }
  my ($smp)=($d=~m/([^\/]+)$/);
  my %mircov; #mature miRNA => read count
  my $fbw="$d/quant_r/reads_$smp.collapsed.fa_mapped2mature.bwt";
  my $fr="$d/reads_$smp.collapsed.fa";
  die("Error: cannot find $fr\n") unless -f $fr;
  die("Error: cannot find $fbw\n") unless -f $fbw;
  if (@mirs==0) {
    my $fpre="$d/quant_r/mature.converted";
    open(PRE, $fpre) || die("Error opening $fpre!\n");
    while (<PRE>) {
      if (m/^>(\S+)/) {
        push(@mirs, $1);
      }
    }
    close(PRE);
    print XCSV "sample\t".join("\t",@mirs)."\n";
    print XCSVNORM "sample\t".join("\t",@mirs)."\n";
  }

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

  #print STDERR 'Total reads : '.sprintf('%10d',$numreads)." ($ncr when collapsed)\n";
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
  #print STDERR "     Mapped : ".sprintf('%10d',$rmapped)." ($pm)\n";
  #print STDERR ' ('.scalar(@rl). " collapsed reads mapped)\n";
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
     $mircov{$m}=sprintf('%.2f',$mrcount);
  }
  open(MC, ">$d/${smp}.mcounts.tsv") || die ("Error creating ${smp}.mcounts.tsv\n");
  print MC "mature_miRNA\t$smp\t${smp}_norm\n";
  #print MC "$m\t".sprintf('%.2f', $mrcount)."\t".\n";
  my @mv;
  my @mvn;
  foreach my $mm (@mirs) {
    my $v= $mircov{$mm};
    my $vn= (length($v)==0) ? 0.00 : sprintf('%.2f', ((1000000.00*$v)/$rmapped));
    print MC join("\t", $mm, $v, $vn )."\n";
    push(@mv, $v);
    push(@mvn, $vn);
  }
  close(MC);
  print XCSV "$smp\t".join("\t", @mv)."\n";
  print XCSVNORM "$smp\t".join("\t", @mvn)."\n";
  open(ST, ">$smp/$smp.mcounts.stats")  || die("Error creating file $d/$smp.stats\n");
  print ST "sample\ttotal_reads\tmapped\n";
  print ST join("\t", $smp, $numreads, $rmapped)."\n";
  close(ST);
  print STATS join("\t", $smp, $numreads, $rmapped)."\n";

}
close(XCSV);
close(XCSVNORM);
close(STATS);

