#!/usr/bin/perl
use strict;
my $USAGE=q/Usage:
 get_data.pl <subdir>
/;

die("$USAGE\n") unless @ARGV>0;
my @mirs;
my $fxcsv='all_samples.xcsv';
my $fxcsvnorm='all_samples.norm.xcsv';
my $fstats_all='all_samples.stats';
open(XCSV, ">$fxcsv") || die("Error creating file $fxcsv\n");
open(XCSVNORM, ">$fxcsvnorm") || die("Error creating file $fxcsvnorm\n");
open(STATS, ">$fstats_all") || die("Error creating file $fstats_all\n");
print STATS "sample\ttotal_reads\tmapped\n";

foreach my $d (@ARGV) {
  next unless $d eq '.' || -d $d;
  die("Error: $d/quant_r does not exist!\n") unless -d "$d/quant_r";
  my $totalMapped=0;
  my %mircov; #precursor miRNA => read count
  my @fbw=<$d/quant_r/reads_*.collapsed.fa_mapped.bwt>;
  die("Error: cannot find the read-to-precursor mappings file!\n") unless @fbw>0;
  print STDERR "File is $fbw[0]\n";
  my ($smp)=($fbw[0]=~m/\/reads_(\w+)\./);
  die("Error cannot retrieve sample name!\n") unless $smp;
  print STDERR "Sample name is <$smp>\n";
  my $f="$d/quant_r/reads_${smp}.collapsed.fa.converted";
  open(FA, $f) || die("Error opening $f !\n");
  open(F, $fbw[0]) || die("Error opening $fbw[0]\n");
  my $fpre="$d/quant_r/precursor.converted";
  if (@mirs==0) {
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
  my $rMapped=0;
  my $totalReads=0;
  my %rds;
  while (<F>) {
    next if m/^#/;
    chomp;
    my @t=split(/\t/);
    my ($n)=($t[0]=~m/_x(\d+)$/);
    $mircov{$t[2]}+=$n;
    $totalMapped+=$n;
    next if $rds{$t[0]};
    $rds{$t[0]}=1;
    $rMapped+=$n;
  }
  close(F);

  while(<FA>) {
    if (m/^>(\S+)/) {
      $_=$1;chomp;
      my ($n)=(m/_x(\d+)$/);
      $totalReads+=$n;
    }
  }
  close(FA);
  #my @mlst = sort { $mircov{$b}<=>$mircov{$a} } (keys(%mircov));

  open(CSV, ">$d/$smp.csv") || die("Error creating file $d/$smp.csv\n");
  print CSV "pre_miRNA\t$smp\t$smp.norm\n";
  my @mv;
  my @mvn;
  foreach my $mm (@mirs) {
    #$totalMapped+=$mircov{$mm};
    my $v=int($mircov{$mm});
    my $vn= ($v==0) ? 0.00 : sprintf('%.2f', ((1000000.00*$v)/$totalMapped));
    print CSV join("\t", $mm, $v, $vn )."\n";
    push(@mv, $v);
    push(@mvn, $vn);
  }
  close(CSV);
  print XCSV "$smp\t".join("\t", @mv)."\n";
  print XCSVNORM "$smp\t".join("\t", @mvn)."\n";
  open(ST, ">$d/$smp.stats")  || die("Error creating file $d/$smp.stats\n");
  print ST "sample\ttotal_reads\tmapped\n";
  print ST join("\t", $smp, $totalReads, $rMapped)."\n";
  close(ST);
  print STATS join("\t", $smp, $totalReads, $rMapped)."\n";
}
close(XCSV);
close(XCSVNORM);
close(STATS);

#print STDERR "Total reads  =$totalReads\n";
#my $pm=sprintf('%.2f', (100.00*$rMapped)/$totalReads);
#print STDERR "Mapped reads =$rMapped (${pm}%)\n";
#print STDERR "totalMapped =$totalMapped\n";
