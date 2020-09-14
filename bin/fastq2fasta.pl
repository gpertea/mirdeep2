#!/usr/bin/perl

# miRDeep2 FASTQ-to-FASTA perl script
# Copyright (C) 2008 - 2011  Marc Friedländer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

use strict;

my $file = $ARGV[0] or die "Usage: \n\n $0 fastq file > outputfile\n";;

#open IN,"<$file" or die "File $file not found error here\n";
my $processed = 0;
my $fh=open_fastq($file) || die "File $file error!\n";;
while(1){
  my ($rname, $rseq, $rquals)=getFastq($fh);
  last unless $rname;
  $processed++;
  print ">seq_$processed\n$rseq\n";
}
close($fh);
exit;
#print STDERR "$processed reads processed\n\n";

sub open_fastq {
 my ($fn)=@_;
 my $funz;
 my $fh;
 if ($fn=~m/\.bzi?p?2$/) {
  $funz='bzip2';
 }
 elsif ($fn=~m/.g?zi?p?$/) {
  $funz='gzip';
 }
 if ($funz) {
  open($fh, $funz." -cd '".$fn."'|") || 
    die("Error creating decompression pipe: $funz -cd '$fn' !\n");
 }
 else {
    open($fh, $fn) || die("Error opening file $fn ! \n");
 }
 return $fh;
}

sub getFastq {
 my $fh=$_[0]; # file handle
 #parses next FASTQ record
 #returns ($readname, $readseq, $readquals)
 my ($rname, $rseq, $rquals);
 #while (<$fh>) {
 #  ($rname)=(m/^@(\S+)/);
 #  last if $rname;
 #}
 $_=<$fh>;
 return ($rname, $rseq, $rquals) if length($_)<3;
 ($rname)=(m/^@(\S+)/);
 die("Error reading FASTQ record:\n$_\n") unless $rname;
 while (<$fh>) {
   last if m/^\+/;
   chomp;
   $rseq.=$_;
   die("Error: sequence $rname too large?\n") if length($rseq>100000);
 }
 if ($_) {
   while (<$fh>) {
     chomp;
     $rquals.=$_;
     last if (length($rseq)<=length($rquals));
   }
   die("Error: seq and qual lengths differ for $rname!\n") 
     unless length($rseq)==length($rquals);
 }
 return ($rname, $rseq, $rquals);
}
