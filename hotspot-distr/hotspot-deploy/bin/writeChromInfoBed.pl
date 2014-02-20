#!/usr/bin/perl -w
use strict;

if (@ARGV < 1) {
    die "Usage:  $0  fasta file(s)\nWrites chromInfo.bed file listing chromosome sizes.\n";
}

my %chrom2size = ();

foreach my $fasta (@ARGV) {
    my $seqname  = "";
    my $sequence  = "";
    my $inContig  = 0;
    open IN, $fasta or die "$!";
    while (<IN>) {
        chomp;
        if (/^>\s*([^\s]+)/) {
            # new contig
            if ($inContig) {
                # flush out old one
                $chrom2size{$seqname} = length($sequence);
            }
            $inContig = 1;
            $seqname = $1;
            $sequence = "";
        } else {
            $sequence .= $_;
        }
    }
    if ($inContig) {
        # flush out last one
        $chrom2size{$seqname} = length($sequence);
    }
}

open OUT, "> chromInfo.bed" or die "$!";
foreach my $seqname (sort keys %chrom2size) {
    my $size = $chrom2size{$seqname};
    print OUT "$seqname\t0\t$size\t$seqname\n";
}
close OUT;

