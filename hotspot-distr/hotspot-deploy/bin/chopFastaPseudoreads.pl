#!/usr/bin/perl -w
use strict;

# Bowtie won't report hits with an N in the reference,
# so set this flag if mapping uniqueness with Bowtie.
my $exclude_reads_with_Ns = 1;

my $read_length = 36;
if (@ARGV < 1) {
    die "Usage:  chopFastaPseudoreads.pl <input fasta file> [mersize=$read_length]\n";
}
my $file = $ARGV[0];
if ($file =~ /.*\.gz$/) {
    open INFH, "zcat $file |" or die "Failed to read $file\n";
} else {
    open INFH, $file or die "Failed to read $file\n";
}
if (@ARGV > 1) {
    $read_length = 0 + $ARGV[1];
}
my $chromosome = "";
my $seq = "";
while (<INFH>) {
    chomp;
    my $line = $_;
    if (/^>(.*)/) {
        my $name = $1;
        if (length($seq) > 0) {
            writeFragments();
        }
        $chromosome = $name;
    } else {
        $seq .= uc( $line );
    }
}
writeFragments();
close INFH;

# Write read-length substrings of chromosome
sub writeFragments {
    my $len = length( $seq );
    for (my $i = 0; $i <= $len - $read_length; ++$i) {
        my $readSeq = substr( $seq, $i, $read_length );
        my $readName = "$chromosome:$i";
        if (($readSeq =~ /.*N.*/) and $exclude_reads_with_Ns) {
            #warn "Skipping $readSeq\n";
        } else {
            print ">$readName\n$readSeq\n";
        }
    }
    $seq = "";
}

