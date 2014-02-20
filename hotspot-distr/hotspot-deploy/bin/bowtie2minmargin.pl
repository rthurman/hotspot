#!/usr/bin/perl -w
use strict;

# Custom bowtie output filter for uniqueness, instead of using -m 1

# maybe we want to exclude reads with a 1-mismatch-worse next-best mapping
# that could cause trouble at a SNP
# So throw out a 0-mismatch read if it has a 1-mismatch position somewhere else,
#    throw out a 1-mismatch read if it has a 2-mismatch position somewhere else.

#FC3014N-1:2:3:408	+	chr14	30796590	CAGGCTGGAGTGCAATGGCGTGATCTC	BBBBBBBBBBBBBBBBBBBABBBBABB	2068	
#FC3014N-1:2:3:408	+	chr1	53012609	CAGGCTGGAGTGCAATGGCGTGATCTC	BBBBBBBBBBBBBBBBBBBABBBBABB	2068	
#FC3014N-1:2:3:528	+	chr10	29462672	CCACCACCTCCCAGGTTCGAGTGGTTC	97B<ABABABBBCABBBBBBB47A7=B	0	
#FC3014N-1:2:3:528	+	chr16	88720266	CCACCACCTCCCAGGTTCGAGTGGTTC	97B<ABABABBBCABBBBBBB47A7=B	1	18:A>G

my $prevName = undef;
my @matches = ();

# Depends on bowtie parameters.  Can be fixed after the fact if set too high,
# output will lack results of (maxDetectableMargin - 1)
my $maxDetectableMargin = 3; # maybe 3 or 2?
if (@ARGV > 0) {
    $maxDetectableMargin = 0 + $ARGV[0];
}

while (<STDIN>) {
    chomp;
    my $line = $_;
    my $mapping = parseBowtie( $line );
    if ((not defined $prevName) or $prevName ne $mapping->{"readName"}) {
        if (defined $prevName) {
            printMinMargin();
        }
        $prevName = $mapping->{"readName"};
    }
    push @matches, $mapping;
}
if (defined $prevName) {
    printMinMargin();
}

# Uses global @matches containing one or two reads
sub printMinMargin {
    my $best = $matches[0];
    #warn $best->{"line"} . "\n";
    my $margin = 0;
    if (@matches == 1) {
        $margin = $maxDetectableMargin;
    } else {
        my $lesser = $matches[1];
        my $bestMM = $best->{"mismatches"};
        my $moreMM = $lesser->{"mismatches"};
        $margin = $moreMM - $bestMM;
    }
    my ($chr,$pos) = split( /:/, $prevName );
    my ($min,$max) = ($pos,$pos+1);
    my $line = "$chr\t$min\t$max\t$margin\n";
    print $line;
    @matches = ();
    return $margin;
}

sub parseBowtie {
    my $line = shift;
    # reserved field may be number of equally-good alternate positions
    my ($readName,$strand,$chromosome,$position,$sequence,$quality,$reserved,$substitutions) = split( /\t/, $line );
    my @substDescriptors = split( /,/, $substitutions );
    foreach my $descriptor (@substDescriptors) {
        if ($descriptor =~ /^([\d]+)\:([ACGTN])\>([ACGTN])$/) {
            my ($readOffset,$refBase,$readBase) = ($1,$2,$3);
            #warn "$refBase/$readBase at $readOffset\n";
        } else {
            die "Failed to parse ($descriptor)\n";
        }
    }
    my $mismatches = @substDescriptors;
    my %hash = ( 
        line => $line,
        readName => $readName,
        mismatches => $mismatches,
        chromosome => $chromosome,
        position => $position,
        strand => $strand,
        sequence => $sequence
    );
    return \%hash;
}

