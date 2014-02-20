#!/usr/bin/perl
use strict;

if (@ARGV < 2 ) {
    die "Usage:  $0 read_length path_to_genomic_bowtie_index chromosome(s).fasta\n"
    . "Can be run on a whole genome fasta, or in parallel on single chromosomes\n";
}

my ($merSize, $bowtieRefSeq, $fastaFile ) = @ARGV;
chomp( my $binDir = `dirname $0` );

my $testExistBowtieIndexFile = "$bowtieRefSeq.1.ebwt";
unless (-f $testExistBowtieIndexFile) {
    die "Failed to find bowtie index file $testExistBowtieIndexFile\n"; 
}

# Fasta format
my $READFORMAT='-f';
# Params to find top two matches.  Further filtering will decide if the top one stands alone.
my $UNIQUENESSMARGIN=1;
my $MATCHOPTIONS=" -n $UNIQUENESSMARGIN -v $UNIQUENESSMARGIN -k 2 ";

#  mersize-length FIFO queue for seeing if a previous uniquely mappable K-mer ends just prior to this index.
my @remembMer = ();
initFIFOqueue( $merSize );
my $prevIndex = -1;
my $prevChrom = undef;

open IN, "$binDir/chopFastaPseudoreads.pl $fastaFile $merSize "
     .  " | bowtie --mm $MATCHOPTIONS $READFORMAT $bowtieRefSeq - "
     .  " | $binDir/bowtie2minmargin.pl |" or die "$!";
while (<IN>) {
    chomp;
    my ($chromosome,$min0,$max1,$code) = split /\t/;
    if (defined($prevChrom) and ($prevChrom ne $chromosome)) {
        # fill to the end, in case there's remembMer contains non-zero entries
        zeroFillToEnd( $prevChrom, $merSize, $prevIndex);
        # reset the queue
        initFIFOqueue( $merSize );
        $prevIndex = -1;
    }
    $prevChrom = $chromosome;
    # fill in zeroes up to this point
    while ($prevIndex < $min0 - 1) {
        # read a 0
        ++$prevIndex;
        processUniquenessCode( $chromosome, $prevIndex, 0 );
    }
    # now do this input line
    for (my $i = $min0; $i < $max1; ++$i) {
        if ($code > 0) {
            $code = 1; # either 0 non-unique or 1 unique
        }
        processUniquenessCode( $chromosome, $i, $code );
        $prevIndex = $i;
    }
}
close IN;
zeroFillToEnd( $prevChrom, $merSize, $prevIndex);

# a mapped K-mer can indicate a cut on either end depending on which genomic strand it maps to 
sub initFIFOqueue {
    my $merSize = shift;
    @remembMer = ();
    for (my $i = 0; $i < $merSize; ++$i) { 
        push( @remembMer, 0 );
    }
}

# Fill to the end to get the last reverse-strand matches, in case the queue contains non-zero entries
# Goes just one bp shy of the end, to keep tag coordinates within chromosome coordinates 
sub zeroFillToEnd {
    my ($chromosome,$merSize,$prevIndex) = @_;
    #for (my $i = 0; $i < $merSize; ++$i) {
    for (my $i = 0; $i < $merSize - 1; ++$i) {
        ++$prevIndex;
        processUniquenessCode( $chromosome, $prevIndex, 0 );
    }
}

## Process one position.
#    index:  0-based chromosome index
#    code:   0 for non-unique, 1 for uniquely-mappable mer starting at this position 
sub processUniquenessCode {
    my ($chromosome,$index,$code) = @_;
    my $prevMerCode = shift( @remembMer );
    push( @remembMer, $code ); 
    my $mappableLevel = $code + $prevMerCode;
    if ($mappableLevel > 0) {
        print join("\t", ($chromosome,$index,$index+1) ) . "\n";
    }
}

