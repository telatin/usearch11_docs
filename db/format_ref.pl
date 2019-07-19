#!/usr/bin/env perl
# ABSTRACT: A program to format GreenGenes/SILVA FASTA + taxonomy in USEARCH compatible DB

use 5.012;
use warnings;
use FASTX::Reader;
my ($fa, $tax) = @ARGV;

open my $TAX, '<', "$tax" || die "Unable to open taxonomy $tax\n";
my %annot = ();
while (my $line = readline($TAX) ) {
        chomp($line);
        my ($id, $taxonomy) = split /\t/, $line;
        $taxonomy =~s/k__/d__/g;
        $taxonomy =~s/\s//g;
        $taxonomy =~s/__/:/g;
        $taxonomy =~s/;/,/g;
        $annot{$id} = $taxonomy;
}
print STDERR "Taxonomy parsed ($tax)\n";

my $f = FASTX::Reader->new({ filename => $fa });
while ( my $s= $f->getRead() ) {
        die "Unable to get annotation for $s->{name}\n" if (not defined $annot{ $s->{name} });
        my $name = $s->{name} . ';tax=' . $annot{ $s->{name} } . ";";
        print ">$name\n", $s->{seq}, "\n";
}
__END__
Desired format:
>AB008314;tax=d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus;

INPUT 1: Silva Taxonomy:
GY203941         k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella 7; s__?
GY324971         k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Photobacterium; s__?
JQ771130         k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__Acinetobacter nectaris

INPUT 2: Silva FASTA:
>GY203941
CTCAGGATGAACGCTGGCTACAGGCTTAACACATGCAAGTCGAGGGGAAACGGCATTGAGTGCTTGCACTGAATGGACGT
