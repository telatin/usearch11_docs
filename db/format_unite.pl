#!/usr/bin/env perl
# ABSTRACT: A script to format UNITE reference in USEARCH format


use 5.012;
use warnings;
use FASTX::Reader;
my ($fa) = @ARGV;

my $f = FASTX::Reader->new({ filename => $fa });
while ( my $s= $f->getRead() ) {
        my ($species, $id, $sh, undef, $tax) = split /\|/, $s->{name};
        $tax =~s/__/:/g;
        $tax =~s/;/,/g;
        my $name = "$id:$sh;tax=$tax";
        print ">$name\n", $s->{seq}, "\n";
}
__END__
GOAL:
>AB008314;tax=d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus;

INPUT (UNITE):
>Symbiotaphrina_buchneri|DQ248313|SH1641879.08FU|reps|k__Fungi;p__Ascomycota;c__Xylonomycetes;o__Symbiotaphrinales;f__Symbiotaphrinaceae;g__Symbiotaphrina;s__Symbiotaphrina_buchneri
