#!env perl
use strict;
use warnings;
my $cutoff = 1e-50;
my %seen;
while(<>) {
    next if /^\#/;
    my @row = split(/\s+/,$_);
    my $query = $row[0];
    my ($sp,$q) = split(/\|/,$query);
    my $t = $row[3];
    my $evalue = $row[6];
    next if $evalue > $cutoff;
    if( exists $seen{$sp}->{$t} ) {
	next;
    }
    $seen{$sp}->{$t} = $query;
}
for my $sp ( keys %seen ) {
    for my $t ( keys %{$seen{$sp}} ) {
	print join("\t",$sp,$seen{$sp}->{$t},$t),"\n";
    }
}
