#!env perl
use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::SeqIO;
my $db = shift || die "need a db";
my $dbh = Bio::DB::Fasta->new($db);
my %capture;
while(<>) {
    my ($sp,$acc,$gene) = split;
    push @{$capture{$gene}}, $dbh->get_Seq_by_acc($acc);
}
for my $gene ( keys %capture ) {
    my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$gene.hits.fas");
    $out->write_seq(@{$capture{$gene}});
}
