#!env perl
use strict;
use warnings;
use Bio::SeqIO;
my $csv = 'strain_names_lookup.csv';
my $dir = shift || 'marker_aln';
open(my $fh => $csv) || die $!;
my $i =0;
my $line = <$fh>;
chomp($line);
my %header  = map { $_ => $i++ } split(/,/,$line);

my %lookup;
while(<$fh>) {
    chomp;
    my @row = split(/,/,$_);
    my ($sn, $nrrl, $sp) = map { $row [ $header{$_} ] } ('Sample ID',
							 'NRRL','Name');
    $sp =~ s/\s+$//;
    $sp =~ s/\s+/_/g;
    $lookup{$sn} = sprintf("%s_NRRL%s",$sp,$nrrl);
}

opendir(DIR,$dir) || die $!;
for my $file( readdir(DIR) ) {
    next unless $file =~ /\.(cds|aa)$/;
    my $in = Bio::SeqIO->new(-format => 'fasta',
			     -file   => "$dir/$file");

    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$dir/$file.rename");
    while( my $seq = $in->next_seq ) {
	my $id = $seq->id;
	$id=~ s/\.([^\.]+)\.scaffold//;
	if( $lookup{$id} ) {
	    $id = $lookup{$id};
	} else {
	    warn("no $id found\n");
	}
	$seq->id($id);
	$out->write_seq($seq);
    }
}
