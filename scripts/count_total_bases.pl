#!env perl
# count total number of bases in the datasets of reads
use strict;
use warnings;
use Bio::SeqIO;
my $csv = 'strain_names_lookup.csv';
my $dir = shift || '.';
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
print join("\t", qw(SPECIES NUM_READS TOTAL_BASES)),"\n";
for my $dir( readdir(DIR) ) {
    next unless $dir =~ /(\S+).PHYling/;
    my $id = $1;
#    if( $lookup{$id} ) {
#	warn($lookup{$id},"\n");
#    }
    open(my $fh => "seqstat $dir/$id.fasta | " ) || die $!;
    my ($n,$len);
    while(<$fh>) {
	if( /Total # residues:\s+(\d+)/ ) {
	    $len = $1;
	} elsif( /Number of sequences:\s+(\d+)/ ) {
	    $n = $1;
	}
    }
    print join("\t", $lookup{$id},$len,$n),"\n";
}
