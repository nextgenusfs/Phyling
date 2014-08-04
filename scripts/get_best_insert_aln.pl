#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Bio::SeqIO;
use Bio::DB::Fasta;

my @domtbl;
# = "orthologs-vs-nr_90.HMMSEARCH.domtbl";
my $pepdir = 'pep';
my $pepext = 'pep.fa';
my $odir = 'updated_seqs';
my $db; #='CLAT.nr_90.fasta';
my $cutoff = 1e-25;
GetOptions(
    'tbl:s'    => \@domtbl,
    'db:s'     => \$db,
    'o|out:s'  => \$odir,
    'c|cutoff:s' => \$cutoff,
    );

mkdir($odir) unless -d $odir;
my %add;
my $dbh = Bio::DB::Fasta->new($db);
if( ! $dbh ) { 
 warn "cannot open $db\n";
 exit;
}
for my $file ( @domtbl ) {
 my $last_model;
 open(my $fh => $file) || die $!;

 while(<$fh>) {
    next if /^\#/;
    my ($gene,$desc,$genelen,$model,$modeldesc,$modellen,$evalue) = split;
    next if $evalue > $cutoff;
    if( ! $last_model || $last_model ne $model) {
        $add{$model}->{$gene}++;
    }
    $last_model = $model;
}
for my $model ( keys %add ) {
	my $in = Bio::SeqIO->new(-format => 'fasta',
				     -file   => sprintf("%s/%s.%s",
							$pepdir,
							$model,
							$pepext));
	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => sprintf(">%s/%s.fa",$odir,$model));
	my %seqs;
	while( my $seq = $in->next_seq ) {
	    my ($sp) = split(/\|/,$seq->id);
	    if( $seqs{$sp} ) {
		if( $seqs{$sp}->length < $seq->length ) {
		    $seqs{$sp} = $seq;
		}
	    } else {
		$seqs{$sp} = $seq;
	    }
	}
	for my $s ( values %seqs  ) {
	    $out->write_seq($s);
	}
        for my $gene ( keys %{$add{$model}} ) {
	 my $s = $dbh->get_Seq_by_acc($gene);
	 if ( ! $s ) { warn("cannot find $gene\n"); next }
	 $out->write_seq($s);
        }
    }

}
