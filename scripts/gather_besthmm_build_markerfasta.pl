#!env perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;
my $cutoff = 1e-40;
my $idir = 'out';
my $odir = 'DB/markers/fungi/marker_files';
my $ext = 'domtbl';
my $dbfile = 'DB/genomes/fungi';
my $debug = 0;
GetOptions('db:s'     => \$dbfile,
	   'v|debug!' => \$debug,
	   'o|out:s'  => \$odir,
	   'i|in:s'   => \$idir,
	   'ext:s'    => \$ext,
    );
my $db = Bio::DB::Fasta->new($dbfile);

opendir(DIR, $idir) || die "cannot open $idir: $!";
my %hits_by_query;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.\Q$ext\E$/;
    my $p = $1;
    my $hits = &get_best_hit_domtbl("$idir/$file");
    for my $h ( keys %$hits ) {
	$hits_by_query{$h}->{$p} = $hits->{$h};
    }
}
for my $marker ( keys %hits_by_query ) {
    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$odir/$marker.fa");
    for my $sn ( values %{$hits_by_query{$marker}} ) {
	my $s = $db->get_Seq_by_acc($sn);
	if( ! $sn ) {
	    warn("cannot find $sn ($marker) \n");
	    next;
	}
	$out->write_seq($s);
    }
}

sub get_best_hit_domtbl {
    my $file = shift;
    open(my $fh => $file) || die $!;
    my $seen;
    while(<$fh>) {
	next if /^\#/;
	my @row = split(/\s+/,$_);
	my $t = $row[0];
	my $q = $row[3];
	my $evalue = $row[6];
	next if $evalue > $cutoff;
	if( exists $seen->{$q} ) {
	    next;
	}
	$seen->{$q} = $t;
    }
    close($fh);
    return $seen;
}
