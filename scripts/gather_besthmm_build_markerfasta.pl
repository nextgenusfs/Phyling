#!env perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;
my $cutoff = 1e-40;
my @idirs;
my $force = 0;
my $odir = 'DB/markers/fungi/marker_files';
my $ext = 'domtbl';
my @dbfile; # = 'DB/genomes/fungi';
my $debug = 0;
my @cdsdbfile; # optional for the extraction of CDS
my $basename = 0;
GetOptions('db|pep:s'    => \@dbfile,
	   'cds:s'    => \@cdsdbfile,
	   'v|debug!' => \$debug,
	   'o|out:s'  => \$odir,
	   'i|in:s'   => \@idirs,
	   'ext:s'    => \$ext,
	   'basename!' => \$basename,
	   'cutoff|e|evalue:s' => \$cutoff,
	   'force!'   => \$force,
    );

my @dbs;

for my $n ( @dbfile ) {
  push @dbs, Bio::DB::Fasta->new($n);
}

mkdir($odir) unless -d $odir;
mkdir("$odir/aa") unless -d "$odir/aa";
mkdir("$odir/cds") unless -d "$odir/cds";

my @cdsdb;
if( @cdsdbfile ) {
    for my $n ( @cdsdbfile ) {
     push @cdsdb, Bio::DB::Fasta->new($n);
    }
}

my %hits_by_query;
for my $idir ( @idirs ) {
 opendir(DIR, $idir) || die "cannot open $idir: $!";
 for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.\Q$ext\E$/;
    my $p = $1;
    $p =~ s/\.TFASTX//;
    $p =~ s/\.v\d+//;
    my $hits = &get_best_hit_domtbl("$idir/$file");
    for my $h ( keys %$hits ) {
	$hits_by_query{$h}->{$p} = $hits->{$h};
    }
 }
}
for my $marker ( keys %hits_by_query ) {
    next unless ($force || ! -f "$odir/$marker.aa");
    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$odir/aa/$marker.aa");
    my $cdsout;
    if( @cdsdb ) {
	$cdsout = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">$odir/cds/$marker.cds");

    }
    while( my ($sp,$sdata) = each %{$hits_by_query{$marker}} ) {
	my ($sn,$sevalue,$hmmstart,$hmmend) = @$sdata;
	my $s;
	for my $db ( @dbs ) {
	  $s = $db->get_Seq_by_id($sn);
	  last if $s;
        }
	if( ! $s ) {
	    warn("cannot find $sn ($marker) in AA db \n");
	    next;
	}
	if( $basename ) {
	    $s = Bio::Seq->new(-display_id => $sp,
			       -seq        => $s->seq,
			       -desc       => sprintf("%s E=%s alnlen=%d hstart=%d hend=%d",$s->display_id,$sevalue,abs($hmmend-$hmmstart),$hmmstart,$hmmend));
	}

	$out->write_seq($s);
	if( @cdsdb ) {
	    for my $cdsdb ( @cdsdb ) {
	     $s = $cdsdb->get_Seq_by_id($sn);
	     last if $s;
            }
	    if( ! $s ) {
		warn("cannot find $sn ($marker) in CDS db\n");
		next;
	    }
	    if( $basename ) {
		$s = Bio::Seq->new(-display_id => $sp,
				   -seq        => $s->seq,
				   -desc       => $s->display_id);
	    }
	    $cdsout->write_seq($s);
	}
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
	warn("evalue = $evalue, start = ",$row[15], " end = ",$row[16],"\n") if $debug;
	$seen->{$q} = [$t,$evalue,$row[15],$row[16]];
    }
    close($fh);
    return $seen;
}
