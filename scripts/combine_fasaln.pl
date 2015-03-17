#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Getopt::Long;
my $iformat = 'fasta';
my $oformat = 'nexus';
my $outfile = 'allseq.nex';
my $ext = 'fasaln.trim';
my $dir;
my @expected;
my $expected_file;
my $debug ;
my $skip;
GetOptions('d|dir:s'   => \$dir,
	   'ext:s'     => \$ext,
	   'if:s'       => \$iformat,
	   'of:s'       => \$oformat,
	   'expected:s' => \$expected_file,
	   'o|out:s'   => \$outfile,
	   'v|debug!' => \$debug,
           's|skip:s' => \$skip,
	   );

die("need a dir") unless $dir && -d $dir;
my %skip_aln;
if( $skip ) {
 open(my $fh => $skip ) || die $!;
 while(<$fh>) { 
   my ($id) = split;
   $skip_aln{$id}++;
  }
}
opendir(DIR, $dir) || die"$dir: $!";

if( $expected_file && open(my $fh => $expected_file) ) {
    while(<$fh>) {
	chomp;
	push @expected, $_;
    }
}
my (%matrix);

for my $file (sort readdir(DIR) ) {
    next if $file eq $outfile;
    next unless ($file =~ /(\S+)\.\Q$ext\E$/);
    my $stem = $1;
    if( $skip_aln{$stem} ) {
	warn("skipping $stem\n");
	next;
    }
    my $in = Bio::AlignIO->new(-format => $iformat,
			       -file   => "$dir/$file");
    warn($file,"\n") if $debug;
    if( my $aln = $in->next_aln ) {
	my %seen;
	for my $seq ( $aln->each_seq ) {
	    my $id = $seq->id;
	    if( $id =~ /(\S+)\|/) { 
		$id = $1;
	    }
		my $s = $seq->seq;
		$s =~ s/\./-/g;
	    $matrix{$id} .= $s;
	    $seen{$id}++;
	}
	for my $exp ( @expected ) {
	    if( ! $seen{$exp} ) {
		$matrix{$exp} .= '-' x $aln->length;
	    }
	}
    }
}

my $bigaln = Bio::SimpleAlign->new;

while( my ($id,$seq) = each %matrix ) {
 warn("$id and ",length($seq),"\n");
    $bigaln->add_seq(Bio::LocatableSeq->new(-id  => $id, #-end => length($seq),
					    -seq => $seq));
}
$bigaln->set_displayname_flat(1);

my $out = Bio::AlignIO->new(-format => $oformat, -idlength=>40,
			    -file   => ">$outfile");
$out->write_aln($bigaln);


