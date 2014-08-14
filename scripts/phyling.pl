#!env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Env qw(PHYLINGHOME);
use File::Spec;
use Getopt::Long;

use IO::String;
use Bio::SeqIO;
my @EXPECTED_APPS = qw(FASTQ_TO_FASTA HMMALIGN HMMSEARCH TRANSEQ
CDBFASTA CDBYANK PHRAP GENEWISEDB SREFORMAT TRIMAL FASTTREE MUSCLE);

my $app_conf;
if( $PHYLINGHOME) {
    $app_conf = File::Spec->catfile($PHYLINGHOME, "lib","apps.conf");
} elsif ($Bin) {
    $app_conf = File::Spec->catfile($Bin, "..","lib","apps.conf");
}

my $debug = 1;
my $qual_offset = 33;
my $evalue_cutoff = '1e-10';

my ($hmm2_models,$hmm3_models,$marker_hmm);
my $clade = 'fungi';
my $cpus = 1;
my $cleanup = 0;
my $force = 0; # force re-processing files even if there is cached intermediate
my $tmpdir;
my $prefix;
my $seqprefix;
my $rDNA_hmm;
GetOptions('ac|app|appconf:s' => \$app_conf,
	   'v|debug!'         => \$debug,
	   'q|qual:i'         => \$qual_offset,
	   't|temp|tmpdir:s'  => \$tmpdir,
	   'cleanup!'         => \$cleanup,
	   'p|prefix:s'       => \$prefix,
	   'sp|prefix:s'      => \$seqprefix,
	   
	   'force!'           => \$force,
	   'c|clade:s'        => \$clade,
	   'cpus:i'           => \$cpus,
	   'hmm:s'            => \$marker_hmm,
	   'hmm2:s'           => \$hmm2_models,
	   'hmm3:s'           => \$hmm3_models,
	   'rDNA!'            => \$rDNA_hmm,
    );

$hmm3_models ||= File::Spec->catdir($Bin,'..','DB','markers',$clade,"HMM3");

if( ! -d $hmm3_models ) {
    die("$hmm3_models is not a valid directory with HMM3 models\n");
}

$marker_hmm ||= File::Spec->catfile($hmm3_models, 'markers_3.hmmb');

if( ! -f $marker_hmm ) {
    die("need a valid marker protein HMM for the markers to extract from the reads\n");
}

$hmm2_models ||= File::Spec->catdir($Bin,'..','DB','markers',$clade,
				     'HMM2');
if( ! -d $hmm2_models ) {
    die("$hmm2_models is not a valid directory with HMM2 models\n");
}
$rDNA_hmm ||= File::Spec->catfile($Bin,'..','DB','markers',$clade,
				    'rDNA_3.hmmb');
if( ! -f $rDNA_hmm ) {
    $rDNA_hmm = undef;
}

my $error = 0;
debug("app conf is $app_conf\n");
if( ! $app_conf ) {
    die("Must provide app config file via $PHYLINGHOME or --app --appconf\n");
}
my $paths = &parse_config($app_conf);
for my $p ( @EXPECTED_APPS ) {
    if( ! exists $paths->{$p} || ! -x $paths->{$p} ) {
	warn("cannot find Application $p ",$paths->{$p},"\n");
	$error = 1;
    }
}
die("Error in App config file\n") if $error;

my $in_file = shift @ARGV;

my (undef,$dir,$fname) = File::Spec->splitpath($in_file);

if( ! $tmpdir ) {
    $tmpdir = $prefix.".PHYling";
} 

if(  ! -d $tmpdir ) {
    mkdir($tmpdir);
}

my $fasta_file;

if( $fname =~ /(\S+)\.(fastq|fq)/ ) {
    $prefix = $1 if ! defined $prefix;
    $fasta_file = File::Spec->catfile($tmpdir,"$1.fasta");
    if( ! -f $fasta_file || $force ) {
	my $cmd = sprintf("%s -Q %d -i %s -o %s",$paths->{'FASTQ_TO_FASTA'},
			  $qual_offset, $in_file,$fasta_file);
	debug("CMD: $cmd\n");
	`$cmd`;
    }

} elsif( $fname =~ /(\S+)\.(seq|fasta|fa|fas)/ ) {
    $prefix = $1 if ! defined $prefix;
    $fasta_file = $in_file;
} else {
    warn("unknown extension in file $fname ($in_file)\n");
    exit;
}

$seqprefix ||= $prefix;

# GET RID OF PROBLEMATIC CHARACTERS IN THE READ IDs
# some caching, if we already 
my $tmp_rename = File::Spec->catfile($tmpdir,$prefix.".fasta.fix");
if( $force || ! -f $tmp_rename ) {
    &fix_read_ids($fasta_file,$tmp_rename);
}
$fasta_file = $tmp_rename;

#INDEX THE READS FILE FOR LATER PROCESSING
&index_file($fasta_file);

# MAKE 6 FRAME TRANSLATION PROTEIN FILE
my $aafile = File::Spec->catfile($tmpdir,$prefix.".6frame.faa");
if( -x $paths->{TRANSEQ} ) {
    &make_6frame_transeq($fasta_file,$aafile);
} else {
    &make_6frame($fasta_file,$aafile);
}

# RUN PROTEIN HMMSEARCH OF MARKERS AGAINST TRANSLATION
my $match_pref = File::Spec->catfile($tmpdir,$prefix);

my ($marker_table,$marker_report) = &hmmsearch_markers($aafile,
						       $marker_hmm,
						       $match_pref);

my $reads_per_marker = &parse_hmmtable($marker_table);

my @trim_files;
for my $marker ( keys %$reads_per_marker ) {
    my @reads = keys %{$reads_per_marker->{$marker}};
    my $reads_file = File::Spec->catfile($tmpdir,$prefix.".$marker.r1.fasta");
    if( $force || ! -f $reads_file || 
	-M $reads_file > -M $marker_table) {
	&retrieve_reads($fasta_file,\@reads,$reads_file);
    }
    my $contigs = &assemble_reads_phrap($reads_file);
    my $pepfile = File::Spec->catfile($tmpdir,$prefix.".$marker.candidate.pep");
    debug("contigs file: $contigs\n");
    if( -f $contigs ) {
	&genewise_contigs(File::Spec->catfile($hmm2_models,$marker.".hmm"),
			  $contigs,$pepfile);
    }
    if( -f $pepfile && ! -z $pepfile ) {
	debug("pepfile: $pepfile\n");
	my $pepseq = &extract_peptide_genewise_output($pepfile);
	if( $pepseq ) {
	    warn("got a pepseq of length ", $pepseq->length, "\n");
	    my $pepresult = File::Spec->catfile($tmpdir,$prefix.".$marker.1.pep");
	    $pepseq->display_id($seqprefix);
	    Bio::SeqIO->new(-format => 'fasta',
			    -file   => ">$pepresult")->write_seq($pepseq);
	    my $marker_alnfile = File::Spec->catfile($tmpdir,
						     $prefix.".$marker.msa");
		
	    &hmmalign(File::Spec->catfile($hmm2_models,$marker.".hmm"),
		      $pepresult,$marker_alnfile);
	    my $marker_trimfile = &trim_aln($marker_alnfile);
	    if( $marker_trimfile ) {
		# build tree
		&build_NJtree($marker_trimfile);
		push @trim_files, $marker_trimfile;
	    }
	}
    } else {
	debug("no pepfile for $contigs\n");
    }
}

if( $rDNA_hmm ) {
    
}

END {
    if( $cleanup ) {	
	warn("rm -rf $tmpdir\n");
    }
}

sub trim_aln {
    my ($infile) = shift;
    my $outfile = $infile . '.trim';
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {
	my $cmd = sprintf("%s -in %s -out %s -automated1 -fasta
    }
}

sub hmmalign {
    my ($hmm_model, $infile, $outfile) = @_;
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {
	my $cmd = sprintf("%s --trim --amino %s %s > %s",
			  $paths->{HMMALIGN},
			  $hmm_model, $infile,$outfile.".stk");
	debug("CMD: $cmd\n");
	`$cmd`;
	$cmd = sprintf("%s clustal %s > %s",
		       $paths->{SREFORMAT},$outfile.".stk",
		       $outfile);
	debug("CMD: $cmd\n");
	`cmd`;
    }
    
}
sub extract_peptide_genewise_output {
    my ($infile) = shift;
    open(my $fh => $infile) || die $!;
    my $ready;
    my $seq;
    while(<$fh>) {
	next if /^\#/ || /^\s+$/;
	if( /^>Results/ ) {
	    while(<$fh>) {
		next if /^Making/;
		last if /^\/\//;
		$seq .= $_;
	    }
	}
    }
    my $pseq;
    if( $seq && $seq =~ /^>/ ) {
	my $iostring = IO::String->new($seq);
	my $in = Bio::SeqIO->new(-format => 'fasta',
				 -fh   => $iostring);
	
	while( my $s = $in->next_seq ) {
	    $pseq = $s if ! $pseq || $pseq->length < $s->length;
	}
    } else {
	debug("No seq for $infile\n");
    }
    $pseq;
}

sub genewise_contigs {
    my ($hmm_model,$infile,$outfile) = @_;
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {
	my $cmd = sprintf("%s -hmmer -pep -splice flat -init local -silent -quiet %s %s > %s",
			  $paths->{GENEWISEDB},
			  $hmm_model,$infile,$outfile);
	debug("CMD: $cmd\n");
	$rc = `$cmd`;
    }
    $rc;
}

sub assemble_reads_phrap {
    my $infile = shift;
    my $contigs = $infile.".contigs";
    if( $force || 
	! -f $contigs ||
	-M $infile > -M $contigs ) {	
	my $cmd = sprintf("%s %s",$paths->{PHRAP},$infile);
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    $contigs;
}

sub retrieve_reads {
    my ($infile,$reads_ar,$outfile) = @_;
    # index_file($infile); # would be redundant to call this many times
    # so should just assume it is indexed?
    my $cmd = sprintf("| %s %s.cidx > %s",
		      $paths->{CDBYANK},$infile,$outfile);
    debug("CMD: $cmd\n");
    open(my $to_cdbyank => $cmd) || die "Cannot open $cmd: $!\n";
    my $i = 0;
    for my $read ( @$reads_ar ) {
	print $to_cdbyank $read,"\n";
	$i++;
    }
    $i;
}
sub parse_hmmtable {
    my ($infile) = @_;
    open(my $fh => $infile) || die "cannot open $infile: $!";
    my $seen;
    while(<$fh>) {
	next if /^\#/;
	my @row = split(/\s+/,$_);
	my $t = $row[0];
	my $q = $row[3];
	# this may be unnecessary
	my $evalue = $row[6];
	next if $evalue > $evalue_cutoff;
	my $id = $t;
	$id =~ s/_[0-6]$//;
	$seen->{$q}->{$id}++;
    }
    $seen;
}


sub hmmsearch_markers {
    my ($seqdb, $markerdb,$hmmfilepref) = @_;
    my $table = $hmmfilepref.".hmmsearch.domtbl";
    my $rpt = $hmmfilepref.".hmmsearch.out";
    my $rc = 1;
    if( $force ||
	! -f $table ||
	-M $table > -M $seqdb ) {
	my $cmd = sprintf("%s -E %s --cpu %d --domtblout %s %s %s > %s ",
			  $paths->{HMMSEARCH}, 
			  $evalue_cutoff,$cpus,
			  $table,$markerdb,$seqdb,$rpt);
	debug("CMD: $cmd\n");
	$rc = `$cmd`;
    }
    ($table,$rpt);
}

sub make_6frame {
    my ($infile,$outfile) = @_;
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {

	my $in = Bio::SeqIO->new(-format => 'fasta', -file => $infile);
	my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile");
	while( my $s = $in->next_seq ) {
	    my $rseq = $s->revcom;
	    my $id = $s->display_id;
	    for my $frame ( 0,1,2) {
		my $tseq = $s->translate(-frame=> $frame,
					 -terminator => 'X');
		$tseq->display_id(sprintf("%s_%d",
					  $id,$frame+1));
		$out->write_seq($tseq);
		$tseq = $rseq->translate(-frame=> $frame,
					 -terminator => 'X');
		$tseq->display_id(sprintf("%s_%d",
					  $id,$frame+4));
		$out->write_seq($tseq);	    
	    }
	    $rc++;
	}
    }
    $rc;
}

sub make_6frame_transeq {
    my ($infile,$outfile) = @_;
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {
	my $cmd = sprintf("%s -trim -clean -frame 6 %s %s",
			  $paths->{TRANSEQ},$infile, $outfile);
	debug("CMD: $cmd\n");
	$rc=`$cmd`;
    }
    $rc;
}

sub index_file {
    my $infile = shift;
    my $rc = 1;
    if( $force ||
	! -f $infile.".cidx" || 
	-M $infile.".cidx" > -M $infile) {
	my $cmd = sprintf("%s %s",$paths->{CDBFASTA},$infile);
	debug("CMD: $cmd\n");
	$rc = `$cmd`;
    }
    $rc;
}

sub fix_read_ids {
    my ($infile,$outfile) = @_;
    my $rc = 1;
    open(my $fh => $infile) || die "Cannot open $infile: $!\n";
    open(my $ofh => ">$outfile") || die("Cannot open $outfile: $!\n");
    while(<$fh>) {
	if( /^>(\S+)/ ) {
	    my $id = $1;
	    $id =~ s/[-:\/]/_/g;
	    $_ = ">$id\n";
	    $rc++;
	} 
	print $ofh $_;
    }    
    close($fh);
    close($ofh);
    $rc;
}

sub parse_config {
    my $config = shift;
    my $apps = {};
    open(my $fh => $config) || die "cannot open $config: $!";
    while(<$fh>) {
	chomp;
	if(/([^=]+)=(\S+)/ ) {	
	    $apps->{$1} = $2;
	} else {
	    warn("cannot parse line $_\n");
	}
    }
    if( $debug ) {
	for my $app ( keys %$apps ) {
	    debug("app is $app with path = ",$apps->{$app},"\n");
	}
    }
    $apps;
}

sub debug {
    my $msg = shift;
    warn($msg) if $debug;
}
