#!env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Copy qw(move);
use Env qw(PHYLINGHOME);
use File::Spec;
use Getopt::Long;
use File::Temp qw(tempfile);

use IO::String;
use Bio::SeqIO;
use Bio::AlignIO;

my $buffer_end_start = 3; # what is the sloppy overhang allowed when bringing in more seqs

my $exonerate_options = '-m p2g --bestn 1 --joinfilter 1 --verbose 0 --ryo ">%ti (%tab - %tae) score=%s rank=%r\n%tcs\n" --showcigar no --showvulgar no --showalignment no --exhaustive';

my %uncompress = ('bz2' => 'bzcat',
		  'gz'  => 'zcat');

my @EXPECTED_APPS = qw(FASTQ_TO_FASTA HMMALIGN HMMSEARCH TRANSEQ 
                       FASTA TFASTY HMMEMIT CAP3
                       CDBFASTA CDBYANK PHRAP SREFORMAT 
                       TRIMAL FASTTREE MUSCLE EXONERATE);

$ENV{WISECONFIGDIR} = '/opt/wise/2.4.0/wisecfg';
my $app_conf;
if( $PHYLINGHOME) {
    $app_conf = File::Spec->catfile($PHYLINGHOME, "lib","apps.conf");
} elsif ($Bin) {
    $app_conf = File::Spec->catfile($Bin, "..","lib","apps.conf");
}

my $debug = 0;
my $qual_offset = 33;
my $hmmer_cutoff = '1e-10';
my $contig_match_cutoff = '1e-8';
my $contig_transsearch_cutoff = '1e-3';
my $Max_rounds = 10; # max iterations
my $scaffold_separator = 'N'x15; # 5 amino acid break in the scaffolded contigs, codon
my $njtree_options = "-boot 1000 -wag -seed 121 -bionj";
my ($hmm2_models,$hmm3_models,$marker_hmm,$marker_fasta_dir);
my $clade = 'AFTOL70';
my $CPUs = 1;
my $cleanup = 0;
my $force = 0; # force re-processing files even if there is cached intermediate
my $tmpdir;
my $prefix;
my $seqprefix;
my $rDNA_hmm;
my $do_MSA = 0;
my $consensus_folder = 'consensus';
my $port = 8001+int rand(100);
GetOptions('ac|app|appconf:s' => \$app_conf,
	   'v|debug!'         => \$debug,
	   'force!'           => \$force,
	   'maxrounds:i'      => \$Max_rounds,
	   'q|qual:i'         => \$qual_offset,
	   't|temp|tmpdir:s'  => \$tmpdir,
	   'cleanup!'         => \$cleanup,
	   'p|prefix:s'       => \$prefix,
	   'sp|prefix:s'      => \$seqprefix,	   

	   'c|clade:s'        => \$clade,
	   'cpus|cpu:i'           => \$CPUs,
	   'hmm:s'            => \$marker_hmm,
	   'hmm2:s'           => \$hmm2_models,
	   'hmm3:s'           => \$hmm3_models,
	   'md|markerdir:s'   => \$marker_fasta_dir,
	   'rDNA!'            => \$rDNA_hmm,
	   'cons|consensus:s' => \$consensus_folder,
	   
	   'msa!'             => \$do_MSA,
	   'port:s'           => \$port,
	   
    );

mkdir($consensus_folder) unless -d $consensus_folder;

$hmm3_models ||= File::Spec->catdir($Bin,'..','DB','markers',$clade,"HMM3");

if( ! -d $hmm3_models ) {
    die("$hmm3_models is not a valid directory with HMM3 models\n");
}

$marker_hmm ||= File::Spec->catdir($Bin,'..','DB','markers',$clade,"markers_3.hmmb");

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

$marker_fasta_dir ||= File::Spec->catdir($Bin,'..','DB','markers',$clade,
				    'marker_files');
#if( ! -d $marker_fasta_dir ) {
#    die("$marker_fasta_dir is not a valid directory with fasta sequences per model\n");
#}

my $error = 0;
debug("app conf is $app_conf\n");
if( ! $app_conf ) {
    die("Must provide app config file via $PHYLINGHOME or --app --appconf\n");
}
my $paths = &parse_config($app_conf);
for my $p ( @EXPECTED_APPS ) {
    if( ! exists $paths->{$p} || ! -x $paths->{$p} ) {
	debug("cannot find Application $p ",$paths->{$p},"\n");
	$error = 1;
    }
}
die("Error in App config file\n") if $error;

my $in_file = shift @ARGV;

my (undef,$dir,$fname) = File::Spec->splitpath($in_file);

my $base = $prefix;
if( ! $prefix && $fname =~ /(\S+)\.(fq|fq|fast\w+|seq)/) {
    $base = $1;
} else {
    $base = $$;
}
if( ! $tmpdir ) {
    $tmpdir = $base.".PHYling";
} 
if(  ! -d $tmpdir ) {
    mkdir($tmpdir);
}


my $fasta_file;

my $data_type;
if( $fname =~ /(\S+)\.(fastq|fq)/ ) {
    $prefix = $1 if ! defined $prefix;
    $fasta_file = File::Spec->catfile($tmpdir,"$1.fasta");
    if( ! -f $fasta_file || $force ) {
	my $cmd;
	if( $fname =~ /\.(bz2|gz)$/ ) {
	    $cmd = sprintf("%s %s | %s -Q %d -o %s",$uncompress{$1},$in_file,
			   $paths->{'FASTQ_TO_FASTA'},$qual_offset,$fasta_file);
	} else {
	    $cmd = sprintf("%s -Q %d -i %s -o %s",$paths->{'FASTQ_TO_FASTA'},
			   $qual_offset, $in_file,$fasta_file);
	}
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    $data_type = 'FASTQ';
} elsif( $fname =~ /(\S+)\.(seq|fasta|fa|fas)/ ) {
    $prefix = $1 if ! defined $prefix;
    $fasta_file = $in_file;
    $data_type = 'FASTA';
} else {
    debug("unknown extension in file $fname ($in_file)\n");
    exit;
}
$seqprefix ||= $prefix;

# GET RID OF PROBLEMATIC CHARACTERS IN THE READ IDs
# some caching, if we already 
my $tmp_rename = File::Spec->catfile($tmpdir,$prefix.".fasta.fix");
if( $force || ! -f $tmp_rename ) {
    &fix_read_ids($fasta_file,$tmp_rename);
} else {
    debug("Read file already renamed\n");
}
$fasta_file = $tmp_rename;

#INDEX THE READS FILE FOR LATER PROCESSING
&index_file($fasta_file);

# make 
my $bitfile = &make_2bit_file($fasta_file);

my $blat_ready;
unless( my $pid = fork() ) {
    $blat_ready = &start_gfServer($port,$bitfile);
} else {
    sleep(60);
    warn("here in parent prcess\n");
    warn("Done with fork test\n");

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
	my $marker_cons = File::Spec->catfile($consensus_folder,"$marker.cons");    
	if( ! -f $marker_cons ) {
	    &make_consensus_HMM(File::Spec->catfile($hmm3_models,$marker.".hmm"),
				$marker_cons);
	} else {
	    debug("consensus HMM for $marker_cons already created\n");
	}

	#my $marker_seqs = &read_marker_refproteins($marker_fasta_dir,$marker);
	my $cdnafile = File::Spec->catfile($tmpdir,$prefix.".$marker.candidate.cdna");
	my $pepfile = File::Spec->catfile($tmpdir,$prefix.".$marker.candidate.pep");
	my $scaffoldfile = File::Spec->catfile($tmpdir,$prefix.".$marker.ord_scaf.fa");
	debug("cdnafile is $cdnafile\n");
	if( -f $cdnafile ) {
	    debug("cdnafile exists\n");
	}
	next if ( ! $force && -f $cdnafile);
	if( ! -f $scaffoldfile || $force ) {
	    my @reads = keys %{$reads_per_marker->{$marker}};
	    my $reads_file = File::Spec->catfile($tmpdir,$prefix.".$marker.r1.fasta");
	    if( $force || ! -f $reads_file || 
		-M $reads_file > -M $marker_table) {
		&retrieve_reads($fasta_file,\@reads,$reads_file);
	    }

#	my $contigsfile = &assemble_reads_phrap($reads_file);
	    my $contigsfile = &assemble_reads_cap3($reads_file);
	    my $contig_count = &seqcount($contigsfile);
	    if( $contig_count == 0 ) {
		warn("No assembled contigs or singlets available\n");
		next;
	    }
	    debug("seqcount for $contigsfile is $contig_count\n");
	    my $change = 1;
	    my $rounds = 0;
	    while( $contig_count > 1 && 
		   $change > 0 && $rounds < $Max_rounds) {
		my $added = &search_and_add($fasta_file,$contigsfile,$reads_file);
		warn("added $added reads to $reads_file\n");
		last if $added == 0;

		#$contigsfile = &assemble_reads_phrap($reads_file);
		$contigsfile = &assemble_reads_cap3($reads_file);
		my $newcount = &seqcount($contigsfile);
		warn("contig count was $contig_count newcount is $newcount\n");
		$change = ($contig_count - $newcount);
		$contig_count = $newcount;
		$rounds++;
	    }
	    warn("performed $rounds of additional read gathering\n");

	    my $updated_contigs = &stitch_order_contigs($marker_cons,$contigsfile);
	    # merge the contigs, in their new order, into one scaffold with some Ns between
	    my $scaff_seq = join($scaffold_separator, (map { defined $_  ? $_->seq : '' } @$updated_contigs));
	    my $scaffold = Bio::Seq->new(-id => "$prefix.$marker.scaffold",
					 -seq => $scaff_seq);
	    Bio::SeqIO->new(-format => 'fasta', -file =>">$scaffoldfile")->write_seq($scaffold);

	    debug("Scaffold file: $scaffoldfile\n");
	}

	if( -f $scaffoldfile ) {
	    &exonerate_best_model($marker_cons,$scaffoldfile,$cdnafile);
	    &translate_cdna($cdnafile,$pepfile);
	} else {
	    warn("no scaffold to process for $marker\n");
	}
	last if $debug;
    }

    if( $rDNA_hmm ) {

    }
}

END {
 &stop_gfServer($port);

}

sub start_gfServer {
    my ($Port,$bitfile) = @_;
    warn("Port is $Port bitfile is $bitfile\n");
    my $cmd = sprintf("%s start localhost %d %s -canStop",
		      $paths->{GFSERVER},$Port,$bitfile);
    debug("CMD: $cmd\n");
    `$cmd`;
    return 1;
}
sub stop_gfServer {
    my ($Port) = @_;
    system($paths->{GFSERVER},'stop','localhost',$Port);    
}

sub make_2bit_file {
    my ($fasta) = shift;
    my ($base) = $fasta;
    $base .= ".2bit";
    if( $force ||
	! -f $base ||
	-M $base > -M $fasta ) {
	my $cmd = sprintf("%s %s %s",$paths->{FATOTWOBIT},$fasta,$base);
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    $base;
}

sub trim_aln {
    my ($infile) = shift;
    my $outfile = $infile . '.trim';
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {
	my $in = Bio::AlignIO->new(-format => 'clustalw',
				   -file   => $infile);

	my $out = Bio::AlignIO->new(-format => 'fasta',
				    -file   => ">$infile.2");
	if( my $aln = $in->next_aln ) {
	    $aln->map_chars('\.','-');
	    $aln->set_displayname_flat(1);
	    $out->write_aln($aln);	    
	}
	move("$infile.2",$infile);
	my $cmd = sprintf("%s -in %s -out %s -automated1 -fasta",
			  $paths->{TRIMAL},$infile,$outfile);
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    $outfile;
}

sub hmmalign {
    my ($hmm_model, $inseqs, $outfile) = @_;
    my $rc = 1;
    debug("outfile is $outfile\n");
    if( $force || ! -f $outfile ) {
	my ($ifh,$infile) = tempfile('tmpXXXX',UNLINK=>1);
	my $out = Bio::SeqIO->new(-format => 'fasta', -fh => $ifh);
	$out->write_seq(@$inseqs);	
	close($ifh);
	my $cmd = sprintf("%s --trim --amino %s %s > %s",
			  $paths->{HMMALIGN},
			  $hmm_model, $infile,
			  $outfile.".stk");
	debug("CMD: $cmd\n");
	`$cmd`;
	$cmd = sprintf("%s clustal %s > %s",
		       $paths->{SREFORMAT},$outfile.".stk",
		       $outfile);
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    
}

=head2 exonerate_best_model

 Title   : exonerate_best_model
 Usage   :
 Function: Finds best protein2genome alignment with exonerate
 Example :
 Returns : Status of run. Creates a CDS out file
 Args    :


=cut


sub exonerate_best_model {
    my ($inpepfile,$contigfile,$outfile) = @_;
    my $rc = 0;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $inpepfile ) {

	my $cmd = sprintf("%s %s %s %s > %s",
			  $paths->{EXONERATE},
			  $inpepfile,$contigfile,
			  $exonerate_options,
			  $outfile);
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
	-M $infile < -M $contigs ) {	
	my $cmd = sprintf("%s %s",$paths->{PHRAP},$infile);
	debug("CMD: $cmd\n");
	`$cmd`;
	`cat $infile.singlets >> $contigs`;
    }
    if( ! -f "$contigs.renum" ||
	-M $contigs < -M "$contigs.renum" ) {	

	my $renumber = Bio::SeqIO->new(-format => 'fasta',
				       -file   => $contigs);

	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">$contigs.renum");
	my $i = 1;
	while( my $s = $renumber->next_seq ) {
	    $s->display_id(sprintf("ctg%d",$i++));
	    $out->write_seq($s);
	}
	$out->close();
	$out = undef;
    }
    "$contigs.renum";
}

sub assemble_reads_cap3 {
    my $infile = shift;
    my $contigs = $infile.".cap.contigs";
    if( $force ||
        ! -f $contigs ||
        -M $infile < -M $contigs ) {
        my $cmd = sprintf("%s %s",$paths->{CAP3},$infile);
        debug("CMD: $cmd\n");
        `$cmd`;
        #`cat $infile.singlets >> $contigs`;
    }
    $contigs;
}

sub retrieve_reads {
    my ($infile,$reads_ar,$outfile) = @_;
    if( ! -f "$infile.cidx" ) {
	&index_file($infile);
    }
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

sub get_read {
    my ($infile,$read_name) = @_;
    if( ! -f "$infile.cidx" ) {
	&index_file($infile);
    }
    my $cmd = sprintf("%s %s.cidx -a %s |",
		      $paths->{CDBYANK},$infile,$read_name);
    debug("CMD: $cmd\n");
    open(my $readseq => $cmd) || die "Cannot open $cmd: $!\n";
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $readseq);
    my @seqs;
    while(my $s = $seqio->next_seq ) {
	push @seqs, $s;
    }
    @seqs;
}

sub read_marker_refproteins {
    my ($dir,$marker_name) = @_;
    my $marker_file = File::Spec->catfile($dir,$marker_name);
    for my $ext ( qw(fa fas fasta pep aa seq) ) {
	if( -f $marker_file .".$ext" ) {
	    $marker_file .= ".$ext";
	    last;
	}
    }
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
				-file   => $marker_file);
    my $seqs = [];
    while( my $seq =$seqio->next_seq ) {
	push @{$seqs}, $seq;
    }
    $seqs;
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
	next if $evalue > $hmmer_cutoff;
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
			  $hmmer_cutoff,$CPUs,
			  $table,$markerdb,$seqdb,$rpt);
	debug("CMD: $cmd\n");
	$rc = `$cmd`;
    }
    ($table,$rpt);
}

sub translate_cdna {
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
	    my $tseq = $s->translate(-frame=> 0,
				     -terminator => 'X');
	    $out->write_seq($tseq);
	}

    }
    $rc;
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
    } else {
	debug("6frame translation already created");
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
    } else {
	debug("6frame translation EMBOSS already run\n");
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
    } else {
	debug("index file already created\n");
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
	    $id =~ s/[-:\/#|]/_/g;
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
	while( my ($app,$path) = each %$apps ) {
	    debug("app is $app with path = $path\n");
	}
    }
    $apps;
}

sub seqcount {
    my $file = shift;
    open(my $fh => "grep -c '^>' $file |") || die $!;
    my $n = <$fh>;
    $n =~ s/\s+//g;
    $n;
}


sub seq_lengths {
    my $file = shift;
    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $file);
    my $res = {};
    while( my $s = $in->next_seq ) {
	debug($s->display_id. " ". $s->length,"\n");
	$res->{$s->display_id} = $s->length;
    }
    $res;
}

=head2 stitch_order_contigs

 Title   : stitch_order_contigs
 Usage   : &stitch_order_contigs($markerpep,$contigs);
 Function: Reorder and scaffold contigs based on a protein query sequence from 
 Returns : Update contigs file
 Args    :


=cut

sub stitch_order_contigs {
    my ($marker_cons,$contigfile) = @_;
    my $cmd = sprintf("%s -T %d -m 8c -E %s %s %s",
		      $paths->{TFASTY}, $CPUs,
		      $contig_transsearch_cutoff,
		      $marker_cons,
		      $contigfile);
    debug("running $cmd\n");
    open(my $run => "$cmd |") || die "cannot run: $cmd\n";
    my @results;
    while(<$run>) {
	next if /^\#/;
	chomp;
	my ($q,$h,$pid,$match,$mismatch,$gap,$qstart,$qend,
	    $tstart,$tend,$evalue,$bits) = split(/\t/,$_);
	my ($tstrand) = (1,1);	
	if( $tstart > $tend ) { 
	    ($tend,$tstart) = ($tstart,$tend);
	    $tstrand = -1;
	}
	debug("result is $qstart,$qend,$h,$tstart,$tend,$tstrand,$evalue,$bits,$pid \n");
	push @results, [$qstart,$qend,$h,$tstart,$tend,$tstrand,$evalue,$bits,$pid];
    }
    my %contigs;
    my $read_contigs = Bio::SeqIO->new(-format => 'fasta', -file => $contigfile);
    while(my $s = $read_contigs->next_seq ) {
	debug("seq id is ", $s->display_id,"\n");
        $contigs{$s->display_id} = $s;
    }
    my $new_order;
    # sort by query protein alignment order
    my %seen;
    for my $res ( sort { $a->[0] <=> $b->[0] } @results ) {
	debug (join("\t", @$res),"\n");
	next if $seen{$res->[2]}++;
	next if ! defined $contigs{$res->[2]};
	if( $res->[5] < 0 ) {
	    push @$new_order, $contigs{$res->[2]}->revcom;
	} else {
	    push @$new_order, $contigs{$res->[2]};
	}
    }

    $new_order; #bizzare love triangle
}

sub search_and_add {
    my ($searchdb,$queryfile,$outputreads) = @_;
    
    my $qlens = &seq_lengths($queryfile);
    my $rlens = &seq_lengths($outputreads);
    my (undef,$searchdir,$fname) = File::Spec->splitpath($searchdb);
    $searchdir = '.';
# replace this with BLAT and near-identity?
#    my $cmd = sprintf("%s -T %d -E %s -m 8c %s %s |",
#		      $paths->{FASTA},$CPUs,$contig_match_cutoff,
#		      $queryfile, $searchdb);
    my $cmd = sprintf("%s %s %s %s %s stdout -out=blast8 |",
		      $paths->{GFCLIENT},'localhost',$port,$searchdir,$queryfile);
    
    debug("CMD: $cmd\n");
    open(my $fasta_res => $cmd) || die "cannot run $cmd\n";
    my @results;
    my @readnames;
    while(<$fasta_res>) {
	next if /^\#/;
	debug($_);
	chomp;
	my ($q,$h,$pid,$match,$mismatch,$gap,$qstart,$qend,
	    $tstart,$tend,$evalue) = split(/\t/,$_);
	next if( exists $rlens->{$h});
	my ($qstrand,$tstrand) = (1,1);
	if( $qstart > $qend ) { 
	    ($qend,$qstart) = ($qstart,$qend);
	    $qstrand = -1;
	}
	if( $tstart > $tend ) { 
	    ($tend,$tstart) = ($tstart,$tend);
	    $tstrand = -1;
	}
	if($tstart > $buffer_end_start ) { # if target alignment start is not 1 or some 
	                                   # number close to 1 (buffer_end_start)
	    push @readnames, $h;
	} elsif( abs($qlens->{$q}-$qend) >= $buffer_end_start ) {
	    # -----|        Query
	    # ------------| Hit
	    debug("$q overhanging $qend vs length: ".$qlens->{$q}."\n");
	    my ($read_seq) = &get_read($searchdb,$h);
	    my $read_len = $read_seq->length;
	    if( $read_len > $tend) { # if the end of read align (tend) after the end of this aln
		push @readnames, $h;
	    }
	}
    }
    debug("readnames are @readnames\n");
    if( @readnames ) {
	&retrieve_reads($searchdb,\@readnames,"$outputreads.add");
	my $in = Bio::SeqIO->new(-format => 'fasta',
				 -file   => "$outputreads.add");
	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">>$outputreads");
	while( my $s = $in->next_seq ) {
	    $out->write_seq($s);
	}
    }
    return scalar @readnames;
}

sub make_consensus_HMM {
    my ($hmmfile,$outfile) = @_;
    my $cmd = sprintf("%s -c %s > $outfile",
		      $paths->{HMMEMIT},$hmmfile);
    `$cmd`;
}

sub debug {
    my @msg = @_;
    warn(join(" ",@msg)) if $debug;
}


END {
    if( $cleanup ) {	
	warn("rm -rf $tmpdir\n");
    }
}
