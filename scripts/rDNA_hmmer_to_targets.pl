#!env perl
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Copy qw(move);
use Env qw(PHYLINGHOME);
my $debug =0;

my $cutoff = 10;
my @EXPECTED_APPS = qw(FASTQ_TO_FASTA HMMALIGN HMMSEARCH TRANSEQ
CDBFASTA CDBYANK PHRAP GENEWISEDB SREFORMAT TRIMAL FASTTREE MUSCLE);

$ENV{WISECONFIGDIR} = '/opt/wise/2.4.0/wisecfg';
my $app_conf;
if( $PHYLINGHOME) {
    $app_conf = File::Spec->catfile($PHYLINGHOME, "lib","apps.conf");
} elsif ($Bin) {
    $app_conf = File::Spec->catfile($Bin, "..","lib","apps.conf");
}

my ($fasta_file,$prefix,$force);

GetOptions(
    'v|debug!'   => \$debug,
    'p|prefix:s' => \$prefix,
    'f|fasta:s'  => \$fasta_file,
    );

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

my $report = shift || die "need a report file";

if( ! $prefix ) {
    ($prefix) = split(/\./,$report);
}

open(my $fh => $report) || die "cannot open $report: $!";

my %seen;
while(<$fh>) {
    next if /^\#/;
    my @row = split(/\s+/,$_);
    my $t = $row[0];
    my $q = $row[3];
    my $evalue = $row[11];
    next if $evalue > $cutoff;
    if( exists $seen{$q} ) {
	next;
    }
    $seen{$t} = $q;
}

my @reads = keys %seen;
my $reads_file = File::Spec->catfile($prefix.".rDNA.fasta");
if( $force || ! -f $reads_file ||
    -M $reads_file > -M $report ) {
    &retrieve_reads($fasta_file,
		    \@reads,$reads_file);
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

sub debug {
    my $msg = shift;
    warn($msg) if $debug;
}
