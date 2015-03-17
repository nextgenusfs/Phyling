#!env perl
use strict;
use warnings;
my $ctg = 'contigs.fa';
my $prefix = shift || die "cannot run without prefix.";
my $dir = shift || ".";
opendir(DIR,$dir) || die"cannot open $dir: $!";
my %stats;
for my $vdir (readdir(DIR) ) {
    next if $vdir =~ /^\./;
    next if ! -d "$dir/$vdir";
    next unless $vdir =~ /^\Q$prefix\E/;
    warn "processing $vdir\n";
    my $ctgsfile = "$dir/$vdir/$ctg";
    if( ! -f $ctgsfile ) {
	warn("no $ctgsfile!");
	next;
    }
    if ( -f "$dir/$vdir/Log" ) {
	open(my $fh => "$dir/$vdir/Log") || die "$vdir/Log: $!";
	while(<$fh>) {
	    if( /Final graph has (\d+) nodes and n50 of (\d+), max (\d+), total (\d+)/) {
		$stats{$vdir} = { contigs => $1,
				  n50     => $2,
				  max     => $3,
				  len   => $4 };
	    }
	}
    }
}
for my $vdir ( sort { $stats{$b}->{max} <=> $stats{$a}->{max} } 
	       keys %stats ) {
    print join("\t", $vdir, map { $stats{$vdir}->{$_} } 
	       qw(contigs n50 max len)),"\n";
}
