#!env perl
use strict;

my $dir = shift || die $! ;
opendir(DIR, $dir) || die $!;

my %d;
for my $subdir (readdir(DIR) ) {
    next if ($subdir !~ /PHYling/);
    opendir(SUB,"$dir/$subdir") || die $!;
    for my $file (readdir(SUB) ) {
	next unless $file =~ /(\S+)\.(\S+)\.1\.pep$/;
	my ($strain,$gene) = ($1,$2);
	push @{$d{$gene}}, [$strain, "$dir/$subdir/$file"];
    }
}
for my $g ( keys %d ) {
    mkdir("$dir/marker_aln");
    open(my $fh => ">$dir/marker_aln/$g.aa") || die $!;
#    print $g,"\n";
    for my $n ( @{$d{$g}} ) {
#	print join("\t", '',@$n), "\n";
	open(my $ifh => $n->[1]) || die $n->[1], " : $!";
	while(<$ifh>) {
	    print $fh $_;
	}
    }
}
