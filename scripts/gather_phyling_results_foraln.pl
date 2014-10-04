#!env perl
use strict;
use warnings;
use Bio::SeqIO;
my $dir = shift || die $! ;
opendir(DIR, $dir) || die $!;

my %d;
for my $subdir (readdir(DIR) ) {
    next if ($subdir !~ /PHYling/);
    opendir(SUB,"$dir/$subdir") || die $!;
    for my $file (readdir(SUB) ) {
	if( $file =~ /(\S+)\.(\S+)\.1\.(pep|cdna)$/ ) {
	    my ($strain,$gene,$type) = ($1,$2,$3);
	    push @{$d{$gene}->{$type}}, [$strain, "$dir/$subdir/$file"];
	} elsif ($file =~ /(\S+)\.(\S+)\.candidate\.(cdna)$/ ) {
	    my ($strain,$gene,$type) = ($1,$2,$3);
	    push @{$d{$gene}->{$type}}, [$strain, "$dir/$subdir/$file"];
	}
    }
}
for my $g ( keys %d ) {
    mkdir("$dir/marker_aln");    
    my %out = ( 'pep' => Bio::SeqIO->new(-format => 'fasta', 
					 -file =>">$dir/marker_aln/$g.aa"),
		'cdna' => Bio::SeqIO->new(-format => 'fasta', 
					  -file =>">$dir/marker_aln/$g.cds"),
	);
    for my $type ( keys %{$d{$g}} ) {
	for my $n ( @{ $d{$g}->{$type} } ) {
	    my $seq = Bio::SeqIO->new(-format => 'fasta',
				      -file   => $n->[1]);
	    my $firstseq = $seq->next_seq;
	    if( $firstseq ) {
		my $id = $firstseq->display_id;
		my @r = split(/\./,$id);
		if( @r  >= 4 ) { 
		    splice(@r,-3);
		    $firstseq->display_id(join(".",@r));
		}
		$out{$type}->write_seq($firstseq);	
	    } else {
		warn("no seq in ",$n->[1],"\n");
	    }
	}
    }
}
