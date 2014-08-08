#!env perl
use strict;
use warnings;
use Getopt::Long;

=head1 NAME 

fasta_name_checker - fix names which a

=head1 USAGE

 -t / --title  specify the read title base (default: Unknown)
 -p / --prefix specify the read prefix (default: UNK)

=head1 AUTHOR Jason Stajich

Jason Stajich - jason[at]bioperl.org

=cut

my $count = 1;
my $title = 'Unknown';
my $prefix = 'UNK';

GetOptions("t|title:s"  => \$title,
	   "p|prefix:s" => \$prefix);
while(<>) {
    if( /^>/ ) {
          # $count keeps track of number of unique IDs for replacement
	if( s/no name/$title\_$count/i ) { $count++ }
	if( ! /\|/ ) {
	    s/>/>$prefix|/;
	}
	s/[:\/]/_/g;
    }
    print;
}
