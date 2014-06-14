#!/usr/bin/perl
use strict;
use warnings;
use Bio::TreeIO;

my ($treefile,$taxon_name) = @ARGV;
my $in = Bio::TreeIO->new(-format => 'newick', -file => $treefile);
my $tree = $in->next_tree;
if( ! $tree ) {
 die "cannot parse treefile $treefile and find a valid tree";
}
my $node = $tree->find_node(-id => $taxon_name);
if( $node ) { 
 my $parent = $node->ancestor;
 my @neighbors;
 for my $c ( grep { $_->is_Leaf } $parent->get_all_Descendents ) {
   push @neighbors, $c->id if $c->id ne $node->id;
 }
 printf "Neighbor(s) of %s are the node(s): %s\n", $node->id, join(",", @neighbors);

} else {
 warn("cannot find node $taxon_name in the tree, make sure you spelled it correctly, this is an exact match search");
}