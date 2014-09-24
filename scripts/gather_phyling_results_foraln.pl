#!env perl
use strict;

my $dir = shift || die ;
opendir(DIR, $dir);

for my $subdir (readdir(DIR) ) {
    next if ($subdir !~ /PHYling/);
    opendir(SUB,"$dir/$subdir")
}
