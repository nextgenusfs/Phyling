#!env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my $hmmbuild3 = '/opt/hmmer/3.1b1/bin/hmmbuild';
my $hmmbuild2 = '/opt/hmmer/2.3.2/bin/hmmbuild2';
my $hmmcalibrate = '/opt/hmmer/2.3.2/bin/hmmcalibrate2';
my $target_hmm_folder = '/shared/stajichlab/projects/Phylogenomics/Phyling/DB/markers';
my $grouping;
my $indir;
my $ext = 'stk';
GetOptions(
    'target:s' => \$target_hmm_folder,
    'g|grouping:s' => \$grouping,
    'i|indir:s'    => \$indir,
    );



die ("no grouping specified\n") unless $grouping;

my $locale2 = File::Spec->catdir($target_hmm_folder,
				 $grouping, 'HMM2');
my $locale3 = File::Spec->catdir($target_hmm_folder,
				 $grouping, 'HMM3');
if( ! -d $locale2 ) {
    mkdir($locale2);
    mkdir($locale3);
}

opendir(IN,$indir)|| die $!;
for my $file ( readdir(IN) ) {
    next unless $file =~ /(\S+)\.\Q$ext\E$/;
    my $stem = $1;
     warn($stem,"\n");
#    next unless $stem eq 'OGFLEX_1065';
   #`sreformat stockholm $indir/$file > $indir/$stem.stk`;
    `$hmmbuild3 --seed 11 --amino $locale3/$stem.hmm $indir/$file`;
    `$hmmbuild2 -g $locale2/$stem.hmm $indir/$file`;
    `$hmmcalibrate --seed 11 $locale2/$stem.hmm`;

}
