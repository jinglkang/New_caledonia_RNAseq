#!/usr/bin/perl
use strict;
use warnings;

#!/usr/bin/perl
use strict;
use warnings;

# ~/software/database/swiss_pro_info.txt
my $uni="/home/jlkang/software/database/uniprot_sprot.fasta";
my %seq; my ($ID);
open UNI, $uni or die "can not open $uni\n";
while (<UNI>) {
        chomp;
        if (/>/) {
                s/\>//;
                my @a=split;
                $ID=$a[0];
        } else {
                $seq{$ID}.=$_;
        }
}

my @blasts=<final_blast_orth_group/*\.blastp_result>;
foreach my $blast (@blasts) {
        (my $orth)=$blast=~/final_blast_orth_group\/(.*)\.blastp_result/;
        open ORTH, $blast or die "can not open $blast\n";
        while (<ORTH>) {
                chomp;
                my @a=split /\t/;
                print ">$orth\n$seq{$a[0]}\n";
                last;
        }
}
