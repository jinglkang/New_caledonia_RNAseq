#!/usr/bin/perl
use strict;
use warnings;

# ~/software/database/swiss_pro_info.txt
my $ann="/home/jlkang/software/database/swiss_pro_info.txt";
my %anno;
open ANN, $ann or die "can not open $ann\n";
while (<ANN>) {
        chomp;
        my @a=split /\t/;
        (my $des)=$a[1]=~/(.*?)\s+\[/;
        $anno{$a[0]}=$des;
}

my @blasts=<final_blast_orth_group/*\.blastp_result>;
foreach my $blast (@blasts) {
        (my $orth)=$blast=~/final_blast_orth_group\/(.*)\.blastp_result/;
        open ORTH, $blast or die "can not open $blast\n";
        while (<ORTH>) {
                chomp;
                my @a=split /\t/;
                my ($uni, $gene)=$a[0]=~/sp\|(.*?)\|(.*?)\_/;
                print "$orth\t$a[0]\t$gene\t$anno{$uni}\n";
        }
}
