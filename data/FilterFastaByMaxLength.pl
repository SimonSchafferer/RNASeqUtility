#!/usr/bin/perl
use strict;
use warnings;

my $maxLength = shift or die "Please provide 'maxLength' paramter\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen <= $maxLength);
    }
    local $/="\n";
}
