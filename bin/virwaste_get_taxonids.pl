#!/usr/bin/perl
#   
#   The input is a list of Genebank ids.
#             and the taxonomizator output. 
#   This script will generate 2 files:
#   - list of taxonids corresponding to the given ids.
#   - List of taxonids from ALL species in the same genus as given ids species (+ the genus taxonids)
# USAGE:
#   virwaste_get_taxonids.pl  . . . . . . > . . . 

use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils qw(uniq);

my $taxonfl = $ARGV[0];
my $idslist  = $ARGV[1];
my $specoutfl  = $ARGV[2];
my $genusoutfl = $ARGV[3];

print STDERR "### READING GENEBANK IDS FILE\n";
open(IDF, $idslist);
my @gbids;
my $id;
while (<IDF>) {
    chomp;
    my @a= split '\.', $_;
    $id=$a[0];
    push(@gbids, $id);
};
close(IDF);

print STDERR "### READING TAXONOMY FILE\n";
open(TAX, $taxonfl);
my %bygenustx; #This hash will be filled with genus taxids as keys and species taxids as values (list).
my @spectxids;
my @genustxids;
while (<TAX>) {
    chomp;
    my @ln=split /\t/o, $_;
    my @b=split '\.', $ln[0];
    $id=$b[0];
    my $stx="NA";
    my $gtx="NA";
    my @g;
    ## Store in  variables genus and specie taxon ids.
    for (my $i=0; $i < scalar(@ln); $i++) {
            if ($ln[$i] =~ /^species/o) {
               @g = split /:/o, $ln[$i];
               $stx = $g[2];
            }elsif($ln[$i] =~ /^genus/o){
                @g = split /:/o, $ln[$i];
                $gtx = $g[2];
           };
    };
     #print STDERR $id;
     #print STDERR " ### ";
     #print STDERR $stx;
     #print STDERR " ### ";
     #print STDERR $gtx;
     #print STDERR " ###\n ";
    
    ## Append to hash.
    if (exists $bygenustx{$gtx}) {
        push @{$bygenustx{$gtx}}, $stx;
    } else {
        @bygenustx{$gtx}=[ $stx ];
    };

    ## If specie is in the genebank ids list store it in a new list.
    for (my $i=0; $i < @gbids; $i++) {
            if ($id eq $gbids[$i]) {
               push @spectxids, $stx;
               push @genustxids, $gtx;
               last;
            };
    };
};
close(TAX);
@genustxids = uniq @genustxids;
@spectxids = uniq @spectxids;

#Create a list with all the species included in the genus of interest.
my @biglist;
foreach my $g (@genustxids) {
    
    for (my $i=0; $i<length($bygenustx{$g}); $i++){
        push @biglist, $bygenustx{$g}[$i];
    };
};
@biglist = uniq @biglist;

## print files

print STDERR "### WRITTING INFORMATION\n";
open ( ROUT, ">", $specoutfl );
foreach (@spectxids) { print ROUT join("", $_, "\n"); };
close(ROUT);

print STDERR "### WRITTING INFORMATION\n";
open ( GOUT, ">", $genusoutfl );
print GOUT "#GENUS_IDS\n";
foreach (@genustxids) { print GOUT join("", $_, "\n"); };
print GOUT "#SPECIES_IDS (species included in the previous genus)\n";
foreach (@biglist) { print GOUT join("", $_, "\n"); };
close(GOUT);
