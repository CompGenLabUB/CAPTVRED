#!/usr/bin/perl
#
# USAGE:
# ./taxonomizator.pl merged.dmp names.dmp nodes.dmp columids.genus.tmp.tbl
use strict;
use warnings;
use Data::Dumper;

##
my %updatetx = ();
my %ids      = ();
my %names    = ();
my %refnames = ();
my %tree     = ();
my @O        = ();
my %genus    = ();


#
# READING MERGED.DMP
# 1.  old_taxid      -- depracated taxonomy identifier
# 2.  new_taxid    -- updated taxonomy identifier
#
print STDERR "### READING MERGED IDS\n";
open(MG, shift @ARGV);
my ($n, $c) = (0, 0);
while (<MG>) {
    our (@ln, $id);
    chomp;
    $_ =~ s/\t\|//og;
    @ln = split /\t/o, $_;
    $updatetx{$ln[0]}=$ln[1];
    print STDERR "old: $ln[0]  --  new: $ln[1]\n";
    &counter(++$n,\$c);
};
close(MG);
print STDERR "\n# Total MERGED: $n\n";


#
# READING NAMES.DMP
# 1.  tax_id      -- the id of node associated with this name
# 2.  name_txt    -- name itself
# 3.  unique name -- the unique variant of this name if name not unique
# 4.  name class  -- scientific name, synonym, common name, ...
#
print STDERR "### READING NAMES IDS\n";
open(NI, shift @ARGV);
($n, $c) = (0, 0);
while (<NI>) {
    our (@F, $id);
    chomp;
    $_ =~ s/\t\|//og;
    @F = split /\t/o, $_;
    $id = $F[2] eq q{} ? $F[1] : $F[2];
    $id =~ s/ /_/og;
    exists($ids{$F[0]}) || ($ids{$F[0]} = []);
    push @{ $ids{$F[0]} }, $id;
    $names{$id} = $F[0]; 
    $refnames{$F[0]} = $id if $F[3] =~ /scientific name/io;
    # print STDERR "n";
    &counter(++$n,\$c);
};
close(NI);
print STDERR "\n# Total NAMES: $n\n";
# print STDERR "\n",Data::Dumper->Dump([ \%names, \%refnames ],[ qw( *names *refnames) ]),"\n";

#
# READING NODES.DMP
# 1.  tax_id        -- node id in GenBank taxonomy database
# 2.  parent tax_id -- parent node id in GenBank taxonomy database
# 3.  rank          -- rank of this node (superkingdom, kingdom, ...) 
print STDERR "### READING NODES TREE\n";
open(TR, shift @ARGV);
($n, $c) = (0, 0);
while (<TR>) {
    our (@F);
    $_ =~ s/\t\|//og;
    @F = split /\t/o, $_;
    $F[2] =~ s/ /_/og;
    $tree{$F[0]} = [ $F[1], $F[2] ];
    # print STDERR "t";
    &counter(++$n,\$c);
};
close(TR);
print STDERR "\n# Total NODES: $n\n";
# print STDERR "\n";

#
# READING SPECIES/GENUS LIST
### 1. species/genus id (whitespaces as underscores)
# 1. SeqID
# 2. TaxID
print STDERR "### READING SEQ IDs LIST\n";
open(DT, shift @ARGV);
($n, $c) = (0, 0);
while (<DT>) {
    our ($sid, $tx, $nm, $f, $ori);
    ++$n;
    chomp;
    ($sid, $tx, undef) = split /\s+/, $_, 3;
    $nm = exists($refnames{$tx}) 
            ? $refnames{$tx}
            : ( exists($updatetx{$tx}) && exists($refnames{$updatetx{$tx}})
                ? $refnames{$updatetx{$tx}}
                : 32644 );  #32644 is the TaxonomyID for "unidentified species"
    $f = 1;
    $ori = $nm;
    printf STDERR "%08d %s : %d %s : %s : ", $n ,$sid, $tx, exists($refnames{$tx}) ? $refnames{$tx} : "---", $nm;
    # while (!exists($names{$nm}) && $f == 1) {
    # 	$f = $nm =~ s/_[^_]*$//;
    #	print STDERR " >>$nm<<";
    # };
    push @O, [ $nm, $sid, $tx ];
    exists($genus{$nm}) || do {
        $genus{$nm} = {};
        &taxonomy($nm,$genus{$nm},1);
    };
    print STDERR "\n";
};
close(DT);
# print STDERR "\n";

#
# PRINTING FULL TAXONOMY FOR GIVEN GENUS
print STDERR "### FILTERING TAXONOMY\n";
foreach my $ary (@O) {
    my ($i, $seqid, $txcode, @rnkids);
    ($i, $seqid, $txcode) = @$ary;
    print STDERR ">> $i $seqid $txcode\n";
    $names{$i} = "NA" unless exists($names{$i});
    print "$seqid\t$names{$i}\t$i";
    @rnkids = sort { $b->[0] <=> $a->[0] }
               map { [ split /:/, $_, 2 ] }
              keys %{ $genus{$i} };
    print "\tNA\:Unknown" if scalar(@rnkids) == 0;
    foreach my $ary (@rnkids) {
        my ($j,$d) = @$ary;
        my $J=$j.":".$d;
        print "\t$d\:$genus{$i}{$J}";
    };
    print "\n";
};
print STDERR "### DONE...\n";

exit(0);

##
sub counter() {
    my ($N, $C) = @_;
    my @D = qw( - / | \ + * );
    ($N %    500 == 0) && do {
        print STDERR ($$C > 0 ? "\b" : q{}).$D[$$C];
        $$C++;
    };
    ($N % 250000 == 0) && do {
        printf STDERR "[%08d]\r", $N;
    };
    $$C = 0 if $$C > 5;
} # counter

sub taxonomy() {
    my ($id,$ref,$lvl) = @_;
    ($id ne "root" && exists($names{$id})) || return;
    my $code = $names{$id};
    exists($tree{$code}) || return;
    my ($par,$rnk) = @{ $tree{$code} };
    exists($refnames{$par}) || ($refnames{$par} = "Unknown");
    my $ncode = $refnames{$par};
    $ref->{$lvl.":".$rnk} = "$id\:$code";
    # (my $tid    = $id   ) =~ s/ /_/og;
    # (my $trnk   = $rnk  ) =~ s/ /_/og;
    # (my $tncode = $ncode) =~ s/ /_/og;
    print STDERR " $code/$id/$rnk/$par/$ncode";
    # ($rnk eq 'phylum' || $ncode eq $id) && return;
    &taxonomy($ncode,$ref,++$lvl);
} # taxonomy
