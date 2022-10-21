#!/usr/bin/perl
#
# coverageblasttbl.pl
#
#   Calculates the overall coverage for all queries and targets
#   from BLAST output, accounting for overlaps among HSPs within
#   each single query and target sequences.
#
# ####################################################################
#
#        Copyright (C) 2013 - Josep Francesc ABRIL FERRANDO
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# ####################################################################
#
# $Id$
#
# USAGE:
#   coverageblasttbl.pl [options] blast.output6.tbl > coverage.tbl
#
# NCBI-BLAST2+ requested format (short tbl):
#
# -outfmt '6 qseqid qlen sseqid slen qstart qend sstart send length
#            score evalue bitscore pident nident mismatch positive gapopen gaps ppos qframe sframe qseq sseq'
#
use strict;
use warnings;
use Getopt::Long;
#
# VARS
my ($Qaaflg,$Taaflg,$program,$debugflg) = (undef,undef,undef,undef);
my %BLAST = (  # fixing coords to compute coverages...
	       #   Q  T
    'BLASTN'  => [ 0, 0 ], # ratio 1 nuc to 1 nuc
    'BLASTP'  => [ 0, 0 ], # ratio 1 aa  to 1 aa ~ 1 nuc to 1 nuc
    'BLASTX'  => [ 0, 1 ], # ratio 3 nuc to 1 aa
    'TBLASTN' => [ 1, 0 ], # ratio 1 aa  to 3 nuc
    'TBLASTX' => [ 0, 0 ]  # ratio 1 nuc to 1 nuc
    );

GetOptions( "Qaa"    => \$Qaaflg,
            "Taa"    => \$Taaflg,
            "prog=s" => \$program,
            "debug"  => \$debugflg
    );

defined($Qaaflg)  || ($Qaaflg = 0);
defined($Taaflg)  || ($Taaflg = 0);
$debugflg = defined($debugflg) ? 1 : 0;
defined($program) && do {
    exists($BLAST{uc($program)}) || ($program = 'blastn');
    ($Qaaflg,$Taaflg) = @{ $BLAST{uc($program)} };
    print STDERR "## PARSING TBL FILE from ",uc($program)," [$Qaaflg$Taaflg]\n";
};

my $CSTR = ("\b" x 9).'%09d';
my $RPTC = 1000;

my %Q = ();
my %T = ();

#
print STDERR "## PARSING TBL FILE... 000000000";
my $z = 0;
while (<>) {
    my (@F, $qi, $qs, $qe, $qr, $ti, $ts, $te, $tr, $E, $S, $L, $ql, $tl,
        $saln, $pidn, $ppos);

    next if /^\#/o;
    next if /^\s*$/o;
    
    chomp;

    @F = split /\s+/o, $_; 

    # GFF -> 0,3,4,6,9,10,11,12,15,21,29
    # TBL       -> 0,3,4,7,8,9,10,11,12,14,15
    # 6 qseqid qgi qacc qlen sseqid sgi sacc slen qstart qend sstart send length score evalue bitscore pident nident mismatch positive gapopen gaps ppos qframe sframe qseq sseq
    # SHORT-TBL -> 0,1,2,3,4,5,6,7,8,10,11 (+ 9,12,18 for output purposes only)
    # 6 qseqid qlen sseqid slen qstart qend sstart send length score evalue bitscore pident nident mismatch positive gapopen gaps ppos qframe sframe qseq sseq
    ($qi, $ql, $ti, $tl, $qs, $qe, $ts, $te, $L, $E, $S, $saln, $pidn, $ppos) = @F[0,1,2,3,4,5,6,7,8,10,11,9,12,18];
    $qr = ($qs > $qe) ? '-' : '+';
    $tr = ($ts > $te) ? '-' : '+';
    $qs > $qe && (($qs, $qe) = ($qe, $qs));
    $ts > $te && (($ts, $te) = ($te, $ts));
    # $ql = $qe - $qs + 1;
    # $tl = $te - $ts + 1;

    $Qaaflg && do {
        ($qs, $qe) = ($qs * 3 - 2, $qe * 3);
        $ql = $ql * 3;
    };
    $Taaflg && do {
        ($ts, $te) = ($ts * 3 - 2, $te * 3);
        $tl = $tl * 3;
    };

    exists($Q{$qi}) || do {
	$Q{$qi} = {
	    "L" => [ $ql, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
                   # length NnucF % NnucR % NnucB % MaxFR % BestHit % alnsco idnpct pospct
	    "+" => [],
	    "-" => [],
	    "." => [],
            "b" => 'NA', # best hit ID
            "B" => []    # to store best hit coords
	};
    };
    exists($T{$ti}) || do {
	$T{$ti} = {
	    "L" => [ $tl, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
	    "+" => [],
	    "-" => [],
	    "." => [],
            "b" => 'NA',
            "B" => []
	};
    };

    push @{ $Q{$qi}{$qr} }, [ $qs, $qe, $ti, $S, $E, $qr.$tr, $saln, $pidn, $ppos ];
    push @{ $T{$ti}{$tr} }, [ $ts, $te, $qi, $S, $E, $tr.$qr, $saln, $pidn, $ppos ];

    # print join("\t", $qi, $qs, $qe, $qr, $ti, $ts, $te, $tr, $ql, $tl, $S, $E),"\n";

    print STDERR sprintf($CSTR,$z) if (++$z % $RPTC == 0);
};
print STDERR sprintf($CSTR,$z) if ($z % $RPTC != 0); print STDERR " lines.\n";

#
print STDERR "## PROCESSING BEST HIT for QUERIES... 000000000";
$z = 0;
foreach my $q (keys %Q) {
    $Q{$q}{'b'} = &get_besthit($Q{$q}{'+'},$Q{$q}{'-'},$Q{$q}{'B'},\%T);
    print STDERR sprintf($CSTR,$z) if (++$z % $RPTC == 0);
};
print STDERR sprintf($CSTR,$z) if ($z % $RPTC != 0); print STDERR " queries.\n";

#
print STDERR "## PROCESSING BEST HIT for TARGETS... 000000000";
$z = 0;
foreach my $q (keys %T) {
    $T{$q}{'b'} = &get_besthit($T{$q}{'+'},$T{$q}{'-'},$T{$q}{'B'},\%Q);
    print STDERR sprintf($CSTR,$z) if (++$z % $RPTC == 0);
};
print STDERR sprintf($CSTR,$z) if ($z % $RPTC != 0); print STDERR " targets.\n";

#
print STDERR "## PROCESSING OVERLAPS for QUERIES... 000000000";
$z = 0;
foreach my $q (keys %Q) {
    &process_overlaps($Q{$q}{'+'});
    &process_overlaps($Q{$q}{'-'});
    push @{ $Q{$q}{'.'} }, @{ $Q{$q}{'+'} }, @{ $Q{$q}{'-'} };
    &process_overlaps($Q{$q}{'.'});
    &process_overlaps($Q{$q}{'B'});
    print STDERR sprintf($CSTR,$z) if (++$z % $RPTC == 0);
};
print STDERR sprintf($CSTR,$z) if ($z % $RPTC != 0); print STDERR " queries.\n";

#
print STDERR "## PROCESSING OVERLAPS for TARGETS... 000000000";
$z = 0;
foreach my $q (keys %T) {
    &process_overlaps($T{$q}{'+'});
    &process_overlaps($T{$q}{'-'});
    push @{ $T{$q}{'.'} }, @{ $T{$q}{'+'} }, @{ $T{$q}{'-'} };
    &process_overlaps($T{$q}{'.'});
    &process_overlaps($T{$q}{'B'});
    print STDERR sprintf($CSTR,$z) if (++$z % $RPTC == 0);
};
print STDERR sprintf($CSTR,$z) if ($z % $RPTC != 0); print STDERR " targets.\n";

#
print STDERR "## PROCESSING COVERAGE for QUERIES... 000000000";
$z = 0;
foreach my $q (keys %Q) {
    $Q{$q}{'L'}[1] = &get_coverage($Q{$q}{'+'});
    $Q{$q}{'L'}[3] = &get_coverage($Q{$q}{'-'});
    $Q{$q}{'L'}[5] = &get_coverage($Q{$q}{'.'});
    $Q{$q}{'L'}[7] = ($Q{$q}{'L'}[3] > $Q{$q}{'L'}[1]) ? $Q{$q}{'L'}[3] : $Q{$q}{'L'}[1];
    $Q{$q}{'L'}[9] = &get_coverage($Q{$q}{'B'});
    if ($Q{$q}{'L'}[0] != 0) {
	$Q{$q}{'L'}[2] = sprintf("%6.3f", ($Q{$q}{'L'}[1] / $Q{$q}{'L'}[0]) * 100);
	$Q{$q}{'L'}[4] = sprintf("%6.3f", ($Q{$q}{'L'}[3] / $Q{$q}{'L'}[0]) * 100);
	$Q{$q}{'L'}[6] = sprintf("%6.3f", ($Q{$q}{'L'}[5] / $Q{$q}{'L'}[0]) * 100);
	$Q{$q}{'L'}[8] = sprintf("%6.3f", ($Q{$q}{'L'}[7] / $Q{$q}{'L'}[0]) * 100);
	$Q{$q}{'L'}[10] = sprintf("%6.3f", ($Q{$q}{'L'}[9] / $Q{$q}{'L'}[0]) * 100);
    } else {
	$Q{$q}{'L'}[2] = $Q{$q}{'L'}[4] = $Q{$q}{'L'}[6] = $Q{$q}{'L'}[8] = $Q{$q}{'L'}[10] = "NA";
    };

    print STDOUT join("\t", "Q", $q, @{ $Q{$q}{'L'} }, $Q{$q}{'b'}),"\n";

    print STDERR sprintf($CSTR,$z) if (++$z % $RPTC == 0);
};
print STDERR sprintf($CSTR,$z) if ($z % $RPTC != 0); print STDERR " queries.\n";

#
print STDERR "## PROCESSING COVERAGE for TARGETs... 000000000";
$z = 0;
foreach my $q (keys %T) {

    $T{$q}{'L'}[1] = &get_coverage($T{$q}{'+'});
    $T{$q}{'L'}[3] = &get_coverage($T{$q}{'-'});
    $T{$q}{'L'}[5] = &get_coverage($T{$q}{'.'});
    $T{$q}{'L'}[7] = ($T{$q}{'L'}[3] > $T{$q}{'L'}[1]) ? $T{$q}{'L'}[3] : $T{$q}{'L'}[1];
    $T{$q}{'L'}[9] = &get_coverage($T{$q}{'B'});
    if ($T{$q}{'L'}[0] != 0) {
        $T{$q}{'L'}[2] = sprintf("%6.3f", ($T{$q}{'L'}[1] / $T{$q}{'L'}[0]) * 100);
        $T{$q}{'L'}[4] = sprintf("%6.3f", ($T{$q}{'L'}[3] / $T{$q}{'L'}[0]) * 100);
        $T{$q}{'L'}[6] = sprintf("%6.3f", ($T{$q}{'L'}[5] / $T{$q}{'L'}[0]) * 100);
        $T{$q}{'L'}[8] = sprintf("%6.3f", ($T{$q}{'L'}[7] / $T{$q}{'L'}[0]) * 100);
        $T{$q}{'L'}[10] = sprintf("%6.3f", ($T{$q}{'L'}[9] / $T{$q}{'L'}[0]) * 100);
    } else {
        $T{$q}{'L'}[2] = $T{$q}{'L'}[4] = $T{$q}{'L'}[6] = $T{$q}{'L'}[8] = $T{$q}{'L'}[10] = "NA";
    };

    print STDOUT join("\t", "T", $q, @{ $T{$q}{'L'} }, $T{$q}{'b'}),"\n";

    print STDERR sprintf($CSTR,$z) if (++$z % $RPTC == 0);
};
print STDERR sprintf($CSTR,$z) if ($z % $RPTC != 0); print STDERR " targets.\n";

#
exit(0);

#
sub get_besthit($$$$) {
    my ($F,$R,$B,$M) = @_;
    my $str = '';
    my %K = ();
    my @T = ();

    foreach my $a (@$F, @$R) {
        my ($s,$e,$i,$S,$E,$r, $saln,$pidn,$ppos) = @$a;
        exists($K{$i}) || do {
            $K{$i} = {
                "S" => $S,
                "E" => $E,
                "R" => $r,
                "P" => [ [ $s, $e ] ],
                "X" => [ $saln, $pidn, $ppos ]
            };
            next;
        };
        $K{$i}{"S"} > $S || ($K{$i}{"S"} = $S, $K{$i}{"E"} = $E, $K{$i}{"R"} = $r);
        push @{ $K{$i}{"P"} }, [ $s, $e ];
    };

    @T = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] || $a->[2] <=> $b->[2] } # greater bitscore or smaller evalue
         map { [ $_, $K{$_}{"S"}, $K{$_}{"E"} ] } keys %K;

    $str = $T[0]; # $str = $T[0];
    @$B = map { [ $_->[0], $_->[1] ] } @{ $K{$str}{"P"} };
    $str .= " ".$M->{$T[0]}{"L"}[0]." ".$K{$str}{"S"}." ".
                $K{$str}{"E"}." ".$K{$str}{"R"}." ".join(" ",@{$K{$str}{"X"}});

    # reducing memory footprint
    @$F = map { [ $_->[0], $_->[1] ] } @$F; 
    @$R = map { [ $_->[0], $_->[1] ] } @$R; 
    undef %K;
    undef@T;

    &printcoords("BHE",$B) if $debugflg;
    
    return $str;
} # get_besthit

sub process_overlaps($) {
    my $A = shift;
    my @T = ();

    scalar(@$A) > 1 || return;

    # ensuring coords sorted
    @$A = sort { $a->[0] <=> $b->[0] || $b->[1] <=> $a->[1] } @$A;

    &printcoords("POS",$A) if $debugflg;

    my ($m, $n) = (undef,undef);
    foreach my $a (@$A) {

	my ($s,$e) = @$a;

	defined($m) || ($m=$s, $n=$e, next);
	# $n < $e || next;
	$s > ($n + 1) && do { # otherwise they can be merged into a single hsp
	    push @T, [ $m, $n ];
	    ($m,$n) = ($s,$e);
	    next;
	};
	$n < $e && ($n = $e);
    };
    push @T, [ $m, $n ];

    @$A = @T;

    &printcoords("POE",$A) if $debugflg;

} # process_overlaps

sub get_coverage($) {
    my $A = shift;
    my $sum = 0;
    foreach my $a (@$A) {
        my ($s,$e) = @$a;
        $sum += ($e - $s + 1);
    };
    return $sum;
} # get_coverage

sub printcoords($$) {
    my ($l,$A) = @_;

    printf STDERR "$l : ";

    foreach my $a (@$A) {

        printf STDERR "%s-%s ", $a->[0], $a->[1];

    };

    printf STDERR "\n";

} # printcoords
