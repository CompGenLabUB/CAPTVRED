#!/usr/bin/perl
#
# 
#
#
#
#
#
#
# USAGE:
# 
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use global qw( :Counter );
$global::_verbose{RAW} = 1;
$_cntN = 50000;
my $debug = 0; 
my $datab_info; # = $ARGV[0];
my $blastcov;  # = $ARGV[1];
my $clas_ids;  # = $ARGV[2];
my $unclas_ids;  # = $ARGV[3];
my $out_prefix;  # = $ARGV[4];
my $mincov;
my $petbl;   ## Tabular file obtained with samtools idxstats: each row is a contig, 3rd column is #pe-reads mapped 
my $setbl;   ## Tabular file obtaiden with samtools idxstats: each row is a contig, 3rd column is #se-reads mapped 

GetOptions (
    'i|db_info=s'       => \$datab_info,
    'B=s'               => \$blastcov,
    'mincov:i'          => \$mincov,
    'c|class_ids:s'     => \$clas_ids,
    'u|unclass_ids:s'   => \$unclas_ids,
    'o|out_prefix=s'    => \$out_prefix,
    'pe=s'              => \$petbl,
    'sg=s'              => \$setbl,
    'debug'             => sub { $debug = 1 } 
);


my %tids   = ();   #Targets info Hash (HoL)
my %qids   = ();   #Queries info Hash (HoL)
my %rids   = ();   #Refseqs - contigs relation Hash (HoL)  
                     #%rids = { refseq1=[contig1, ..., contigN], refseq1=[contig1, ..., contigN], ...  }
#my %kids   = ();
my %specs   = ();   #Species info Hash (HoL)
#my %tcount = ();
my %kcount = ();
#my %bcount = ();
my %apf = ();   #Approach used to found the hit (Kaiju or blast) + score. (HoL)  %apf={ contig1=[ "K", score], contig1=[ "B", "NA"] } 

# 3 dicts used structure will be:
  #%qids={ contig ID = [length, BestHitlen, BestHitCoverage, Reference sequence, taxon_approach], ContigID2 =[...], ... }

defined($mincov) || ($mincov=70);
my $ci=defined($clas_ids) ? $clas_ids : "Not provided";
my $ui=defined($unclas_ids) ? $unclas_ids : "Not provided";
my $pf=defined($petbl) ? $petbl : "Not provided";
my $sf=defined($setbl) ? $setbl : "Not provided";

print STDERR "\n=========================================\n";
print STDERR "### RUNNING PARAMS ###  \n";
print STDERR "=========================================\n";
print STDERR "DB  info  file:\t$datab_info\n";
print STDERR "Coverage  file:\t$blastcov\n";
print STDERR "Min.  coverage:\t$mincov\n";
print STDERR "Class. ids  fl:\t$ci \n";
print STDERR "Unc.  ids   fl:\t$ui \n";
print STDERR "#pe-reads map :\t$pf\n";
print STDERR "#se-reads map :\t$sf\n";
print STDERR "Outfile prefix:\t$out_prefix\n";
print STDERR "=========================================\n\n";

# --------------------------------- Under construction ----#
if (defined($clas_ids)){
    print STDERR "### READING CLASSIFIED IDS FILE\n";
    my @inf;
    open(CL, $clas_ids);
    while (<CL>) {
        chomp;
        # @inf=  ["K", "NA"];
        # $apf{$_} = $@inf;
        @apf{$_} = ["K", "NA"];  #"NA" must be substituted by kaiju score
    }
    close(CL);
}

if (defined($unclas_ids)){
    print STDERR "### READING UNCLASSIFIED IDS FILE\n";
    open(UC, $unclas_ids);
    while (<UC>) {
        chomp;
        # @inf= ["B", "NA"];
        # $apf{$_} = $@inf;
        @apf{$_} = ["B", "NA"];
    }
    close(UC);
}
# -------------------------- End of Under construction ----#

my %rcounts = ();
if (defined($petbl)){
    print STDERR "### READING MAPPED PE READS FILE\n";
    my @info;
    open(PE, $petbl);
    while (<PE>) {
        next if /^\s*$/o;
        chomp;
        @info = split /\t/o, $_;
        if (exists($rcounts{$info[0]})) {
            print STDERR "# WARNING # $info[0] is already there...\n";
        } else {
            $rcounts{$info[0]} = 0;
        };
        $rcounts{$info[0]} += $info[2];
    }
    close(PE);
}

if (defined ($setbl)){
    print STDERR "### READING MAPPED SE READS FILE\n";
    my @info;
    open(SE, $setbl);
    while (<SE>) {
        next if /^\s*$/o;
        chomp;
        @info = split /\t/o, $_;
        exists($rcounts{$info[0]}) || ($rcounts{$info[0]} = 0);
        $rcounts{$info[0]} += $info[2];
    }
    close(SE);
}

print STDERR Dumper \%rcounts if $debug;

## Process blast coverage out file to get info by sequence and complete info by genome.
print STDERR "### READING BLAST COVERAGE FILE\n";
open(COV, $blastcov);
$c = "B"; $n = 0;
# my %qrcnts;    # n reads maped in each query
while (<COV>) {
    next if /^\s*$/o;
    chomp;
    my (@obs, $cid, $rf, $tag, $ks, $nrd, @t, @r);
    @obs = split /\t/o, $_;
    $cid = $obs[1];

    if ( $obs[0] eq "Q" ) {

         #Query lines -> Save in hash:
         #%qids={ contig ID = [length, BestHitlen, BestHitCoverage, Reference sequence, taxon_approach, #reads], ContigID2 =[...], ... }
         #exists($qids{$obs[1]} ) && next;
         #print STDOUT $obs[1]."\n";
         (undef, undef, $rf, undef) = split /\|/o, $obs[13];   #refseq info is split to get only the seq ID.
         $tag = (defined($clas_ids) ? $apf{$cid}[0] : "B");
         $ks  = (defined($clas_ids) ? $apf{$cid}[1] : "NA");
         #exists ($rcounts{$obs[1]}) || ($rcounts{$obs[1]} = 1);  ## singletons
         $nrd = (exists($rcounts{$cid}) ? $rcounts{$cid} : 1);
         $qids{$cid} = [ $tag, $obs[2], $obs[11], $obs[12], $rf, $ks, $nrd ];
                     #  0.TAG   1.len   2.BHlen   3.BHCov  4.RefSqId  5.KaiScore  6.n.reads 
         exists($rids{$rf}) || ($rids{$rf} = { '##SUM##' => 0, '##NCTGS##' => 0 });
         $rids{$rf}{$cid} = $nrd;
         $rids{$rf}{'##SUM##'} += $nrd;
         $rids{$rf}{'##NCTGS##'}++;
         # push @{ $rids{$rf} }, $cid;
         # mapped reads summatory:
       #  exists($qrcnts{$rf[2]}) || ($qrcnts{$rf[2]} = 0);
         #exists ($rcounts{$obs[1]}) || ($rcounts{$obs[1]} = 1);  ## singletons
       #  $qrcnts{$rf[2]} += $rcounts{$obs[1]};
         # . . . 
         #$kids{$rf[2]}{$obs[1]}=$rcounts{$obs[1]};
         #$kids{$rf[2]}{"tot"}+=$rcounts{$obs[1]};
        # exists($bcount{$rf[2]}) || ($bcount{$rf[2]} = 0);
        # $bcount{$rf[2]}++;

    } else { # if ( $obs[0] eq "T") {

        # Target lines -> add to %tids  legnth (in nts) and coverage of the total coverage:
        @t = split /\|/o, $cid;
        @r = split /\s/o, $obs[13];
        # my $sid=$t[2];

        $tids{$t[2]} = [ $obs[2], $obs[7], $obs[8], $r[6] ];
            # SEQID      SEQLEN   NUCALN   COVPCT   BHIDENT
    }
} continue {
    $n++;
    &counter($n,$c) if ($n % 1000 == 0);
}
&counter_end($n);
close(COV);

# print STDOUT Dumper \%rids;

print STDERR "### READING GENOMES INFO FILE\n";
open(GIN, $datab_info);
my ($sid, $tid, $spc, $fam,$k, $nreads, $ncontigs);
$c = "G"; $n = 0;
while (<GIN>) {
    next if /^\s*$/o;
    chomp;
    #All reference sequences-> Save in hash:
         
    my @gn = split /\t/o, $_;
    $sid=$gn[0];
    exists($tids{$sid}) && do {
        $fam = 'NA';
        $tid=$gn[1]; # it can be "NA"
        $spc=(defined($gn[2]) && $gn[2] ne '') ? $gn[2] : "NA";
        $spc =~ s/\W+/_/og; # cleaning up species name
        $spc =~ s/_+/_/og;
        $spc =~ s/(_\b|\b_)//og;
        for (my $i=3; $i < @gn; $i++) {
            if ($gn[$i] =~ /^family/o) {
               my @g = split /:/o, $gn[$i];
               $fam = $g[1];
               last;
            };
        };
        # $nreads=$qrcnts{$sid};  # number of reads mapped 

       # push @{$tids{$sid}}, (exists($bcount{$sid}) ? $bcount{$sid} : 0);
        #push @{$tids{$sid}}, (exists($kcount{$tid}) ? $kcount{$tid} : 0);
        
        ## Add n.contigs i n.reads to targets HoL.
        $ncontigs = exists($rids{$sid}{'##NCTGS##'}) ? $rids{$sid}{'##NCTGS##'} : 0;
        $nreads   = exists($rids{$sid}{'##SUM##'})   ? $rids{$sid}{'##SUM##'}   : 0;
        # $nreads=0;
        foreach my $contigid (keys %{ $rids{$sid} }) {
            next if $contigid =~ /##(SUM|NCTGS)##/;
        #     $nreads += $qids{$contigid}[6];  ## nReads sumatory
            push @{ $qids{$contigid} }, $tid, $spc, $fam;   ## Add species, family and taxonID information to queries HoL.
        };
        push @{$tids{$sid}}, $tid, $spc, $fam, $ncontigs, $nreads;

        print STDERR $sid." --> >".$ncontigs."<\n" if $debug;
        ## Add n.reads of each contig to queries HoH.
     #   foreach my $contigid (keys %{ $kids{$sid} }) {    #sid=refseq_id
            # if ($qids{$key}[4] eq $sid){
      #          push @{ $qids{$contigid} }, $tid, $spc, $fam, $kids{$sid}{$contigid};
            # }                  #Taxonid    specie    family    nreads_maped_to_the_contig
       #}
    }; # exist $tid{$sid}
    
} continue {
    $n++;
    &counter($n,$c) if ($n % 1000 == 0);
}; # while GIN
&counter_end($n);
close(GIN);


print STDERR "### WRITTING GENOMES INFORMATION\n";
open ( ROUT, ">", join( "_", $out_prefix, "taxonomysum_bysequence.tbl" ));
print ROUT "SEQUENCE_ID\tTAXON_ID\tSPECIES\tFAMILY\tSEQ_LENGTH\tNUCS_ALN\tCOVERAGE_PCT\tBHIT_IDENTITY\tBESTHSP_COUNT\tNREADS_MAPED\n";
$c = "T"; $n = 0; 
my $m = 0; my %SC = ();
my (@k,$kn,$taxon,$coverage,$pidentity,$nconts,$nreadsmp,$i,$info);
my %FNDTID = ();
foreach my $j (keys %tids){  # 0.SEQLEN 1.NUCALN 2.COVPCT 3.BHIDENT 4.TXNID 5.SPC 6.FAM 7.NCONTIGS 8.NREADS
#    $kn = scalar @{ $tids{$j} } - 1;
        $taxon    =$tids{$j}[4];
        $coverage =$tids{$j}[2];
        $pidentity=$tids{$j}[3];
        $nconts   =$tids{$j}[7];
        $nreadsmp =$tids{$j}[8];
    exists($SC{$taxon}) || ($SC{$taxon} = [ 0, 0, 0, 0 ]);
    $SC{$taxon}[0] += $nconts;
    $SC{$taxon}[1] += $nreadsmp;
    $m++;
    ($coverage >= $mincov) || next;
    ($nconts > 0) || next;
        $FNDTID{$j}=$nconts;
        $SC{$taxon}[2] += $nconts;
        $SC{$taxon}[3] += $nreadsmp;
   # if (@{$tids{$j}}[4] > 0 ) { 
        print ROUT join("\t", $j, ( map { defined($tids{$j}[$_]) ?  $tids{$j}[$_] : "NA" } (4,5,6,0,1,2,3,7,8) )), "\n";
        if (exists($specs{$taxon})) {
            #Update coverage info:
            $coverage < $specs{$taxon}[2] && ($specs{$taxon}[2]=$coverage); #cov min
            $coverage > $specs{$taxon}[3] && ($specs{$taxon}[3]=$coverage); #cov max
            $specs{$taxon}[4] += $coverage; # add cov to coverage summatory
            #Update identity perc info:
            $pidentity < $specs{$taxon}[5] && ($specs{$taxon}[5]=$pidentity); #pid min
            $pidentity > $specs{$taxon}[6] && ($specs{$taxon}[6]=$pidentity); #pid max
            $specs{$taxon}[7]  += $pidentity; # add pid to Pident sumatory
            #Update counters:
            $specs{$taxon}[8]++; # add +1 to refseqs sumatory
            $specs{$taxon}[9]  += $nconts;  # add n contigs to NContigs summatory
            $specs{$taxon}[10] += $nreadsmp;  # add n contigs to NContigs summatory
            # Update info field:
            $info=";SEQ:$j,LEN:$tids{$j}[0],NUCALN:$tids{$j}[1],COV:$coverage,BHIDENT:$pidentity,N:$nconts";
            $specs{$taxon}[11].=$info;  # Update info field.
        } else {
            $info="SEQ:$j,LEN:$tids{$j}[0],NUCALN:$tids{$j}[1],COV:$coverage,BHIDENT:$pidentity,N:$nconts";
            $specs{$taxon} = [ $tids{$j}[5], $tids{$j}[6], $coverage, $coverage, $coverage, $pidentity, $pidentity, $pidentity, 1, $nconts, $nreadsmp, $info ];
        }
    # }
    #print ROUT join("\t", $j,  $tids{$j}[5], $tids{$j}[6], $tids{$j}[7], $tids{$j}[0], $tids{$j}[1], $tids{$j}[2], $tids{$j}[3], $tids{$j}[4] ), "\n";
    $n++;
    &counter($n,$c) if ($n % 100 == 0);
}; # foreach $j
&counter_end($n);
close(ROUT);
$debug && do {
    my @sc = ( 0, 0, 0, 0 );
    printf STDERR "#># %-9s %9s %9s %9s %9s\n", qw/ TXN Ncnt Nrds Ncnt* Nrds* /;  
    foreach my $kk (sort { $a <=> $b } keys %SC) {
        printf STDERR "#># %-9s %9d %9d %9d %9d\n", $kk, @{$SC{$kk}};
        $sc[0] += $SC{$kk}[0];
        $sc[1] += $SC{$kk}[1];
        $sc[2] += $SC{$kk}[2];
        $sc[3] += $SC{$kk}[3];
    };
    printf STDERR "#># %-9s %9d %9d %9d %9d : ALL taxons %d FILTered taxons %d @sc\n", "TOTAL", @sc, $m, $n;  
};

print STDERR "### WRITTING SPECIES INFORMATION\n";
open(SOUT, ">", join( "_", $out_prefix, "taxonomysum_byspecie.tbl" ));
print SOUT "TAXON_ID\tSPECIES\tFAMILY\tMIN_COV\tMAX_COV\tMEAN_COV\tMIN_PID".
           "\tMAX_PID\tMEAN_PID\tN_SEQS\tN_CONTIGS\tNREADS_MAPED\tINFO\n";
$c = "S"; $n = 0;
foreach my $sp (keys %specs){
    $specs{$sp}[4] = $specs{$sp}[4]/$specs{$sp}[8];
    $specs{$sp}[7] = $specs{$sp}[7]/$specs{$sp}[8];
    print SOUT "$sp\t".join( "\t", @{ $specs{$sp} })."\n";
    $n++;
    &counter($n,$c) if ($n % 10 == 0);
}; # foreach $sp
&counter_end($n);
close(SOUT);

print STDERR "### WRITTING QUERIES INFORMATION\n";
open(QOUT, ">", join( "_", $out_prefix, "taxonomysum_byread.tbl" ));
print QOUT "READ_ID\tTAG\tREAD_LENGTH\tBESTHIT_LEN\tBESTHIT_COVERAGE\tREFSEQ_ID".
           "\tKAIJU_SCORE\tNREADS_MAPED\tTAXONID\tSPECIES\tFAMILY\n";
$c = "Q"; $n = 0;
my %found = (); $m = 0;
foreach my $i (keys %qids){
    #scalar @{$qids{$i}} == 10 || next;
 #   if ( @{$qids{$i}} != 10 ) { print STDOUT join("\t", $i, @{$qids{$i}}),"\n"; }
    #@{$qids{$i}}[3] >= $mincov || next;
    my $tid = $qids{$i}[4];
    $m++;
    exists($found{$tid}) || ($found{$tid} = [ 0, 0, 0, 0 ]);
    $found{$tid}[0]++;
    $found{$tid}[1]+=$qids{$i}[6];
    exists($FNDTID{$tid}) || next;
    $found{$tid}[2]++;
    $found{$tid}[3]+=$qids{$i}[6];
    print QOUT join("\t", $i, @{$qids{$i}}), "\n";
          # ID   # TAG   len    BHlen   BHCov  RefSqId  KaiScore NREADS TXNID SPC FAM 
    $n++;
    &counter($n,$c) if ($n % 1000 == 0);
}; # foreach $i
&counter_end($n);
close(QOUT);
$debug && do {
    my @sc = ( 0, 0, 0, 0 );
    printf STDERR "#*# %-9s %9s %9s %9s %9s\n", qw/ RID Rcnt Rrds Rcnt* Rrds* /;  
    foreach my $tx (keys %found) {
        printf STDERR "#*# %-9s %9d %9d %9d %9d\n", $tx, @{$found{$tx}};
        $sc[0] += $found{$tx}[0];
        $sc[1] += $found{$tx}[1];
        $sc[2] += $found{$tx}[2];
        $sc[3] += $found{$tx}[3];
    };
    printf STDERR "#*# %-9s %9d %9d %9d %9d : ALL %d FILTered %d\n", "TOTAL", @sc, $m, $n;  
};

exit(0);
