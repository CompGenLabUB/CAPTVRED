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
#use global qw( :Counter );
$global::_verbose{RAW} = 1;
my $debug = 0; 
my $sampid;
my $kainames; 
my $out_prefix;  


GetOptions (
    'samp=s'            =>\$sampid,
    'K=s'               => \$kainames,
    'o|out_prefix=s'    => \$out_prefix,
    'debug'             => sub { $debug = 1 } 
);


my %tids   = ();   #Targets info Hash (HoL)
my$totrds = 0;
my $totspecs = 0;

print STDERR "#   STARTING... SAMPLES ID is $sampid   #  \n";
print STDERR "\n=========================================\n";
print STDERR "### RUNNING PARAMS ###  \n";
print STDERR "=========================================\n";
print STDERR "Kaiju outfile:\t$kainames\n";
print STDERR "Outfile prefix:\t$out_prefix\n";
print STDERR "=========================================\n\n";



print STDERR "### READING INFILE AND WRITING BY READ\n";
## Process kaiju names.out file to get info by sequence and taxon.
open(NMS, $kainames);
open(SOUT, ">", join( "_", $out_prefix, "taxonomysum_byread.tbl" ));
print SOUT "CONTIG_ID\tTAXON_ID\tSPECIES\tFAMILY\tKAIJU_SCORE";
while (<NMS>) {
    next if /^\s*$/o;
    chomp;
    my (@obs, @tax, $cid, $tid, $sco, $fam, $spc);
    @obs = split /\t/o, $_;

    $obs[0] eq "C" || next;
    $cid = $obs[1];
    $tid = $obs[2];
    $sco = $obs[3];
    @tax = split /;\s|;/o, $obs[7];
    $fam = $tax[1];
    $spc = $tax[2];


    $totrds ++;
    print SOUT "$cid\t$tid\t$spc\t$fam\t$sco\n";
    if ( exists($tids{$tid}) ) {

        $tids{$tid}[2] +=1;
        $tids{$tid}[3] = "$tids{$tid}[3];$cid";

    } else {

        $tids{$tid} = [$spc, $fam, 1, $cid];

    }
}
close(SOUT);

print STDERR "### WRITTING BY TAXON\n";
open(QOUT, ">", join( "_", $out_prefix, "taxonomysum_bytaxon.tbl" ));
foreach my $i (keys %tids){
    print QOUT join("\t", $i, @{$tids{$i}}), "\n";
    $totspecs ++;
}
close(QOUT);

print STDERR "### WRITTING STATS\n";
open(STSOUT, ">", join( ".", $out_prefix, "stats.out" ));
printf STSOUT "SAMPId:$sampid\tNREADS:$totrds\tNSPECS:$totspecs\n";
close(STSOUT);

exit(0);

=comment
    if ( $obs[0] eq "C" ) {

         #Query lines -> Save in hash:
         #%qids={ contig ID = [length, BestHitlen, BestHitCoverage, Reference sequence, taxon_approach, #reads], ContigID2 =[...], ... }
         #exists($qids{$obs[1]} ) && next;
         #print STDOUT $obs[1]."\n";
         (undef, undef, $rf, undef) = split /\|/o, $obs[13];  
            #refseq info is split to get only the seq ID;
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
        $spc='NA';
        $tid='NA';
      #  $spc=(defined($gn[2]) && $gn[2] ne '') ? $gn[2] : "NA";
      #  $spc =~ s/\W+/_/og; # cleaning up species name
      #  $spc =~ s/_+/_/og;
      #  $spc =~ s/(_\b|\b_)//og;
        for (my $i=3; $i < @gn; $i++) {
            ($gn[$i] =~ /^family/o) && do {
                      my @g = split /:/o, $gn[$i];
                      $fam = $g[1]; };
            ($gn[$i] =~ /^species/o) && do {
                      my @g = split /:/o, $gn[$i];
                      $spc = $g[1];
                      $tid = $g[2]; };
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
            push @{ $qids{$contigid} }, $tid, $spc, $fam if scalar(@{ $qids{$contigid} }) < 10;   ## Add species, family and taxonID information to queries HoL.
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
my $totsqs=$n;
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
my $totspecs=$n;
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
my $totrds=$n;
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

print STDERR "### WRITTING STATS\n";
open(STSOUT, ">", join( ".", $out_prefix, "stats.out" ));
#my $totrds=scalar keys %qids;
#my $totspecs=scalar keys %specs;
#my $totsqs=scalar keys %tids;
#my $seqsfromkit=
#my $specsfromkit=
printf STSOUT "SAMPId:$sampid\tNREADS:$totrds\tNSEQS:$totsqs\tNSPECS:$totspecs\n";
close(STSOUT);
exit(0);
=end
