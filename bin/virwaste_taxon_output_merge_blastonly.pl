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

my $datab_info; # = $ARGV[0];
my $blastcov;  # = $ARGV[1];
my $clas_ids;  # = $ARGV[2];
my $unclas_ids;  # = $ARGV[3];
my $out_prefix;  # = $ARGV[4];
my $mincov;

GetOptions (
    'i|db_info=s'       => \$datab_info,
    'B=s'               => \$blastcov,
    'mincov:i'          => \$mincov,
    'c|class_ids:s'     => \$clas_ids,
    'u|unclass_ids:s'   => \$unclas_ids,
    'o|out_prefix=s'  => \$out_prefix
);


my %tids   = ();
my %qids   = ();   
my %kids   = ();
my %specs   = ();
#my %tcount = ();
my %kcount = ();
my %bcount = ();
my %apf = ();

# 3 dicts used structure will be:
  #%qids={ contig ID = [length, BestHitlen, BestHitCoverage, Reference sequence, taxon_approach], ContigID2 =[...], ... }

defined( $mincov ) || ( $mincov=70);
my $ci=defined($clas_ids) ? $clas_ids : "Not provided";
my $ui=defined($unclas_ids) ? $unclas_ids : "Not provided";

print STDERR "\n=========================================\n";
print STDERR "### RUNNING PARAMS ###  \n";
print STDERR "=========================================\n";
print STDERR "DB  info  file:\t$datab_info\n";
print STDERR "Coverage  file:\t$blastcov\n";
print STDERR "Min.  coverage:\t$mincov\n";
print STDERR "Class. ids  fl:\t$ci \n";
print STDERR "Unc.  ids   fl:\t$ui \n";
print STDERR "Outfile prefix:\t$out_prefix\n";
print STDERR "=========================================\n\n";

if (defined ($clas_ids)){
    print STDERR "### READING CLASSIFIED IDS FILE\n";
    my @inf;
    open(CL, $clas_ids);
    while (<CL>) {
        chomp;
        # @inf=  ["K", "NA"];
        # $apf{$_} = $@inf;
        @apf{$_} = ["K", "NA"];
    }
    close(CL);
}

if (defined ($unclas_ids)){
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


## Process blast coverage out file to get info by sequence and complete info by genome.
print STDERR "### READING BLAST COVERAGE FILE\n";
open(COV, $blastcov);
$c = "B"; $n = 0;
while (<COV>) {
    next if /^\s*$/o;
    chomp;
    my @obs = split /\t/o, $_;

    if ( $obs[0] eq "Q" ) {

         #Query lines -> Save in hash:
         #%qids={ contig ID = [length, BestHitlen, BestHitCoverage, Reference sequence, taxon_approach], ContigID2 =[...], ... }
         #exists($qids{$obs[1]} ) && next;
         my @rf = split /\|/o, $obs[13];   #refseq info is split to get only the seq ID.
         my $tag=(defined ($clas_ids) ? $apf{$obs[1]}[0] : "B");
         my $ks=(defined ($clas_ids) ? $apf{$obs[1]}[1] : "NA");
         $qids{$obs[1]} = [ $tag, $obs[2], $obs[11], $obs[12], $rf[2], $ks ];
                       #    TAG   len      BHlen     BHCov   RefSqId  KaiScore
         $kids{$rf[2]}{$obs[1]}++;
         exists($bcount{$rf[2]}) || ($bcount{$rf[2]} = 0);
         $bcount{$rf[2]}++;

    } else { # if ( $obs[0] eq "T") {

        # Target lines -> add to %tids  legnth (in nts) and coverage of the total coverage:
        my @t= split /\|/o, $obs[1];
        my @r= split /\s/o, $obs[13];
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



print STDERR "### READING GENOMES INFO FILE\n";
open(GIN, $datab_info);
my ($sid, $tid, $spc, $fam,$k);
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
    
        push @{$tids{$sid}}, (exists($bcount{$sid}) ? $bcount{$sid} : 0);
        #push @{$tids{$sid}}, (exists($kcount{$tid}) ? $kcount{$tid} : 0);
        push @{$tids{$sid}}, $tid, $spc, $fam;
        $k = 1;
        foreach my $key (keys %{ $kids{$sid} }) {
            # if ($qids{$key}[4] eq $sid){
                push @{ $qids{$key} }, $tid, $spc, $fam;
            # } 
       }
    }; # exist $tid{$sid}
    
} continue {
    $n++;
    &counter($n,$c) if ($n % 1000 == 0);
}; # while GIN
&counter_end($n);
close(GIN);

print STDERR "### WRITTING QUERIES INFORMATION\n";
open ( QOUT, ">", join( "_", $out_prefix, "taxonomysum_byread.tbl" ));
print QOUT "READ_ID\tTAG\tREAD_LENGTH\tBESTHIT_LEN\tBESTHIT_COVERAGE\tREFSEQ_ID\tKAIJU_SCORE\tTAXONID\tSPECIES\tFAMILY\n";
$c = "Q"; $n = 0;
foreach my $i (keys %qids){
    scalar @{$qids{$i}} == 9 || next;
    @{$qids{$i}}[3] >= $mincov || next; 
    print QOUT join("\t", $i, @{$qids{$i}}), "\n";
                             # TAG   len    BHlen   BHCov  RefSqId  KaiScore  TXNID SPC FAM
} continue {
    $n++;
    &counter($n,$c) if ($n % 1000 == 0);
}; # while GIN
&counter_end($n);
close(QOUT);

print STDERR "### WRITTING GENOMES INFORMATION\n";
open ( ROUT, ">", join( "_", $out_prefix, "taxonomysum_bysequence.tbl" ));
print ROUT "SEQUENCE_ID\tTAXON_ID\tSPECIES\tFAMILY\tSEQ_LENGTH\tNUCS_ALN\tCOVERAGE_PCT\tBHIT_IDENTITY\tBESTHSP_COUNT\n";
$c = "T"; $n = 0;
foreach my $j (keys %tids){  # 0.SEQLEN 1.NUCALN 2.COVPCT 3.BHIDENT 4.COUNTSblast 5.TXNID 6.SPC 7.FAM
    my (@k,$kn);
#    $kn = scalar @{ $tids{$j} } - 1;
    ( @{$tids{$j}}[2] >= $mincov and @{$tids{$j}}[4] > 0) || next;
   # if (@{$tids{$j}}[4] > 0 ) { 
        print ROUT join("\t", $j, ( map { defined($tids{$j}[$_]) ?  $tids{$j}[$_] : "NA" } (5,6,7,0,1,2,3,4) )), "\n";
        my $taxon=@{$tids{$j}}[5];
        my $coverage=@{$tids{$j}}[2];
        my $pidentity=@{$tids{$j}}[3];
        my $nconts=@{$tids{$j}}[4];
        if (exists $specs{$taxon}) {
            #Update coverage info:
            if ( $coverage < $specs{$taxon}[2] ) { $specs{$taxon}[2]=$coverage; }#cov min
            if ( $coverage > $specs{$taxon}[3] ) { $specs{$taxon}[3]=$coverage; }#cov max
            $specs{$taxon}[4] += $coverage; # add cov to coverage summatory
            #Update identity perc info:
            if ( $pidentity < $specs{$taxon}[5] ) { $specs{$taxon}[5]=$pidentity; } #pid min
            if ( $pidentity > $specs{$taxon}[6] ) { $specs{$taxon}[6]=$pidentity; }#pid max
            $specs{$taxon}[7] += $pidentity; # add pid to Pident sumatory
            #Update counters:
            $specs{$taxon}[8] += 1; # add +1 to refseqs sumatory
            $specs{$taxon}[9] += $nconts;  # add n contigs to NContigs summatory
            # Update info field:
            my $i="SEQ:$j,LEN:$tids{$j}[0],NUCALN:$tids{$j}[1],COV:$coverage,BHIDENT:$pidentity,N:$nconts";
            my $info= "$specs{$taxon}[10];$i";
            $specs{$taxon}[10]=$info;  # Update info field.
        } else {
            my $info="SEQ:$j,LEN:$tids{$j}[0],NUCALN:$tids{$j}[1],COV:$coverage,BHIDENT:$pidentity,N:$nconts";
            $specs{$taxon} = [ $tids{$j}[6], $tids{$j}[7], $coverage, $coverage, $coverage, $pidentity, $pidentity, $pidentity, 1, $nconts, $info ];
        }
    # }
    #print ROUT join("\t", $j,  $tids{$j}[5], $tids{$j}[6], $tids{$j}[7], $tids{$j}[0], $tids{$j}[1], $tids{$j}[2], $tids{$j}[3], $tids{$j}[4] ), "\n";
} continue {
    $n++;
    &counter($n,$c) if ($n % 1000 == 0);
}; # while GIN
&counter_end($n);
close(ROUT);


print STDERR "### WRITTING SPECIES INFORMATION\n";
open ( SOUT, ">", join( "_", $out_prefix, "taxonomysum_byspecie.tbl" ));
print SOUT "TAXON_ID\tSPECIES\tFAMILY\tMIN_COV\tMAX_COV\tMEAN_COV\tMIN_PID\tMAX_PID\tMEAN_PID\tN_SEQS\tN_CONTIGS\tINFO\n";
foreach my $sp (keys %specs){ 
    print SOUT "$sp\t".join( "\t", $specs{$sp}[0], $specs{$sp}[1], 
                                    $specs{$sp}[2], $specs{$sp}[3],
                                    $specs{$sp}[4]/$specs{$sp}[8],
                                    $specs{$sp}[5], $specs{$sp}[6],
                                    $specs{$sp}[7]/$specs{$sp}[8],
                                    $specs{$sp}[8], $specs{$sp}[9],
                                    $specs{$sp}[10])."\n";
};
close(SOUT);
exit(0);
