// zegrep '^>' refseqs/rvdb_nt/C-RVDBvCurrent.fasta.gz | gawk -vFS="|" '{print $3}' > refseqs/rvdb_nt/C-RVDBvCurrent.ids
//efetch -db nuccore -input refseqs/rvdb_nt/C-RVDBvCurrent.ids -format gb > refseqs/rvdb_nt/C-RVDBvCurrent.gb

// gawk -vOFS="\t" '$1 ~ /^VERSION$/ { ID=$2 } match($0, /^ *\/db_xref="taxon:(.*)" *$/, a) { print ID, a[1] }' refseqs/rvdb_nt/C-RVDBvCurrent.gb > refseqs/rvdb_nt/C-RVDB_genebank_to_taxon.ids


process get_taxonids () {
    input:
        val odir
        val fasta_fl
        val name
        val accession2taxid

    output:
        val txid_fl

    script:
        idsfl="${odir}/${name}_genomes.ids"
        gbfl="${odir}/${name}_missingtax_genomes.gb"
        missids="${odir}/${name}_missingtax.txt"
        txid_fl="${odir}/${name}_genebank_to_taxon.gb"

        """
        zegrep '^>' $fasta_fl | awk '{print \$1}' - | sed 's/>//' > $idsfl; 
        zcat $accession2taxid   | \ 
            gawk 'BEGIN{
                while(getline<ARGV[1]>0) idslist[\$1]; ARGV[1]=""
            } (\$2 in idslist) { 
                 print \$2, \$3 > "$txid_fl"; delete idslist[\$2] 
            } END{ 
                 for(i in idslist) {print i > "$missids" }
            }' $idsfl - ;

        efetch -db nuccore -input $missids -format gb > $gbfl;
        gawk -vOFS="\\t" '
            \$1 ~ /^VERSION\$/ { 
                ID=\$2 
            } match(\$0, /^ *\\/db_xref="taxon:(.*)" *\$/, a) { 
                print ID, a[1] 
            }' $gbfl >> $txid_fl;
        """
}

process get_taxonids_rvdb () {
    input:
        val odir
        val fasta_fl
        val name      // full path
        val accession2taxid

    output:
        val txid_fl

    script:
        idsfl=name.replaceAll(".fasta.gz","_genomes.ids")
        gbfl=name.replaceAll(".fasta.gz","_missingtax_genomes.gb")
        missids=name.replaceAll(".fasta.gz", "_missingtax.txt")
        txid_fl=name.replaceAll(".fasta.gz","_genebank_to_taxon.gb")

        """
        zegrep '^>' $fasta_fl | awk -vFS="|" '{print \$3}' -  > $idsfl; 
        zcat $accession2taxid   |  gawk 'BEGIN{
                                        while(getline<ARGV[1]>0) idslist[\$1]; ARGV[1]=""
                                    } (\$2 in idslist) { 
                                        print \$2, \$3 > "$txid_fl"; delete idslist[\$2] 
                                    } END{ 
                                        for(i in idslist) {print i > "$missids" }
                                    }' $idsfl - ;

        efetch -db nuccore -input $missids -format gb > $gbfl;
        gawk -vOFS="\\t" '
            \$1 ~ /^VERSION\$/ { 
                ID=\$2 
            } match(\$0, /^ *\\/db_xref="taxon:(.*)" *\$/, a) { 
                print ID, a[1] 
            }' $gbfl >> $txid_fl;
        """
}

process taxonomizator () {
    input:
        val odir
        val txidfl
        val bindir
        val n_dir

    output:
        val taxfl  // Taxonomy file for set

    script:
            taxfl=txidfl.replaceAll("_genebank_to_taxon.gb", ".tax.gz")
            """
            $bindir/taxonomizator_byid.pl          \
                $n_dir/merged.dmp               \
                $n_dir/names.dmp                \
                $n_dir/nodes.dmp                \
                $txidfl                          |\
                gzip -9 -  > $taxfl;
            """

}


