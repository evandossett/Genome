class SNP {
    public required string id;
    public required string id_pref;
    public required string chromosome;
    public required string position;
    public required char[]? genotype;
}

//1 snp_id(e.g., rs12345)                                                   <----------------------------
//2 snp_class(e.g., snp, in -del, etc.)
//3 allele(e.g., C / T)                                                     <----------------------------
//4 global_maf(the global minor allele frequency)
//5 global_maf_allele(which allele is the minor)
//6 chromosome(e.g., 7)                                                     <----------------------------
//7 chr_pos(GRCh37 or GRCh38 position, depending on your query parameters)
//8 cytogenetic(cytogenetic band)
//9 gene_symbol(closest or overlapping gene)
//10 gene_id(NCBI GeneID)
//11 fxn_class(functional class, e.g., missense_variant)
//12 weight                                                                 <----------------------------
//13 validation (e.g., by1000genomes, byHapMap, etc., or empty)
//14 avg_heterozygosity
//15 heterozygosity_se (standard error of heterozygosity)
//16 clinical_significance (e.g., Pathogenic, Benign)                       <----------------------------
//17 phenotype (if available)

class db_Chromosome {
    public List<db_SNP> db_SNPs { get; }
    public string? id { get; set; }

    public db_Chromosome(List<db_SNP> data) => db_SNPs = data;

    public db_Chromosome() => db_SNPs = new List<db_SNP>();
}

class db_SNP {
    public required string ID { get; set; }
    public required string Chromosome { get; set; }
    public required string Position { get; set; }
    public required string SNP_Class { get; set; }
    public List<Tuple<string, string>>? alleles { get; set; }
    public List<Tuple<string, string>>? MAFs { get; set; }
    public List<Tuple<string, string>>? Genes { get; set; }
    public string? Population { get; set; }
    public string? Clinical_Significance { get; set; }

}

//ID = <SNP_ID>...</SNP_ID>
//Chromosome = <CHR>...</CHR>
//Position = <CHRPOS_PREV_ASSM>...</CHRPOS_PREV_ASSM>
//SNP_Class = <SNP_CLASS>...</SNP_CLASS>
//MAFs = <GLOBAL_MAFS>...</GLOBAL_MAFS> THEN <MAF>...</MAF> THEN <STUDY>...</STUDY> AND <FREQ>...</FREQ>
//Genes = <GENES>...</GENES> THEN <GENE_E>...</GENE_E> THEN <NAME>...</NAME> AND <GENE_ID>...</GENE_ID>
//Population = <GLOBAL_POPULATION>...</GLOBAL_POPULATION>
//Clinical_Significance = <CLINICAL_SIGNIFICANCE>...</CLINICAL_SIGNIFICANCE>