# pipeline-for-miRNA-prediction

This pipeline is used for prediciton for miRNA in genome, invoking blastn, miRNAFold and RNAfold

```
/Software/blast/bin/makeblastdb -dbtype nucl -in genome.hardmasked.fa -parse_seqids

/Software/blast/bin/blastn -task blastn-short -db genome.hardmasked.fa -query miRNA.mature.fa -num_threads 4 -outfmt 6 -out miRNA.mature.fa.blastn.m8 -ungapped -penalty -1 -reward 1

perl filter_m8.pl miRNA.mature.fa.blastn.m8 miRNA.mature.fa >miRNA.mature.fa.blastn.filter.m8

perl get_uniq_seq.pl Genome_with_miRNAs.fa miRNA.mature.fa.blastn.filter.m8 >miRNA.mature.fa.blastn.filter.m8.nt110.fa

perl calNgap.pl miRNA.mature.fa.blastn.filter.m8.nt110.fa --whole | awk '$6>5' >miRNA.mature.fa.blastn.filter.m8.nt110.Nsize.gt5%.id

perl select_seq.pl miRNA.mature.fa.blastn.filter.m8.nt110.Nsize.gt5%.id miRNA.mature.fa.blastn.filter.m8.nt110.fa --except >miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.fa

perl predict_hairpin.pl miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.fa

perl with_gene_info.pl miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.xls genome.gene.filter.gff >miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls

perl select_seq.pl miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.xls --except | awk '{print $0"\tintergenic"}' >>miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls

less miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls | grep -v 'exon' | sort -k2,2 -k3n >miRNA.final.xls
