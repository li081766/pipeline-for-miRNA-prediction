perl get_uniq_seq.pl Genome_with_miRNAs.fa miRNA.mature.fa.blastn.filter.m8 >miRNA.mature.fa.blastn.filter.m8.nt110.fa

perl calNgap.pl miRNA.mature.fa.blastn.filter.m8.nt110.fa --whole | awk '$6>5' >miRNA.mature.fa.blastn.filter.m8.nt110.Nsize.gt5%.id

perl select_seq.pl miRNA.mature.fa.blastn.filter.m8.nt110.Nsize.gt5%.id miRNA.mature.fa.blastn.filter.m8.nt110.fa --except >miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.fa

perl predict_hairpin.pl miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.fa

perl with_gene_info.pl miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.xls genome.gene.filter.gff >miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls

perl select_seq.pl miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.xls --except | awk '{print $0"\tintergenic"}' >>miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls

less mmiRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls | grep -v 'exon' | sort -k2,2 -k3n >miRNA.final.xls
