perl get_uniq_seq.pl Genome_with_miRNAs.fa migratory_locust.miRNA.mature.fa.blastn.filter.m8 >migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.f
a
perl /home/jiali/pl_script/seq_tool/calNgap.pl  migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.fa --whole | awk '$6>5' >migratory_locust.miRNA
.mature.fa.blastn.filter.m8.nt110.Nsize.gt5%.id
fishInWinterTf migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.Nsize.gt5%.id migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.fa --excep
t >migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.fa
perl predict_hairpin.pl migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.fa
perl with_gene_info.pl migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.xls SchistoV2.gene.filter.gff >migratory_lo
cust.miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls
fishInWintert migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls migratory_locust.miRNA.mature.fa.blast
n.filter.m8.nt110.filter_N.premature.best.info.xls --except | awk '{print $0"\tintergenic"}' >> migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110
.filter_N.premature.best.info.with_gene.xls
le migratory_locust.miRNA.mature.fa.blastn.filter.m8.nt110.filter_N.premature.best.info.with_gene.xls | grep -v 'exon' | sort -k2,2 -k3n >desert_locust
.miRNA.final.xls
