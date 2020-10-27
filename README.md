# pipeline-for-miRNA-prediction

This pipeline is used for prediciton for miRNA in genome, invoking blastn, miRNAFold and RNAfold

```
/Software/blast/bin/makeblastdb -dbtype nucl -in genome.hardmasked.fa -parse_seqids
/Software/blast/bin/blastn -task blastn-short -db genome.hardmasked.fa -query migratory_locust.miRNA.mature.fa -num_threads 4 -outfmt 6 -out miRNA.mature.fa.blastn.m8 -ungapped -penalty -1 -reward 1
perl ../filter_m8.pl migratory_locust.miRNA.mature.fa.blastn.m8 migratory_locust.miRNA.mature.fa >migratory_locust.miRNA.mature.fa.blastn.filter.m8
