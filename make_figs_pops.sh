fembot -r1 -t'Making figures.'

gs_ddRAD2015 ld_figures --ld-pickle /home/gus/remote_mounts/louise/data/genomes/glossina_fuscipes/annotations/SNPs/vcftools_out/ddrad58_populations/individuals/tsetseFINAL_14Oct2014_f2_53.recode.renamed_scaffolds.maf0_05.KG_indiv.geno.ld.pkl --out-dir /home/gus/gs_2015_ld/KG --contig-length /home/gus/Dropbox/uganda_data/data_repos/genome_info/assembly_info/contig_name_length.csv --formats png --save-tables

fembot -r1 -t'done with KG'

gs_ddRAD2015 ld_figures --ld-pickle /home/gus/remote_mounts/louise/data/genomes/glossina_fuscipes/annotations/SNPs/vcftools_out/ddrad58_populations/individuals/tsetseFINAL_14Oct2014_f2_53.recode.renamed_scaffolds.maf0_05.NB_indiv.geno.ld.pkl --out-dir /home/gus/gs_2015_ld/NB --contig-length /home/gus/Dropbox/uganda_data/data_repos/genome_info/assembly_info/contig_name_length.csv --formats png --save-tables

fembot -r1 -t'done with NB'

gs_ddRAD2015 ld_figures --ld-pickle /home/gus/remote_mounts/louise/data/genomes/glossina_fuscipes/annotations/SNPs/vcftools_out/ddrad58_populations/individuals/tsetseFINAL_14Oct2014_f2_53.recode.renamed_scaffolds.maf0_05.MS_indiv.geno.ld.pkl --out-dir /home/gus/gs_2015_ld/MS --contig-length /home/gus/Dropbox/uganda_data/data_repos/genome_info/assembly_info/contig_name_length.csv --formats png --save-tables

fembot -r1 -t'done with MS'

gs_ddRAD2015 ld_figures --ld-pickle /home/gus/remote_mounts/louise/data/genomes/glossina_fuscipes/annotations/SNPs/vcftools_out/ddrad58_populations/individuals/tsetseFINAL_14Oct2014_f2_53.recode.renamed_scaffolds.maf0_05.OT_indiv.geno.ld.pkl --out-dir /home/gus/gs_2015_ld/OT --contig-length /home/gus/Dropbox/uganda_data/data_repos/genome_info/assembly_info/contig_name_length.csv --formats png --save-tables

fembot -r1 -t'done with ALL'
