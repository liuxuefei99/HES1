

samtools='/data3/liuxuefei/software/conda21/envs/celescope1/bin/samtools'
sambamba='/data3/liuxuefei/software/conda21/bin/sambamba'
bowtie2='/data3/liuxuefei/software/bowtie2-2.4.5-linux-x86_64/bowtie2'
macs2='/data3/liuxuefei/software/conda21/bin/macs2'
bedtools='/data3/liuxuefei/software/conda21/bin/bedtools'
bed2bigwg='/data3/liuxuefei/software/bedGraphToBigWig'
fasta='/data3/dell/mamba/envs/homer/share/homer/data/genomes/hg38/genome.fa'

bowtie2.mm10.index='/data3/liuxuefei/database/human38/bowtie2/hg38'
chome.size='/data3/liuxuefei/database/human38/bowtie2/hg38.chrom.sizes.txt'
fai='/data3/liuxuefei/database/human38/all.fa/GRCh38.fa.fai'
gtf='/data3/liuxuefei/database/human38/gencode.v25.annotation.gtf'
bamCoverage='/data3/liuxuefei/software/conda21/envs/deeptools/bin/bamCoverage'
computeMatrix='/data3/liuxuefei/software/conda21/envs/deeptools/bin/computeMatrix'
bamCompare='/data3/liuxuefei/software/conda21/envs/deeptools/bin/bamCompare'
plotHeatmap='/data3/liuxuefei/software/conda21/envs/deeptools/bin/plotHeatmap'

findMotifsGenome='/data3/dell/mamba/envs/homer/bin/findMotifsGenome.pl'
annotatePeaks='/data3/dell/mamba/envs/homer/bin/annotatePeaks.pl'




files=list.files('/data2/liuxf/Mela_mouse_human_New/Cuttag_all/raw') 
files=files %>% str_remove('_raw_\\d.fq.gz') %>% unique()


readpath='/data2/liuxf/Mela_mouse_human_New/Cuttag_all/'
TNtable=data.frame(row.names = files[files %>% str_count('^HES1|^JUN|^KLF4') %>% as.logical()] ,
                   V1=files[files %>% str_count('^HES1|^JUN|^KLF4') %>% as.logical()],
                   V2= c('Ctrl_A375_Invirto' , 'Ctrl_MEL167_Invirto' ,'Ctrl_MEL167_Invirto' , 'Ctrl_MEL167_KidenyM','Ctrl_MEL167_KidenyM',
                         "Ctrl_MEL167_LungM","Ctrl_MEL167_LungM" ,'Ctrl_MEL167_OvaryM' ,'Ctrl_MEL167_OvaryM' ,"Ctrl_A375_Invirto" ,"Ctrl_MEL167_Invirto_ForJUN" ,"Ctrl_MEL167_Invirto_ForJUN",
                         'Ctrl_MEL167_KidenyM','Ctrl_MEL167_KidenyM', "Ctrl_MEL167_LungM","Ctrl_MEL167_LungM" ,"Ctrl_A375_Invirto_ForKLF4" ,"Ctrl_A375_Invirto_ForKLF4" , 'Ctrl_MEL167_Invirto',
                         'Ctrl_MEL167_Invirto' , 'Ctrl_MEL167_KidenyM','Ctrl_MEL167_KidenyM', "Ctrl_MEL167_LungM","Ctrl_MEL167_LungM" ) )




#--------1.clean fastq--------
for(i in files){
  sink(paste0(readpath,'/program/1.fastp.',i))
  cat(paste0('mkdir -p ',readpath,'1.clean/',i,'\n'))
  cat(paste0('mkdir -p ',readpath,'Run_log\n'))
  cat(paste0('/data3/liuxuefei/software/fastp -i ',readpath,'raw/',i,'_raw_1.fq.gz -I ',readpath,'raw/',i,'_raw_2.fq.gz -o ',readpath,'/1.clean/',i,'/',i,'_clean_R1.fq.gz -O ',readpath,'/1.clean/',i,'/',i,'_clean_R2.fq.gz -w 4'))
  sink()
}

sink(paste0(readpath,'/program/fastp.sh'))
for(i in files)cat(paste0('bash ',readpath,'/program/1.fastp.',i,'\n'))
sink()

#-----2.bowtie2比对--------
for(i in files){
  sink(paste0(readpath,'/program/2.bowtie',i))
  cat(paste0('mkdir -p ',readpath,'2.bowtie/',i,'\n'))
  cat(paste0('mkdir -p ',readpath,'Run_log\n'))
  cat(paste0('gzip -d ',readpath,'/1.clean/',i,'/',i,'_clean_R1.fq.gz\n'))
  cat(paste0('gzip -d ',readpath,'/1.clean/',i,'/',i,'_clean_R2.fq.gz\n'))
  cat(paste0(bowtie2,' -t -p 15 -x ',bowtie2.mm10.index,' -1 ',readpath,'/1.clean/',i,'/',i,'_clean_R1.fq -2 ',readpath,'/1.clean/',i,'/',i,'_clean_R2.fq | ',samtools,' sort -O bam -o ',readpath,'/2.bowtie/',i,'/',i,'_sort.bam\n'))
  cat(paste0(samtools,' index -b ',readpath,'/2.bowtie/',i,'/',i,'_sort.bam\n'))
  cat(paste0(samtools,' view -Sb -F 4 ',readpath,'/2.bowtie/',i,'/',i,'.sam -o ',readpath,'/2.bowtie/',i,'/',i,'.bam\n'))
  cat(paste0(samtools,' sort ',readpath,'/2.bowtie/',i,'/',i,'.bam -o ',readpath,'/2.bowtie/',i,'/',i,'_sort.bam\n'))
  cat(paste0(sambamba,' markdup -t 2 ',readpath,'/2.bowtie/',i,'/',i,'_sort.bam ',readpath,'/2.bowtie/',i,'/',i,'_sort_markdup.bam\n'))
  cat(paste0(samtools,' sort ',readpath,'/2.bowtie/',i,'/',i,'_sort_markdup.bam -o ',readpath,'/2.bowtie/',i,'/',i,'_sort_markdup_sort.bam\n'))
  cat(paste0(samtools,' index -b ',readpath,'/2.bowtie/',i,'/',i,'_sort_markdup_sort.bam\n'))
  
  #cat(paste0('rm ',readpath,'2.bowtie/',i,'/',i,'.sam ',readpath,'2.bowtie/',i,'/',i,'_sort.bam ',readpath,'2.bowtie/',i,'/',i,'.bam '))
  sink()
}

sink(paste0(readpath,'/program/bowtie2.sh'))
for(i in files)cat(paste0('bash ',readpath,'/program/2.bowtie',i,'\n'))
sink()



#-----3.macs2找peak-------
for(i in rownames(TNtable)){
  sink(paste0(readpath,'/program/3.macs2',i))
  cat(paste0('mkdir -p ',readpath,'3.macs2/',i,'\n'))
  cat(paste0('mkdir -p ',readpath,'Run_log\n'))
  cat(paste0(macs2,' callpeak -t ',readpath,'/2.bowtie/',i,'/',i,'_sort_markdup_sort.bam -c ',readpath,'/2.bowtie/',TNtable[i,2],'/',TNtable[i,2],'_sort_markdup_sort.bam -f BAM -g hs -n ',i,' -B -p 1e-4 --outdir ',readpath,'/3.macs2/',i,'/\n'))
  cat(paste0('perl /data3/liuxuefei/software/liruing/m6A_annotate_forGTF_xingyang_v2.pl ',gtf,' ',readpath,'/3.macs2/',i,'/',i,'_summits.bed ',readpath,'/3.macs2/',i,'/GTFpeak\n'))
  sink()
}

sink(paste0(readpath,'/program/macs2.sh'))
for(i in files)cat(paste0('bash ',readpath,'/program/3.macs2',i,'\n'))
sink()



#-------4.deeptools------
for(i in rownames(TNtable)){
  sink(paste0(readpath,'/program/4.deeptools',i))
  cat(paste0('mkdir -p ',readpath,'4.deeptools/',i,'\n'))
  cat(paste0('mkdir -p ',readpath,'4.deeptools/',TNtable[i,2],'\n'))
  cat(paste0('mkdir -p ',readpath,'Run_log\n'))
  cat(paste0(bamCoverage,' -b ',readpath,'/2.bowtie/',i,'/',i,'_sort_markdup_sort.bam -o ',readpath,'/4.deeptools/',i,'/',i,'bamcoverage.RPKM.bw --binSize 10 -p 6 --normalizeUsing RPKM \n'))
  cat(paste0(bamCoverage,' -b ',readpath,'/2.bowtie/',i,'/',i,'_sort_markdup_sort.bam -o ',readpath,'/4.deeptools/',i,'/',i,'bamcoverage_NoneNorma.bw --binSize 10 -p 6 --normalizeUsing None \n'))
  cat(paste0(bamCoverage,' -b ',readpath,'/2.bowtie/',TNtable[i,2],'/',TNtable[i,2],'_sort_markdup_sort.bam -o ',readpath,'/4.deeptools/',TNtable[i,2],'/',TNtable[i,2],'bamcoverage.RPKM.bw --binSize 10 -p 6 --normalizeUsing RPKM\n'))
  cat(paste0(bamCoverage,' -b ',readpath,'/2.bowtie/',TNtable[i,2],'/',TNtable[i,2],'_sort_markdup_sort.bam -o ',readpath,'/4.deeptools/',TNtable[i,2],'/',TNtable[i,2],'bamcoverage_NoneNorma.bw --binSize 10 -p 6 --normalizeUsing None\n'))
  cat(paste0(bamCompare,' -b1 ',readpath,'/2.bowtie/',i,'/',i,'_sort_markdup_sort.bam -b2 ',readpath,'/2.bowtie/',TNtable[i,2],'/',TNtable[i,2],'_sort_markdup_sort.bam -o ',readpath,'/4.deeptools/',i,'/',i,'bamComparelog2ratio.bw\n'))
  cat(paste0(computeMatrix,' reference-point --referencePoint center -p 6 -R ',readpath,'/3.macs2/',i,'/',i,'_summits.bed -S ',readpath,'/4.deeptools/',i,'/',i,'bamComparelog2ratio.bw -o ',readpath,'/4.deeptools/',i,'/',i,'callpeak.gz --outFileSortedRegions ',readpath,'/4.deeptools/',i,'/',i,'genes.bed -a 3000 -b 3000\n'))
  cat(paste0(plotHeatmap,' -m ',readpath,'/4.deeptools/',i,'/',i,'callpeak.gz --colorList \'white,#C1292E\' -out ',readpath,'/4.deeptools/',i,'/',i,'callpeak.pdf --heatmapHeight 15 --heatmapWidth 5  --xAxisLabel "" \n'))
  sink()
}




