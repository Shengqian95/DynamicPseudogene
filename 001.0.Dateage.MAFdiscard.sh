#######################################
cd /opt/qians/Pseudogene/Data/UCSC.maf/Mouse
sed 's/^/chr/g' /home/qians/Pseudo/Data/Ref/Mouse/geneMouse.pseudogene.bed > geneMouse.pseudogene.chr.bed
mkdir Human Chimp Rhesus Rat Rabbit Dog Opossum Platypus Chicken XTropicalis Zebrafish
#For Human
cd Human
for i in {1..19} X Y; do wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsHg38/mafSynNet/chr${i}.maf.gz; done
gunzip *gz
#extract align region
for i in {1..19} X Y; do /home/qian/source/mafsInRegion /opt/qians/Pseudogene/Data/UCSC.maf/Mouse/geneMouse.pseudogene.chr.bed chr${i}.mafsInRegion.maf chr${i}.maf; done
rm chr*
#get bed file
cat chr*.mafsInRegion.maf > all.mafsInRegion.maf
rm chr*.mafsInRegion.maf

grep ^s all.mafsInRegion.maf  | sed 's/  */ /g' |sed 's/ /\t/g' | awk -F "\t" '{ if($2~"mm10") print $0"\t"NR+1; else print $0"\t"NR}' > all.mafsInRegion.align
grep mm10 all.mafsInRegion.align | sed 's/mm10.//g' | awk -F "\t" '{print $2"\t"$3"\t"$3+$4"\t"$4"\t"$8"\t"$5}' | sed 's/chr//g' > all.mafsInRegion.align.mm10.bed
#intersect with gene bed file
bedtools intersect -a all.mafsInRegion.align.mm10.bed -b ~/Pseudo/Data/Ref/Mouse/geneMouse.pseudogene.bed -wa -wb > all.mafsInRegion.align.mm10.pseudogene.bed
rm all.mafsInRegion.align all.mafsInRegion.align.mm10.bed

######################################
#For Mouse
###wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

for i in {1..19} X Y; do wget -P Chimp http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsPanTro6/mafSynNet/chr${i}.maf.gz; done
wget -P Rhesus/ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsRheMac10/mm10.rheMac10.synNet.maf.gz
for i in {1..19} X Y; do wget -P Rat http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsRn6/mafSynNet/chr1.maf.gz; done
for i in {1..19} X Y; do wget -P Rabbit http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsOryCun2/mafRBestNet/chr${i}.maf.gz; done
wget -P Dog http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsCanFam4/mm10.canFam4.synNet.maf.gz

cd XTropicalis
ftp hgdownload.cse.ucsc.edu; cd goldenPath/mm10/vsXenTro3/; prompt off; mget -a
#usage http://www.genome.ucsc.edu/goldenPath/help/ftp.html
#wget -P XTropicalis/ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsXenTro3/axtNet/chr1.mm10.xenTro3.net.axt.gz
#wget -P XTropicalis/ http://hgdownload.soe.ucsc.edu/goldenPath/xenTro3/bigZips/xenTro3.chrom.sizes
#usage example: https://groups.google.com/a/soe.ucsc.edu/g/genome/c/9U27kOQqmWs
#/home/qian/source/axtToMaf chr10.mm10.xenTro3.net.axt /opt/qians/Pseudogene/Data/UCSC.maf/Mouse/mm10.chrom.sizes xenTro3.chrom.sizes -tPrefix=mm10. -qPrefix=xenTro3. chr10.mm10.xenTro3.net.axt.maf 
#for i in `ls *net.axt`;do /home/qian/source/axtToMaf $i /opt/qians/Pseudogene/Data/UCSC.maf/Mouse/mm10.chrom.sizes xenTro3.chrom.sizes -tPrefix=mm10. -qPrefix=xenTro3. ${i}.maf; done
grep chr chr*axt >> all.mm10.xenTro3.net.axt
awk -F " " '{print $1"\t"$2"\t"$3"\t"$9"\t"$8"\t"$7}' all.mm10.xenTro3.net.axt.ID > Mouse.mm10.XTropicalis.xenTro3.net.bed
awk -F " " '{print $4"\t"$5"\t"$6"\t"$9"\t"$8"\t"$7}' all.mm10.xenTro3.net.axt.ID > XTropicalis.xenTro3.Mouse.mm10.net.bed



wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsHg38/mafSynNet/chr1.maf.gz
gunzip *gz
sed 's/^/chr/g' /home/qians/Pseudo/Data/Ref/Mouse/geneMouse.pseudogene.bed > geneMouse.pseudogene.chr.bed
/home/qian/source/mafsInRegion geneMouse.pseudogene.chr.bed chr1.out.maf chr1.maf

#https://www.thegeekstuff.com/2010/02/awk-conditional-statements/
grep ^s chr1.out.maf | sed 's/  / /g' |sed 's/ /\t/g' | awk -F "\t" '{ if($2~"hg38") print $0"\t"NR-1; else print $0"\t"NR}' > a.test
grep mm10 a.test | sed 's/mm10.//g' | awk -F "\t" '{print $2"\t"$3"\t"$3+$4"\t"$4"\t"$8"\t"$5}' | sed 's/chr//g' > a.bed



cd /opt/qians/Pseudogene/Data/UCSC.maf/Human
mkdir Chimp Rhesus Mouse Rat Rabbit Dog Opossum Platypus Chicken XTropicalis Zebrafish

nohup /home/qian/source/lastz-distrib-1.04.00/bin/lastz ~/Pseudo/Data/Ref/Human/Homo_sapiens.GRCh38.dna.primary_assembly.fa[multiple,subset=a]  ~/Pseudo/Data/Ref/Mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa[multiple,subset=b] --notransition --step=20 --nogapped --format=maf > Human98.Mouse98.maf &


R filter

wget ftp://ftp.ensembl.org/pub/release-101/gtf/danio_rerio//Danio_rerio.GRCz11.101.gtf.gz
zcat Danio_rerio.GRCz11.101.gtf.gz | grep -w gene | awk -F ";" '{print $1"\t"$(NF-1)"\t"$(NF-3)}'| awk -F "\t" '{print "chr"$1"\t"$4"\t"$5"\t"$9"\t"$10"\t"$11}' | sed 's/gene_id //g' | sed 's/"//g'| sed 's/gene_biotype //g' | sed 's/gene_name //g' > Zebrafish.gene.bed
for i in `grep TRUE Mouse.Zebrafish.pseudogene.ratio.txt| awk -F " " '{print $1}'`; do grep $i Mouse.Zebrafish.block.pseudogene.blockID.bed| awk -F "\t" '{print $11}' ;done >> overlapped.blockID.txt
for i in `cat overlapped.blockID.txt`;do grep -w $i Zebrafish.Mouse.block.bed ;done >> Zebrafish.overlapped.blockID.bed
