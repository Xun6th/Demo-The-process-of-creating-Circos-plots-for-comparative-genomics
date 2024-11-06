### Raw Data Input
## Get Reference genome/gff
cp ~/insect/cg/00.material/2_Trichogramma_dendrolimi/TD.fa ./
cp ~/insect/cg/00.material/2_Trichogramma_dendrolimi/TD.transcripts.gene.gff3 ./
## Get Query genome/gff
cp ~/insect/cg/00.material/3_Trichogramma_ostrinae/TO.fa ./
cp ~/insect/cg/00.material/3_Trichogramma_ostrinae/TO.transcripts.gene.gff3 ./

### Data Preprocess
## Remove scaffold (by seqkit)
seqkit seq -g -m 9999999 TD.fa > Td_chr.fa
seqkit seq -g -m 9999999 TO.fa > To_chr.fa
## Naming chromosomes according to rules
awk '{
    if ($1 ~ /^>/) {
        count++
        if (count <= 99) {
            print ">Td_Chr" count
        } else {
            print $0
        }
    } else {
        print $0
    }
}' Td_chr.fa > 1_Td_chr.fa
awk '{
    if ($1 ~ /^>/) {
        count++
        if (count <= 99) {
            print ">To_Chr" count
        } else {
            print $0
        }
    } else {
        print $0
    }
}' To_chr.fa > 1_To_chr.fa
## Synchronize this chromosomes/naming with the GFF file
sed 's/Chr01/Td_Chr1/g; s/Chr02/Td_Chr2/g; s/Chr03/Td_Chr3/g; s/Chr04/Td_Chr4/g; s/Chr05/Td_Chr5/g' TD.transcripts.gene.gff3 | grep -E 'Td_Chr1|Td_Chr2|Td_Chr3|Td_Chr4|Td_Chr5' > ./1_Td_chr.gff3
sed 's/Chr01/To_Chr1/g; s/Chr02/To_Chr2/g; s/Chr03/To_Chr3/g; s/Chr04/To_Chr4/g; s/Chr05/To_Chr5/g' TO.transcripts.gene.gff3 | grep -E 'To_Chr1|To_Chr2|To_Chr3|To_Chr4|To_Chr5' > ./1_To_chr.gff3

### Genome Alignment (by MUMmer;SyRI)
## Using SyRI to identify genomic rearrangements from whole-genome alignments 
## generated using MUMmer. A .tsv (out.filtered.coords) file is used as the input.
nucmer --mum --maxgap=500 --mincluster=100 ./1_Td_chr.fa ./1_To_chr.fa -prefix Qry_to
delta-filter -m -i 90 -l 100 Qry_to.delta > to_f.delta
show-coords -THrd to_f.delta > to_f.coords
syri -k -c to_f.coords -d to_f.delta -r ./1_Td_chr.fa -q ./1_To_chr.fa

### Tracks:Obtaining SNP\PAV\SV from output

## SNP
# By SyRI
grep "SNP" syri.out > to_syri.snp
awk 'BEGIN{OFS="\t"}{print $6,$7,$8,"1","+"}' to_syri.snp > to_syri_forTBtools.snp
# By MUMmer
delta-filter -1 -q -r Qry_to.delta > to_1qr.delta
show-snps -ClrTH to_1qr.delta > to_mummer.showsnps
awk '$2 == "." || $3 == "." { print > "to.pav" } $2 != "." && $3 != "." { print > "to.snp" }' to_mummer.showsnps
awk 'BEGIN{OFS="\t"}{print $12,$4,$4,"1","+"}' to.snp > to_mummer_forTBtools.snp
# These two types of outputs can be integrated
cat to_syri_forTBtools.snp to_mummer_forTBtools.snp > raw_to_all.snp
awk 'NR==1 || !seen[$0]++' raw_to_all.snp > to_all.snp

## PAV
# By MUMmer
awk '$3 != "." { print > "to_pav.raw" } $3 == "." { print > "td_pav.raw" }' to.pav
awk '{print > $12 "_part_" FILENAME}' to_pav.raw

awk -F'\t' '{print $4}' To_Chr1_part_to_pav.raw > 1
sort -n 1 > 1_sort
python ~/hang.py 1_sort 1_sort_merge
awk -F'\t' 'NF >= 2' 1_sort_merge > 1_sort_merge_filter
awk '{print "To_Chr1" "\t" $0}' 1_sort_merge_filter > 1_chr
cat 1_chr 2_chr 3_chr 4_chr 5_chr > td_pav.bed

## SV
awk '{print > $1 "_part_" FILENAME}' sv.txt
rm \#_part_sv.txt
cat *part* > all.sv
awk -F'\t' '
{
  # 提取第七列作为第一列，第四列和第五列作为第二列和第三列
  new_line = $7 "\t" $4 "\t" $5;
  
  # 检查第二项是否小于第三项，如果不是，则交换
  if ($4 > $5) {
    tmp = $4;
    $4 = $5;
    $5 = tmp;
    new_line = $7 "\t" $4 "\t" $5;
  }
  
  # 检查第二项与第三项的差是否大于99且小于99999
  if (($5 - $4 > 99) && ($5 - $4 < 99999)) {
    print new_line;
  }
}' all.sv > sv.bed

## Fraction of PAVs/SVs per _bp window (by bedtools)
bedtools makewindows -g chr.len -w _(bp) > windows.bed
bedtools sort -i PAVs/SVs.bed > sort.bed
bedtools merge -i sort.bed > merge.bed
bedtools coverage -a windows.bed -b merge.bed > output.cov


### Other Tracks
## Repeat Sequence (by RepeatMasker)
## TEs (by EDTA)
## Genes/GCratio/GCskew/Nratio (by TBtools)

### SyntenyLink (by TBtools:one step MCScanX)


