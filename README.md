# Metagenomic-Analysis-Pipeline 🦍

Here is a metagenomic sequence data analysis pipeline, nothing different with other pipelines. But if you want to get to know the metagenomic analysis pipeline step by step, maybe you can get some details from this repository. And it is suitable for the beginners i think.

## Update[20241127]

- [pipe_metagenome_parallel.py](pipe_metagenome_parallel.py)

1. The concurrent running mode enhances efficiency.
2. It expands the selection of gene quantification methods by incorporating bwa, bowtie2, and salmon.
3. It introduces antibiotic resistance analysis.

## Pipeline overview 🐫

- Raw sequence quality trim
- Host reference sequence remove
- Metaphlan for composition of microbial communities
- Sequence assembly
- Gene prediction
- Remove redundancy gene and build non-redundant geneset
- Function annotation using emapper
- Organize function results table

## Quick start 🦏

```shell
$python pipe_metagenome.py -h
usage: 
=================================================================
python pipe_metagenome.py
	--fastq_list fq.list
	--output_dir result
	--ref ref_bowtie2_index
ref_bowtie2_index:
canis: /root/database/Canis_GCF_000002285.5/Canis_GCF_000002285_5
human: /root/database/hg38_GCF_000001405.40/GCF_000001405.40/hg38
=================================================================

Pipeline of metagenome

optional arguments:
  -h, --help            show this help message and exit
  -l FQLIST, --fastq_list FQLIST
                        raw fq list
  -o OUTDIR, --output_dir OUTDIR
                        result output
  -r REF, --ref REF     ref genome bowtie2 index
```

What you need to do is to provide two input files:

- --fastq_list # sample - fq_R1 - fq_R2 list, format like `fq.list`
- --ref # host reference genome bowtie2 index
  
And set output dir `--output_dir`, all of output results would be included. 

## Output files explanation 🐊

### Output files tree (not show all files)

```shell
├── 00-result                           # most important results in this fold
│   ├── 00_merged_abundance_table.txt   # composition of microbial communities
│   ├── 01_metaphlan_phylum.txt         # communities in phylum level
│   ├── 02_metaphlan_class.txt          # communities in class level
│   ├── 03_metaphlan_order.txt          # communities in order level
│   ├── 04_metaphlan_family.txt         # communities in family level
│   ├── 05_metaphlan_genus.txt          # communities in genus level
│   ├── 06_metaphlan_species.txt        # communities in species level
│   ├── KO_samples.xls                  # KEGG KO gene composition table
│   └── pathway_samples.xls             # KEGG pathway composition table
├── 01-fastp_trim
├── 02-ref_remove
├── 03-metaphlan
├── 04-megahit
├── 05-prodigal
├── 06-cdhit
├── 07-emapper
├── 08-sam_count
├── 09-emapper_kegg
```

### Output important result files explanation

#### metaphlan_diff-levels.txt

- column 1 : communities information
- column 2 ~ : communities abundance percent of samples

```shell
clade_name	C1
k__Bacteria;p__Firmicutes	86.51132
k__Bacteria;p__Actinobacteria	6.52203
k__Bacteria;p__Bacteroidetes	4.24729
k__Bacteria;p__Proteobacteria	1.96422
k__Bacteria;p__Fusobacteria	0.75514
k__Bacteria;p__Tenericutes	0.0
k__Bacteria;p__Spirochaetes	0.0
k__Bacteria;p__Verrucomicrobia	0.0
...
```

#### KO_samples.xls

- column 1 : KO gene name
- column 2 : KO gene description
- column 3 : KO gene id
- column 4 ~ : KO gene number of samples

```shell
KO_name	KO_des	KO	C1
"E1.1.1.1, adh"	alcohol dehydrogenase [EC:1.1.1.1]	K00001	2817.0
"AKR1A1, adh"	alcohol dehydrogenase (NADP+) [EC:1.1.1.2]	K00002	254.0
hom	homoserine dehydrogenase [EC:1.1.1.3]	K00003	2890.0
"BDH, butB"	"(R,R)-butanediol dehydrogenase / meso-butanediol dehydrogenase / diacetyl reductase [EC:1.1.1.4 1.1.1.- 1.1.1.303]"	K00004	20.0
...
```

#### pathway_samples.xls

- column 1 : pathway level 1
- column 2 : pathway level 2
- column 3 : pathway level 3
- column 4 : pathway id
- column 5 ~ : pathway gene number of samples

```shell
level1	level2	level3	pathway	C1
Metabolism	Carbohydrate metabolism	Glycolysis / Gluconeogenesis	ko00010	91005.0
Metabolism	Carbohydrate metabolism	Citrate cycle (TCA cycle)	ko00020	31442.0
Metabolism	Carbohydrate metabolism	Pentose phosphate pathway	ko00030	53905.0
Metabolism	Carbohydrate metabolism	Pentose and glucuronate interconversions	ko00040	21334.0
...
```

## Step by step 🦥

### Step1 Raw sequence quality trim using fastp

```shell
fastp -i sample_1.fastq.gz \
      -o sample_clean.1.fastq.gz \
      -I sample_2.fastq.gz \
      -O sample_clean.2.fastq.gz \
      -w 8 -h sample.html -j sample.json
```

### Step2 Host reference sequence remove

```shell
bowtie2 -x ref_bowtie2_index -1 sample_clean.1.fastq.gz -2 sample_clean.2.fastq.gz -S sample.sam  2>sample.mapping.log
samtools fastq -@ 8 -f 4 sample.sam -1 sample.unmap.1.fastq.gz -2 sample.unmap.2.fastq.gz -s sample.unmap.single.fastq.gz
```

### Step3 Metaphlan for composition of microbial communities

```shell
zcat sample.unmap.1.fastq.gz sample.unmap.2.fastq.gz|metaphlan --input_type fastq --bowtie2out sample_bowtie2.bz2 --output_file sample_metaphlan.tsv --nproc 8
# when you install metaphlan in the system, you will get script 'merge_metaphlan_tables.py', that's for merge different samples metaphlan result in one file, like: 
merge_metaphlan_tables.py *.tsv > 00_merged_abundance_table.txt
grep -E '(p__)|(clade_name)' 00_merged_abundance_table.txt |grep -v 'c__'|sed 's/|/;/g' > 01_metaphlan_phylum.txt
grep -E '(c__)|(clade_name)' 00_merged_abundance_table.txt |grep -v 'o__'|sed 's/|/;/g' > 02_metaphlan_class.txt
grep -E '(o__)|(clade_name)' 00_merged_abundance_table.txt |grep -v 'f__'|sed 's/|/;/g' > 03_metaphlan_order.txt
grep -E '(f__)|(clade_name)' 00_merged_abundance_table.txt |grep -v 'g__'|sed 's/|/;/g' > 04_metaphlan_family.txt
grep -E '(g__)|(clade_name)' 00_merged_abundance_table.txt |grep -v 's__'|sed 's/|/;/g' > 05_metaphlan_genus.txt
grep -E '(s__)|(clade_name)' 00_merged_abundance_table.txt |grep -v 't__'|sed 's/|/;/g' > 06_metaphlan_species.txt
```

### Step4 Sequence assembly and trim contigs which length < 500bp

```shell
megahit -1 sample.unmap.1.fastq.gz -2 sample.unmap.2.fastq.gz -o sample_megahit --out-prefix sample -t 8
seqkit seq -m 500 sample_megahit/sample.contigs.fa --remove-gaps > sample.contigs_500.fa
sed -i 's/>/>sample_/g' sample.contigs_500.fa
```

### Step5 Gene prediction using prodigal

```shell
prodigal -p meta -a sample_prot.faa -m -d sample_nucl.fna -o sample_genes.gff -f gff -s sample.stat -i sample.contigs_500.fa
```

### Step6 Remove redundancy gene and build non-redundant geneset

```shell
cat sample1_prot.faa sample2_prot.faa ... >  prot.faa
cat sample1_nucl.fna sample2_nucl.fna ... >  nucl.fna
cd-hit -i prot.faa -o prot_nonerude.faa -c 0.95 -T 8 -n 5 -d 0 -aS 0.9 -g 1 -sc 1 -sf 1 -M 0
grep '>' prot_nonerude.faa|awk -F ' ' '{print $1}'|sed 's/>//g' > prot_nonerude.list
seqtk subseq nucl.fna prot_nonerude.list > nucl_nonerude.fna
bwa index nucl_nonerude.fna -p geneset_bwa
bioawk -c fastx '{print $name, length($seq)}' nucl_nonerude.fna > geneset_length.txt
```

### Step7 Function annotation using emapper

```shell
# when you install emapper in the system, you will get script 'emapper.py'
emapper.py -i prot_nonerude.faa -o eggnog --cpu 0 --usemem
cut -f1,12 eggnog.emapper.annotations|grep -v "^#"|sed 's/ko://g'|sed '1i gene\tko'|grep -v "-" > KEGG_KO.txt
cut -f1,13 eggnog.emapper.annotations|grep -v "^#"|sed '1i gene\tpathway'|grep -v "-" > KEGG_PATHWAY.txt
```

### Step8 Gene num count

```shell
bwa mem -t 4 geneset_bwa sample.unmap.1.fastq.gz sample.unmap.2.fastq.gz | samtools view -bS - | samtools sort - > sample_mapping_geneset.bam
samtools view -F 4 -F 256 -F 2048 sample_mapping_geneset.bam|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\t"$2; else print $1"\t0"}'|sed '1i gene\tsample' > sample.count
```

### Step9 Organize function results table

```shell
# using /kegg/kegg.py to analysis, like:
$python kegg.py -h
usage: python kegg.py -kk KEGG_KO.txt -kp KEGG_PATHWAY.txt -mt merged_file.txt -ok out_KO.xls -op out_pathway.xls

Merge KO/pathway count table from eggnog result.

optional arguments:
  -h, --help            show this help message and exit
  -kk KEGGKO, --kegg_KO KEGGKO
                        Sample's kegg KO information, such as KEGG_KO.txt
  -kp KEGGPATHWAY, --kegg_pathway KEGGPATHWAY
                        Sample's kegg pathway information, such as
                        KEGG_PATHWAY.txt
  -mt MERGETABLE, --merge_table MERGETABLE
                        Sample's merged gene count table, such as
                        merged_file.txt
  -ok OUTKO, --out_KO OUTKO
                        Output KO result, such as out_KO.xls
  -op OUTPATHWAY, --out_pathway OUTPATHWAY
                        Output pathway result, such as out_pathway.xls
  -t TMP, --tmp TMP     Tmp files dir
```

## Contact 🐖

Wanjin Hu (<wanjin.hu@outlook.com>)

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=wanjinhu/Metagenomic-Analysis-Pipeline&type=Date)](https://www.star-history.com/#wanjinhu/Metagenomic-Analysis-Pipeline&Date)
