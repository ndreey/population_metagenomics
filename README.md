# population_metagenomics
Population metagenomics of a symbiont across its ecologically divergent Tephritid fly host ecotypes


# Required up-stream metagenomic analysis
This code depends on the outputs from the up-stream metagenomics analysis.
https://github.com/ndreey/stammerula2025

# Down-stream analysis of _Stammerula_

## 1. Copy over the files and summarize
As the stammerula2025 pipeline `sys links` the files to save memory space, we need to copy them over to have a hard copy.

`sbatch scripts/resource/copy-nf-stam-results.sh`

We can then summarize the stats!

Number of bins per population
`bash scripts/resource/collect-raw-bin-stats.sh`

Metagenome assembly statistics (requires bbmap)
`bash scripts/resource/metagenome-assembly-stats.sh`

BUSCO
`bash scripts/stats/busco-stats.sh`

BAKTA
`bash scripts/stats/bakta-stats.sh`

CHECKM2
`bash scripts/stats/checkm2-stats.sh`

GTDB-Tk
`bash scripts/stats/gtdbtk-stats.sh`

## 2. Find STAMMERULA

Create BLAST database with all bins and SILVA database
```bash
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz

sbatch scripts/create-blast-db.sh
sbatch scripts/runBLAST.sh
```

## 2. Lets build custom Kraken2 Database

```bash
sbatch scripts/build-k2-db.sh
sbatch scripts/runKraken2.sh

cat 09-KRAKEN2/decon/kraken2_output/*.output | grep "Candidatus Stammerula" | grep "^C" | wc -l
4,030,035

cat 09-KRAKEN2/raw/kraken2_output/*.output | grep "Candidatus Stammerula" | grep "^C" | wc -l
6,105,747


```

## 3. Prepare metagenome and mag
This will create the required index files for downstream tools.

```bash
mkdir -p metagenome mag

gunzip -c 05-metagenomes/01-metamdbg/CHST-pt_042.contigs.fa.gz > metagenome/CHST-pt_042.contigs.fa
cp 06-metaWRAP-refined-bins/CHST-pt_042/bin_refinement/metawrap_50_10_bins/bin.3.fa mag/ctg4994.fa

sbatch scripts/resource/prep-genome.sh metagenome/CHST-pt_042.contigs.fa
```

## map2meta-gvcf
This will map to the metagenome and create GVCFs for each sample.

```bash
sbatch scripts/map2meta-gvcf.sh
```

## makeGATKvcfs.sh
1. Combines the GVCFs
2. Genotypes them to create a final metagenome VCF
3. Subsets the VCF to only have records for MAG

```bash
sbatch scripts/makeGATKvcfs.sh \
	--gvcf-list doc/gvcfs_paths.list \
	--ref metagenome/CHST-pt_042.contigs.fa \
	--mag-headers doc/stamMAGheader.txt \
	--outdir 12-GATK-VCFs \
	--prefix stam
```

## wrangleVCF
This will wrangle the MAG vcf file. Meaning it will create:
- filterd biallelic SNP vcf file
- filtered biallelic SNP vcf file + invariants
- All intermediary files
- Creates summarized stats and validate the vcf files.
  
```bash
sbatch scripts/wrangleVCF.sh \
	--vcf 12-GATK-VCFs/stam/vcf-final/stam.mag.raw.vcf.gz \
	--ref mag/ctg4994.fa \
	--outdir 13-STAM_VCF \
	--prefix stam
```


## MAG COVERAGE

```bash
# Convert to BED
mamba activate bedops
gff2bed < 07-bin-quality-assessment/CHST-pt_042/bakta/bakta_output/bin.3/bin.3.txt.gff3 > 07-bin-quality-assessment/CHST-pt_042/bakta/bakta_output/bin.3/bin.3.bed

# Run qualimap
sbatch scripts/runQualiMap.sh

# Parse the results
bash scripts/parseQualimap.sh

```

## MAG KEGG

```bash
# DRAM KO
cat 15-DRAM/annotations.tsv | cut -f 9 | grep "K" | paste -sd, | > ko_DRAM.list

# mettannotator
cat 16-mettannotator/ctg4994_annotations.gff | grep -o 'kegg=ko:K[0-9]\+' | cut -f 2 -d ":" | paste -sd, > ko_mettannotator.list

# union
cat 15-DRAM/annotations.tsv | cut -f 9 | grep "K" > tmp.ko_DRAM.list
cat 16-mettannotator/ctg4994_annotations.gff | grep -o 'kegg=ko:K[0-9]\+' | cut -f 2 -d ":" > tmp.ko_mettannotator.list
cat bakta_output/user_ko.txt | cut -f 2 | grep -v "^PK" > tmp.ko_bakta.list
cat tmp.ko_* | sort | uniq | paste -sd, > ko_union.list

# Get completeness
give_completeness --input-list ko_DRAM.list -o 15-DRAM/
give_completeness --input-list ko_mettannotator.list -o 16-mettannotator/
give_completeness --input-list bakta_output/bakta_ko.list -o bakta_output/
give_completeness --input-list ko_union.list -o KEGG_UNION

# Summary
 cat bakta_output/summary.kegg_pathways.bakta.tsv | grep -v "module_" | wc -l
153

cat 15-DRAM/summary.kegg_pathways.tsv | grep -v "module_" | wc -l
149

cat 16-mettannotator/summary.kegg_pathways.tsv | grep -v "module_" | wc -l
155

cat KEGG_UNION/summary.kegg_pathways.tsv | grep -v "module_" | wc -l
159

# Make KEGG MAPPER FRIENDLY FILE
cat ko_union.list | tr "," "\n" > KEGG_MAPPER_FRIENDLY.tmp
awk '{print "gene" NR "\t" $0}' KEGG_MAPPER_FRIENDLY.tmp > KEGG_MAPPER_FRIENDLY.tsv

# Pate raw text into text file, grep only the "^M0", then separate so we have one column module, other pathway name.
cat KEGG_UNION/summary.kegg_pathways.tsv | grep -v "module_" | cut -f 1 > keg_pathways.summary.modules.txt
cat KEGG_MAPPER_UNION_RESULTS.raw.txt | grep -vf keg_pathways.summary.modules.txt > MISSING_MODULES.raw.txt

# Manually add these to the summary.results.




# Manually removed all animal, fungi, or plant specific pathways
cat KEGG_UNION/summary.kegg_pathways.tsv | grep -vf non_bact_modules.txt > KEGG_UNION/summary.kegg_pathways.bact.union.tsv
cat 15-DRAM/summary.kegg_pathways.tsv | grep -vf non_bact_modules.txt > 15-DRAM/summary.kegg_pathways.bact.dram.tsv
cat 16-mettannotator/summary.kegg_pathways.tsv | grep -vf non_bact_modules.txt > 16-mettannotator/summary.kegg_pathways
.bact.mettannotator.tsv

# Summary
cat 16-mettannotator/summary.kegg_pathways.bact.mettannotator.tsv | grep -v "module_" | wc -l
145

cat 15-DRAM/summary.kegg_pathways.bact.dram.tsv | grep -v "module_" | wc -l
139

cat KEGG_UNION/summary.kegg_pathways.bact.union.tsv | grep -v "module_" | wc -l
149
```


## Population structure
```bash
for K in {1..24}; do admixture --cv -j12 data/stam.passed.no-sc.no-outlier.bed $K > stam.passed.no-sc.no-outlier.log-$K.out ; done
```


```bash
vcf2phylip.py --input vcfs/stam.filt.bial.snps.passed.no-sc.no-outlier.vcf.gz --output-folder vcf2phylip/ --output-prefix stam.filt.bial.passed.no-sc.no-outlier --fasta --nexus

genovi -i bakta_output/bin.3/bin.3.gbff -s complete --size -cs purple -k -t "Stammerula tephritidis" -te --title_position center -o stam_purple

iqtree2 -s ROARY/core_gene_alignment.aln --mem 5GB --threads 4 --boot 500 -m MFP --prefix stam_roary
```
## Pyani
```
pyani-plus anim mags/ --database stam_mag_comp.db --create-db --name "stam ANIm"

pyani-plus anib mags/ --database stam_mag_comp.db --name "stam ANIb"

pyani-plus dnadiff mags/ --database stam_mag_comp.db --name "stam dnadiff"

pyani-plus fastani mags/ --database stam_mag_comp.db --name "stam fastani"

pyani-plus sourmash mags/ --database stam_mag_comp.db --name "stam sourmash"


pyani-plus list-runs --database stam_mag_comp.db
                        5 analysis runs in stam_mag_comp.db
┏━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━┳━━━━━━┳━━━━━━┳━━━━━━━┳━━━━━━━━┳━━━━━━━━━━━━━━━┓
┃ ID ┃ Date       ┃ Method   ┃ Done ┃ Null ┃ Miss ┃ Total ┃ Status ┃ Name          ┃
┡━━━━╇━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━╇━━━━━━╇━━━━━━╇━━━━━━━╇━━━━━━━━╇━━━━━━━━━━━━━━━┩
│  1 │ 2025-08-22 │ ANIm     │   81 │    0 │    0 │ 81=9² │ Done   │ stam ANIm     │
│  2 │ 2025-08-22 │ ANIb     │   81 │    0 │    0 │ 81=9² │ Done   │ stam ANIb     │
│  3 │ 2025-08-22 │ dnadiff  │   81 │    0 │    0 │ 81=9² │ Done   │ stam dnadiff  │
│  4 │ 2025-08-22 │ fastANI  │   81 │    0 │    0 │ 81=9² │ Done   │ stam fastani  │
│  5 │ 2025-08-22 │ sourmash │   81 │    0 │    0 │ 81=9² │ Done   │ stam sourmash │
└────┴────────────┴──────────┴──────┴──────┴──────┴───────┴────────┴───────────────┘
```


I have booked Stora Konferensrummet in Ekologihuset for 3 pm on September 3rd.
