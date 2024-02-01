# miniSNV

## Introduction
miniSNV is a lightweight SNV calling algorithm that simultaneously achieves high performance and yield. miniSNV takes the known common variants in the population as variation backgrounds and leverages read pileup, read-based phasing, and consensus generation to discover and genotype SNVs for ONT long reads. Benchmarks on real ONT datasets under various error profiles demonstrate that miniSNV has superior sensitivity and comparable accuracy on SNV detection and runs faster with outstanding scaling performance and lower memory than most state-of-the-art variant callers.

---
## Installation
**1.download index folder at** [here](https://drive.google.com/drive/folders/17NFTbnPuZhJ4SWSjrRqyUZ_9pn38vgc9?usp=sharing)  
  
**2.install whatshap :**
```
conda install bioconda::whatshap

or

conda install bioconda/label/cf201901::whatshap
```
**3.download miniSNV** 
```
git clone https://github.com/CuiMiao-HIT/miniSNV.git

cd miniSNV/Release

make -j 12

cd ..
```

---	
## Usage
```
python chr_chunk_task.py \
--fin_ref FIN_REF \
--fin_bam FIN_BAM \
--fin_index FIN_INDEX \
--fin_bed FIN_BED \
--dup_bed DUP_BED \
--workDir WORKDIR
```

### Options
**Required parameters:**  
```
	-b, --fin_bam           BAM file input. The input file must be samtools indexed.  
	-r, --fin_ref           FASTA reference file input. The input file must be samtools indexed.  
	-i, --fin_index         The folder path containing a miniSNV index(five files in the folder).  
	-hb, --homo_bed         HOMO Bed format input.  
	-db, --dup_bed          DUP HOMO Bed format input.  
	-o, --workDir           Work-directory for distributed job.  
```
**Other parameters:**  
```
	--threads(INT)          Number of threads to use.(default : 16)  
	--chrName               list of chrName (contig names) to be processed, separeted by comma without any blank space.  
	--sample                Sample name in vcf file.(default : SAMPLE)  
	--chunkWidth(INT)       Reference length to detect candidate in one loop.(default : 10000000)  
```
## Contact
Please post on [Github Issue](https://github.com/CuiMiao-HIT/miniSNV/issues) or contact cuimiao@stu.hit.edu.cn.
