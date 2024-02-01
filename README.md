# miniSNV

## Introduction
miniSNV is a lightweight SNV calling algorithm that simultaneously achieves high performance and yield. miniSNV takes the known common variants in the population as variation backgrounds and leverages read pileup, read-based phasing, and consensus generation to discover and genotype SNVs for ONT long reads. Benchmarks on real ONT datasets under various error profiles demonstrate that miniSNV has superior sensitivity and comparable accuracy on SNV detection and runs faster with outstanding scaling performance and lower memory than most state-of-the-art variant callers.

---
## Installation
$ download index folder at [here](https://drive.google.com/drive/folders/17NFTbnPuZhJ4SWSjrRqyUZ_9pn38vgc9?usp=sharing)
$ install whatshap : **conda install bioconda::whatshap** or **conda install bioconda/label/cf201901::whatshap**
$ git clone https://github.com/CuiMiao-HIT/miniSNV.git && cd miniSNV/Release && make -j 12

---	

