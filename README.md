# Code for correlation-based gene-peak linking analysis
Integrating the activity of regulatory elements with gene expression is essential for understanding how the genome operates. Therefore, pairing genes that are expressed (as shown by RNA-seq) with peaks identified from ChIP-seq, ATAC-seq, etc., is important.<br>
Here, I have uploaded Python code that performs a correlation-based gene-peak linking analysis. The method of analysis is derived from [**MR Corces et al., Science 2018**](https://www.science.org/doi/10.1126/science.aav1898)<br>
<img src="https://github.com/keun-hong/gene-peak-linking/assets/43947916/a4b27b93-d2fd-44cd-afdd-b1cc9d898b13" width="500"><be>

## Overview of the analysis process
<img src="https://github.com/keun-hong/gene-peak-linking/assets/43947916/2ab19b00-afc9-4f5d-8d22-9bac87c867f1" width="700"><br>

## Required packages or tools
```bash
# 1. Python packages
pip install numpy pandas scipy rpy2

# 2. Bedtools
sudo apt-get install bedtools
```
## Input files
1. Expression table
2. Peak and signal table
Please check the test input files

## Usage
```tex
usage: gene-peak-linking.py <Exp> <Signal> <Type> <Range> <out_prefix> [options]
	Exp:    expression table (x-axis: gene or transcript IDs / y-axis: sample IDs)
	Signal: signal table (x-axis: peak IDs / y-axis: label, chr, start, end and sample IDs)
	Type:   gene or transcript
	Range:  range limit from tss which connect gene and peak (ex, 5000)
	out_prefix: prefix of output file name

Options:
  -h, --help            show this help message and exit
  -p LENGTH             promoter length ex) up,down (default = 1000,100)
  -z ZERO_FILTER        Cutoff ratio of degree containing non-zero in row
                        lines (default = 0.2)
  -v VARIANCE           Cutoff ratio of high variance range 0 - 1 (default =
                        0.75)
  --ae=ADD_RNA_VALUE    Add specific expression value (default = 1)
  --as=ADD_ATAC_VALUE   Add specific signal value (default = 0)
  --null=NULL_MODEL_FILE
                        Null model file (default = False)
  -n NUM_GENES          # of peaks used in null model (default = 500)
  -t THREADS            # of threads (default = 1)
```

