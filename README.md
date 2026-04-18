# SSfSBT
Small Scripts for Small Bioinformatics Tasks.

Collection of python scripts as a backup for me and hopefully a help for you!

For details on usage and examples, please refer to the [wiki](https://github.com/SimonHegele/SSfSBT/wiki).

## Installation

```
conda create -n ssfsbt python>=3.10 # (optional but recommended)
conda activate ssfsbt               # (optional but recommended)

git clone https://github.com/SimonHegele/SSfSBT
cd SSfSBT
pip install .
```

## Contents

### Commandline-Tools

Tool                    | Summary
------------------------|--------
aln_pos2pos             | Finding corresponding positions in multiple sequence alignments
busco_find              | Extracting BUSCO transcript sequences
busco_merge             | Merging multiple reports from BUSCO and compiling a plot
kallisto2nanosim        | Converting expression profiles from [Kallisto](https://github.com/pachterlab/kallisto) for [NanoSim](https://github.com/bcgsc/NanoSim)
lengths                 | Basic sequence length distribution analysis for one or more FASTA/FASTQ-files
lr_lordec_contam_filter | Filtering contamination from long reads corrected by [HALC](https://github.com/lanl001/halc) or LoRDEC with [Kraken2](https://github.com/DerrickWood/kraken2) filtered short reads
fa2fq                   | FASTA to FASTQ conversion with different error rates for upper and lower case denoted nucleotides
fq2fa                   | FASTQ to FASTA conversion
gfa2fa                  | GFA to FASTA conversion
plot_msa                | Plotting multiple sequence alignments
rnaQUASTcompare         | Merging multiple reports from rnaQUAST and compiling a plot
sample                  | Subsampling a fixed number of sequences from FASTA/FASTQ-files
unambiguous_codes       | Replacing ambiguity codes in FASTA/FASTQ-files

Check the [wiki](https://github.com/SimonHegele/SSfSBT/wiki) for detailed description of each script!

### File-Services

SSfSBT provides a variety of file services that can read from and write to various files used in bioinformatics.
They are located in the file_services folder. Each file service is a class providing class methods.<br> 
Their read()-methods are generators, yielding dictionaries.<br>
Their write()-methods accept iterables of dictionaries.

| File type     | Can read | Can write | Additional info |
|---------------|----------|-----------|-----------------|
| FASTA         | ✅       | ✅       | Sequences
| FASTQ         | ✅       | ✅       | Sequences with qualities
| PAF           | ✅       | ✅       | Pairwise sequence alignments from [Minimap2](https://github.com/lh3/minimap2)
| SAM           | ✅       | ✅       | Pairwise sequence alignments from basically any other alignment tool
| BCALM (FASTA) | ✅       | ❌       | De Bruijn Graph from [BCALM](https://github.com/GATB/bcalm)
| FASTG (FASTA) | ✅       | ❌       | De Bruijn Graph from [SPAdes](https://github.com/ablab/spades)
| GFA           | ✅       | ❌       | Assembly graphs, currently only reads segments
