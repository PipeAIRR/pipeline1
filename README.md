**A pipeline for BCR repertoire libraries from the form of - UMI Barcoded Illumina MiSeq 2x250 BCR mRNA. that were produced in the same fashion as those in [sterm et al. 2014](https://pubmed.ncbi.nlm.nih.gov/25100741/)**


<u>Sequence library:</u>

The sequences were amplified using CPrimers and VPrimers. They were sequences with Illumins MiSeq 2x250. 

Each 250 base-pair read was sequenced from one end of the target cDNA, so that the two reads together cover the entire variable region of the Ig heavy chain. The V(D)J reading frame proceeds from the start of read 2 to the start of read 1. Read 1 is in the opposite orientation (reverse complement), and contains a 15 nucleotide UMI barcode preceding the C-region primer sequence.

<u>Input files:</u>

1. The read can downloaded from the NCBI Sequence Read Archive under BioProject accession ID: PRJNA248475 or downloaded or downloaded first 25,000 sequences of sample M12 (accession: SRR1383456) using fastq-dump from the SRA Toolkit
2. The primers sequences available online at the table below.

<u>Output files:</u>

1. sample_collapse-unique.fastq
2. sample_atleast-2.fastq
3. log and log tab file for each step.


<u>Sequence processing:</u>

* Quality control, UMI annotation and primer masking
	1. FilterSeq quality
	2. MaskPrimer score
* Generation of UMI consensus sequences
	3. PairSeq
	4. BuildConsensus
* Paired-end assembly of UMI consensus sequences
	5. PairSeq	
	6. AssemblePairs align
* Deduplication and filtering
	7. ParseHeaders collapse
	8. CollapseSeq
	9. SplitSeq group



**CPrimers**

| Header     | Primer |
| ----------- | ----------- |
| IGHG   | CSGATGGGCCCTTGGTGG       |
| IGHM   |GGGTTGGGGCGGATGCAC        |
| IGHA   | CCTTGGGGCTGGTCGGGG       |
| IGHD   | CATCCGGAGCCTTGGTGG       |
| IGHE   | CGGATGGGCTCTGTGTGG       |


**VPrimers**

| Header     | Primer |
| ----------- | ----------- |
| IGHV_1   | GAGGTGCAGCTGGTGGAGTCCG        |
| IGHV_2   | GAGGTGCAGCTGGTGGAGTCTGG       |
| IGHV_3   | CAGGTCACCTTGAGGGAGTCTGGTCC    |
| IGHV_4   | CAGGTTCAGCTGTTGCAGCCTGG       |
| IGHV_5   | CAGGTGCAGCTACAGCAGTGGGG       |
| IGHV_6   | CAGGTGCAGCTGGTGCAATCTGG       |
| IGHV_7   | CAGGTCACCTTGAAGGAGTCTGGTCC    |
| IGHV_8   | CAGGTCCAGCTGGTACAGTCTGGG      |
| IGHV_9   | CAGGACCAGTTGGTGCAGTCTGGG      |
| IGHV_10  | CAGGTGCAGCTGGTGGAGTCTGG       |
| IGHV_11   | CAGATGCAGCTGGTGCAGTCTGG       |
| IGHV_12   | CAAATGCAGCTGGTGCAGTCTGGG      |
| IGHV_13   | GAAGTGCAGCTGGTGGAGTCTGGG      |
| IGHV_14   | CAGGTGCAGCTGGTGCAGTCTG        |
| IGHV_15   | GAGGATCAGCTGGTGGAGTCTGGG      |
| IGHV_16   | GAGGTCCAGCTGGTACAGTCTGGG      |
| IGHV_17   | CAGCTGCAGCTGCAGGAGTCC         |
| IGHV_20   | CAGGTACAGCTGCAGCAGTCAGGT      |
|IGHV_21	|GAGATGCAGCTGGTGGAGTCTGGG	|
|IGHV_22	|GAGGTGCATCTGGTGGAGTCTGGG	|
|IGHV_25	|CAGGTCCAACTGGTGTAGTCTGGAGC	|
|IGHV_26	|GAGGTGCAGCTGGTGCAGTCTG		|
|IGHV_27	|GAGGTGCAGCTGGTGGAGTCTCG	|
|IGHV_28	|GAGGTTCAGCTGGTGCAGTCTGGG	|
|IGHV_30	|GAAGTGCAGCTGGTGCAGTCTGG	|
|IGHV_31	|GAGGTGGAGCTGATAGAGTCCATAGA	|
|IGHV_32	|ACAGTGCAGCTGGTGGAGTCTGG	|
|IGHV_33	|GAGGTGCAGCTGGAGGAGTCTGG	|
|IGHV_34	|GAGGTACAGCTGGTGGAGTCTGAAGA	|
|IGHV_35	|CAGGTGCAGCTGCAGGAGTCG		|
|IGHV_36	|GAGGTGCAGCTGTTGGAGTCTGGG	|
|IGHV_37	|CAGGTGCAGCTGGGGCAGTC		|
|IGHV_38	|CAGCTGCAGCTGCAGGAGTCG		|
|IGHV_39	|CAGGTTCAGCTGGTGCAGTCTGGA	|
|IGHV_40	|GAGGTGCAGCTGGTAGAGTCTGGG	|
|IGHV_41	|CAGGTCCAGCTTGTGCAGTCTGGG	|
|IGHV_42	|GAGGTGCAGCTGTTGCAGTCTGC	|
|IGHV_43	|GAGGTACAACTGGTGGAGTCTGGGGG	|
|IGHV_44	|CAGATCACCTTGAAGGAGTCTGGTCC	|
|IGHV_45	|CAGGTACAGCTGATGCAGTCTGGGG	|
	



