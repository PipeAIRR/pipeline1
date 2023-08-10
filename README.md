A pipeline for BCR repertoire libraries from the form of - UMI Barcoded Illumina MiSeq 2x250 BCR mRNA. that were produced in the same fashion as those in [sterm et al. 2014](https://pubmed.ncbi.nlm.nih.gov/25100741/)


Library preperation and sequencing method:

The sequences were amplified using specific primers of the C-region and non coding V region of the IGH chain. The generated libraries were then sequenced with Illumins MiSeq 2x250. 


Each read contains 250 base-pair and was sequenced from one end of the target cDNA, so that the two reads together cover the entire the full variable region of the Ig heavy chain. The V(D)J reading frame proceeds from the start of read 2 to the start of read 1. Read 1 is in the reverse complement orientation, and contains a 15 nucleotide UMI barcode preceding the C-region primer sequence.


Input files:

1. Two fastq file of paird sequencing
2. primer files

To test the pipeline:

We recommend to test the pipeline using a small example from the original reads that can download using the fastq-dump command.

```bash
fastq-dump --split-files -X 25000 SRR1383456
```

And upload directly to dolphinnext. 


Output files:

1. {sampleName}_collapse-unique.fastq
2. {sampleName}_atleast-2.fastq
3. log tab file for each steps
4. report for some of the steps


Pipeline container:

* Docker: immcantation/suite:4.3.0


Sequence processing steps:

* Quality control, UMI annotation and primer masking

	1. FilterSeq qality
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


Primers used:

* [Stern2014_CPrimers (R1)](https://bitbucket.org/kleinstein/presto/src/master/examples/Stern2014/Stern2014_CPrimers.fasta)



* [Stern2014_VPrimers (R2)](https://bitbucket.org/kleinstein/presto/src/master/examples/Stern2014/Stern2014_VPrimers.fasta)

