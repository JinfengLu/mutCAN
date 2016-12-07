Program Description
  mutCAN.pl processes the alignment of query files of mutant & wild-type libraries to reference genome,
exports the coverage depth & nucleotide frequency information for each single genomic position, and calls
candidate mutations by comparing the difference between the two libraries, and also provides reads at the
candidate sites for further confirmation.

Platform & Requirement
  Platform: Linux system
  Requirement: Bowtie
  Install: 
  	Unpack and put the following files to a same folder:
  	-------------------
	mutCAN.pl (main)
    	buildcan.pl
    	can_comparison.pl
   	can_confirmation.pl
    	can_generation.pl
    	mapping_pair.pl
    	mapping_single.pl
	-------------------

Contact & Version
  Author: Jinfeng Lu, jinfeng.lu@ucr.edu
  Contact author: Shou-wei Ding, shou-wei.ding@ucr.edu
  Version: 1.0, Date: 08-23-2016

Usage steps
  (0) Run (optional)
	$ perl mutCAN.pl -h/--help
      to show the usage message
  
  (1) Run
  	$ perl mutCAN.pl --generation [options]
      to align reads & generate files of coverage depth & nucleotide frequency (.can) in /CAN/ folder.
      # A log file (log.txt) recording the mapping information will also be generated.

  (2) Run
	$ perl mutCAN.pl --comparison [options]
      to call mutation candidates by comparing the two .can files in /CAN/ folder.
      # This step can be repeated as long as the two .can files exist in /CAN/ folder.
      # When the chromosome names are different between the reference & annotation file, a correction file will be needed.
      # The correction file contains lines for the chromosomes with different names.
      # First column: chromosome name in annotation file
      # Second column: chromosome name in reference file
      # See /test/Chr_correction for an example applied for Arabidopsis TAIR10 reference & TAIR10 GFF3.

  (3) Run
	$ perl mutCAN.pl --confirmation [options]
      to retrieve reads at the candidate site for confirmation.
      # This step can be repeated given the query reads and the reference genome.

Test
  Testing reads, genome reference, annotation file & correction file are provided under /test/ folder.
  The two artifical libraries are supposed to be from bulked wild-type and mutant pools by 3:1 segregant ratio,
  and the causal mutations are a deletion at Chr1 791-800 & an SNP at ChrC 790 (C>G) which are tightly linked.
  To test the program, run the commands as follows:
	$ perl mutCAN.pl --generation -r <path/>ArabidopsisTAIR10_genome_sub.fa -m <path/>TestReads_mutant_pair1.fastq <path/>TestReads_mutant_pair2.fastq -w <path/>TestReads_wildtype_pair1.fastq <path/>TestReads_wildtype_pair2.fastq
        $ perl mutCAN.pl --comparison -a <path/>TAIR10_GFF3_genes_transposons_sub.gff --cor <path/>Chr_correction --smc 10 --swc 10 --smf 0.9 1 --swf 0.5 0.8 --sjud 0 --dmc 0.5 --dwc 10 --imc 0.5 --iwc 10
        $ perl mutCAN.pl --confirmation -r ./test/ArabidopsisTAIR10_genome_sub.fa -w ./test/TestReads_wildtype_pair1.fastq ./test/TestReads_wildtype_pair2.fastq -m ./test/TestReads_mutant_pair1.fastq ./test/TestReads_mutant_pair2.fastq -c Chr1 --reg 794 800