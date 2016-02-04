#MultiGeMS 1.0

##Overview

Multi-sample Genotype Model Selection (MultiGeMS) is a multiple sample single 
nucleotide variant (SNV) caller that works with alignment files of 
high-throughput sequencing (HTS) data. MultiGeMS calls SNVs based on a 
statistical model selection procedure and accounts for enzymatic substitution 
sequencing errors. 

##Compile

MultiGeMS can be compiled using a GCC compiler with C++11 support, by running:

$ make

##Input

MultiGeMS accepts a text file, listing on seperate lines, paths of SAMtools 
pileup format files. To convert a SAM/BAM alignment file into the pileup 
format, users can use the SAMtools mpileup procedure with option -s.

##Filter

Alignment file reads with undesirable characteristics can be filtered before 
running MultiGeMS. For an explanation why filtering may be desirable and for a 
brief tutorial on how to filter SAM/BAM alignment files using SAMtools view, 
please see the PDF document entitled "Pre-Filtering Alignment Files" available 
at https://github.com/cui-lab/multigems.

##Usage

multigems -i pileuplist.txt -o multigems.out [OPTIONS]

## Options

-b INT   minimum base-calling quality score considered, default is 17

-m INT   minimum mapping quality score considered, default is 20

-s FLOAT maximum likelihood computing steps float value, smaller is slower yet 
         more precise, default is 0.01

-e FLOAT EM algorithm convergence threshold, smaller is slower yet more 
         precise, default is 0.001

-M INT   maximum number of bases to be considered from each sample of each 
         site, 0 indicates unbounded, default is 255

-f FLOAT lFDR SNV threshold, between 0 and 1, default is 0.1  

-C INT   number of sites to be analyzed per analysis cycle, smaller is slower 
         in general and uses less RAM, default is 200

-t INT   number of threads, can be used in conjunction with -C option for 
         multi-thread efficiency, default is 1

-n FLOAT site non-reference allele proportion filter (sites where all samples 
         have a non-reference proportion less than this filter are not 
         analyzed), smaller is slower and more susceptible to false positive 
         SNV calls, higher is faster and more susceptible to false negative SNV 
         calls, between 0 and 1, default is 0.2

-l FLOAT sample deletion placeholder proportion filter (samples with pileup 
         deletion placeholder proportions greater than filter are not 
         analyzed), between 0 and 1, default is 0.05 

##Output

The MultiGeMS output is similar to that of the Variant Call Format (VCF) file 
format. Meta-information lines are provided at the beginning of each output. 
Only sites less than the user-selected lFDR SNV threshold are output and are 
considered SNV calls. An lFDR output of "-1" indicates that the site was
analyzed by MultiGeMS, but not enough information was retained post-filtering
to estimate the lFDR value for the site.

##Contact

Xinping Cui

xinping.cui@ucr.edu

https://sites.google.com/a/bioinformatics.ucr.edu/xinping-cui/
