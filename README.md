# Codec suite

## Description

# 1. CODECSuite: Demultiplexing.
Demultiplexing it requires a sample sheet in csv format which contain Sample and Indexbarcode1,Indexbardcode2.

# 2. CODECSuite: Adapter trimming.
The adapter trimming step removes the adapter sequences from the read and output as uBAM (unmapped BAM format).

-1	R1.fastq
-2	R1.fastq
-o	Output prefix
-u	Num of umi bases to be trimmed from 5'end of read1
-U	Num of umi bases to be trimmed from 5'end of read2
-f	Num of bases to be trimmed at the 5'end of a read after adapter removed 
-t	Num of bases to be trimmed at the 3'end of a read after adapter removed 
-s	Read group sample when output bam

# 3. AlignRawTrimmed
Trimmed uBAM is converted into fastq files . This fastq files are aligned to refrence genome and ouput is mapped BAM.

-K	
-t	number of threads
-p	smart pairing	
-Y	use soft clipping for supplementary alignments

# 4. ZipperBamAlignment
Zips together an unmapped and mapped BAM to transfer metadata into the output BAM.

-i	Mapped BAM.
-u	Unmapped BAM.
-o  Output BAM file.

# 5. ReplaceRawReadGroup
Replace read groups in a BAM file. ( input: bam file from ZipperBamAlignment)

I    	input bam
O    	output bam
RGID 	Read Group ID 
RGPL 	Read Group platform
RGPU 	Read Group platform unit
RGSM 	Read Group sample name Required.

# 6.MarkRawDuplicates
Identifies duplicate reads and tag/mark them in a BAM.

I	input bam file
O	output bam file
M	File to write duplication metrics

# 7. CollectInsertSizeMetrics
Metrics for validating library construction including the insert size distribution.

I			input bam file
O			output bam file
H			File to write insert size Histogram chart.
M			When generating the Histogram, discard any data categories (out of FR, TANDEM, RF)
			that have fewer than this percentage of overall reads. (Range: 0 to 1).

W			Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. 
			Also, when calculating mean and standard deviation, only bins <= Histogram WIDTH wil be included.

DEVIATIONS 		Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. 
			This is done because insert size data typically includes enough anomalous values from chimeras 
			and other artifacts to make the mean and sd grossly misleading regarding the real distribution.

# 8. GroupReadsByUmi
Groups reads together that appear to have come from the same original molecule

-i			input BAM file.
-o			output BAM file.
-f			optional output of tag family size counts.
-strategy	The UMI assignment strategy.

# 9. FgbioCollapseReadFamilies
Calls consensus sequences from reads with the same unique molecular tag.

-i		input BAM file.
-o		ouput BAM file
-p		The Prefix all consensus read names (optional)
-consensus-call-overlapping-bases 		Consensus call overlapping bases in mapped paired end reads

# 10. AlignMolecularConsensusReads
uBAM obtained from FgbioCollapseReadFamilies process is converted into fastq files. 
This fastq files are aligned to refrence genome and ouput is mapped BAM.

# 11. MergeAndSortMoleculeConsensusReads
Zips together an uBAM from FgbioCollapseReadFamilies and mapped BAM AlignMolecularConsensusReads into the output BAM.

-i  Mapped BAM.
-u  Unmapped BAM.
-o  Output BAM file
--tags-to-reverse this will tag rev. on unmapped bam
--tags-to-revcomp 
