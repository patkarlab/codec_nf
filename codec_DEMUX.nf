#!/usr/bin/env nextflow
nextflow.enable.dsl=2

"mkdir Coverview".execute()


log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=

Sample list: ${params.input}
Sequences in:${params.sequences}

"""
process demux{
	
	input:
		tuple val (Sample), val (IndexBarcode1), val (IndexBarcode2)
	output:
		tuple val (Sample), file("*1.fastq.gz"), file("*2.fastq.gz")
	script:
	"""
	${params.demux_sheet} ${Sample} ${IndexBarcode1} ${IndexBarcode2}
	${params.codec} demux  -1 ${params.sequences}/*_S1_R1_001.fastq.gz  -2 ${params.sequences}/*_S1_R2_001.fastq.gz -p demux.csv -o demux_outprefix
	"""
}

process trim {
	input:
		tuple val (Sample), file (forward), file(reverse)
	output:
		tuple val (Sample), file("*.trim.bam")
	script:
	"""
	${params.codec} trim -1 ${forward} -2 ${reverse} -o ${Sample}_trim_outprefix -u 3 -U 3 -f 2 -t 2 -s ${Sample}
	sleep 1s
	"""
}

process AlignRawTrimmed {
	input:
		tuple val (Sample), file(trimbam)
	output:
		tuple val (Sample), file("*.aligned.bam")
	script:
	"""
	${params.samtools} fastq ${trimbam} | bwa mem -K 100000000 -t 20 -p -Y ${params.genome} - | ${params.samtools} view -bS - -o ${Sample}.trim.aligned.bam
	sleep 1s 	
	"""

}

process ZipperBamAlignment {
	input:
		tuple val (Sample), file(TrimBam), file(AlignedBam)
	output:
		tuple val (Sample), file("*.ZipperBam.bam"), file("*.ZipperBam.bam.bai")
	script:
	"""
	${params.fgbio} --compression 0 --async-io ZipperBams -i ${AlignedBam} --unmapped ${TrimBam} --ref ${params.genome} -o ${Sample}.zipped 
	${params.samtools} sort ${Sample}.zipped -o ${Sample}.ZipperBam.bam -O BAM -@ 20 
	${params.samtools} index ${Sample}.ZipperBam.bam -@ 20
	sleep 5s
	"""
}

process MergeSplit {
	input:
		tuple val (Sample), file(ZipperBambam), file (ZipperBambambai)
	output:
		tuple val (Sample), file("*.MergeSplit.bam"), file("*.MergeSplit.bam.bai")
	script:
	"""
	${params.samtools} merge -@ 20 ${Sample}.MergeSplit.bam  ${ZipperBambam} 
	${params.samtools} index ${Sample}.MergeSplit.bam  -@ 20 ${Sample}.MergeSplit.bam.bai
	"""
}

process ReplaceRawReadGroup {
	input:
		tuple val (Sample), file (MergeSplitsortedbam),  file(MergeSplitsortedbambai)
	output:
		tuple val (Sample), file("*.raw.replacerg.bam")
	script:
	"""
	${params.java_path}java -Xmx8g -jar ${params.picard_path} AddOrReplaceReadGroups I=${MergeSplitsortedbam} O=${Sample}.raw.replacerg.bam CREATE_INDEX=true RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${Sample}

	"""
}

process MarkRawDuplicates {
	input:
		tuple val (Sample), file(RawReplacergBam)
	output:
		tuple val (Sample), file("*.marked.bam"), file("*.marked.bai")
	script:
	"""
	${params.java_path}java -Xmx8g -jar ${params.picard_path} MarkDuplicates I=${RawReplacergBam}  O=${Sample}.marked.bam  M=${Sample}.marked_dup_metrics.txt CREATE_INDEX=true TAG_DUPLICATE_SET_MEMBERS=true TAGGING_POLICY=All
	"""
}

process CollectInsertSizeMetrics {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.txt'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.pdf'
	input:
		tuple val (Sample), file(markdupbam), file(markdupbambai)
	output:
		tuple val (Sample), file ("*.txt"), file("*.pdf")
	script:
	"""
	${params.java_path}java -Xmx8g -jar ${params.picard_path} CollectInsertSizeMetrics I=${markdupbam} O=${Sample}.InsertSizeMetrics.txt  H=${Sample}.insert_size_histogram.pdf  W=600 DEVIATIONS=100
	sleep 1s
	"""
}

process GroupReadByUMI {
	input:
		tuple val (Sample), file(markdupbam), file(markdupbambai)
	output:
		tuple val (Sample), file ("*.GroupedByUmi.bam")
	script:
	"""
	${params.fgbio} --compression 1 --async-io GroupReadsByUmi -i ${markdupbam} -o ${Sample}.GroupedByUmi.bam -f ${Sample}.umiHistogram.txt -m 0 --strategy=paired
	"""
}

process FgbioCollapseReadFamilies {
	input:
		tuple val (Sample), file (GroupReadsByUmi)
	output:
		tuple  val (Sample), file ("*.mol_consensus.bam")
	script:
	"""
	${params.fgbio} --compression 1 CallMolecularConsensusReads -i ${GroupReadsByUmi} -o ${Sample}.mol_consensus.bam -p ${Sample} --threads 20 --consensus-call-overlapping-bases false -M 1
	"""
}

process AlignMolecularConsensusReads {
	input:
		tuple val (Sample), file (mol_consensusbam)
	output:
		tuple val (Sample),file ("*.mol_consensus.aligned.bam")
	script:
	"""
	${params.samtools} fastq ${mol_consensusbam} | bwa mem -K 100000000 -t 20 -p -Y ${params.genome} - > ${Sample}.mol_consensus.aligned.bam
	"""
}
process MergeAndSortMoleculeConsensusReads {
	input:
		tuple val (Sample), file(mol_consensusbam), file (molconsensus_alignedbam)
	output:
		tuple val (Sample), file("*.sorted.bam"), file("*.sorted.bam.bai")
	script:
	"""
	${params.fgbio} --compression 0 --async-io ZipperBams -i ${molconsensus_alignedbam} --unmapped ${mol_consensusbam} --ref ${params.genome} --tags-to-reverse Consensus --tags-to-revcomp Consensus -o ${Sample}.MoleculeConsensusRead.bam
	
	${params.samtools} sort ${Sample}.MoleculeConsensusRead.bam -o ${Sample}.MergeConsensusRead.sorted.bam -O BAM -@ 20
	${params.samtools} index ${Sample}.MergeConsensusRead.sorted.bam -@ 20
	sleep 1s
	"""

}

process Call {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.variants_called.txt'
	input:
		tuple val (Sample), file(markdupbam), file(markdupbambai)
	output:
		tuple val (Sample), file ("*.variants_called.txt")
	script:
	"""
	${params.codec} call -b ${markdupbam}  -L ${params.bedfile}  -r ${params.genome}  -p lenient -o ${Sample}
	sleep 1s
	"""
}
process MolConsReadsCall {
	input:
		tuple val (Sample), file(MolConsReadssortbam), file(MolConsReadssortbambai)
	output:
		tuple val (Sample), file ("*.variants_called.txt")
	script:
	"""
	${params.codec} call -b ${MolConsReadssortbam}  -L ${params.bedfile}  -r ${params.genome}  -p lenient -o ${Sample}
	"""

}
workflow CODEC {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.set { Sample }

	main:
		demux(Sample) 
		trim (demux.out)
		AlignRawTrimmed(trim.out)
		ZipperBamAlignment(trim.out.join(AlignRawTrimmed.out))
		MergeSplit(ZipperBamAlignment.out)
		ReplaceRawReadGroup(MergeSplit.out)
		MarkRawDuplicates(ReplaceRawReadGroup.out)

		CollectInsertSizeMetrics(MarkRawDuplicates.out)
		GroupReadByUMI(MarkRawDuplicates.out)
		FgbioCollapseReadFamilies(GroupReadByUMI.out)
		AlignMolecularConsensusReads(FgbioCollapseReadFamilies.out)
		MergeAndSortMoleculeConsensusReads(FgbioCollapseReadFamilies.out.join(AlignMolecularConsensusReads.out))
		
		Call(MarkRawDuplicates.out)
		//MolConsReadsCall(MergeAndSortMoleculeConsensusReads.out)
}
workflow.onComplete {
		log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
