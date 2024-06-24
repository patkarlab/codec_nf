#!/usr/bin/env bash


nextflow -c /home/diagnostics/pipelines/codecPipeline/nextflow.config run codec_DEMUX.nf -entry CODEC \
--bedfile /home/diagnostics/pipelines/codecPipeline/bedfile/Probes-XGEN_9857C1725CE84D37803C380036834C74_g \
-resume -bg 

