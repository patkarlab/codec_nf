#!/usr/bin/env bash


nextflow -c /home/diagnostics/pipelines/codecPipeline/nextflow.config run codec_DEMUX.nf -entry CODEC -resume -bg
