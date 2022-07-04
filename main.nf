#! /usr/bin/env nextflow
nextflow.enable.dsl=2



include { initvars } from './settings.nf'

include { bbduk; fastqc } from './cleanreads.nf'

