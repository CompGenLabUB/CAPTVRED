# virwaste
Analysis pipeline for the VIRWASTE viral metagenomics project.

## Run from cmdline:
```{.sh}
source projectvars.sh
cd $BDIR
/usr/local/bin/nextflow run main.nf 
```

##Files:
* projectvars.sh -> this file is not part of the nextflow environment, is used to set the local directories before launchong nextflow.
* nextflow.config -> Environment variables. This file is the same for all the runs. This file is read by default by nextflo (must be placed in the project directory or in the base directory).
* main.nf 
* 

note/ not used by now: RUNXX.config -> Variables specific for each run. This file mus be modied for each run and placed to the base directory for the run.
note/virwaste.config is no longer used it has been split into nextflow.config and RUNXX.config




