#!/usr/bin/env bash

# e.g. bash BacWGSpipe.sh --configfile config/config.yaml -r -p --cores 32 --singularity-args "-B /lvdata/lyw:/lvdata/lyw" -n

#####################################################################################################################

display_help() {
    echo "BacWGSpipe - a Snakemake workflow for a complete analysis of bacterial whole-genome sequencing data"
    echo ""
        echo "Usage: bash BacWGSpipe.sh --configfile configfile.txt [Snakemake options]"
        echo "Note: if --singularity-args is used, args after it should be quoted"
        echo
    echo " Useful Snakemake parameters:"
        echo "    -n,     --dryrun              do not execute anything"
        echo "    -p,     --printshellcmds      print out the shell commands that will be executed"
        echo "    -t,     --timestamp           add a timestamp to all logging output"
        echo "    -c N,   --cores N             use at most N cores in parallel"
        echo "    --ri,   --rerun-incomplete    re-run all jobs the output of which is recognized as incomplete"
        echo "            --singularity-args    pass additional args to singularity (e.g. bind path by '-B /lvdata/lyw:/lvdata/lyw')"
        echo "    -q,     --quiet               do not output certain information"
        echo "            --cleanup-shadow      cleanup old shadow directories which have not been deleted due to failures or power loss"
        echo "            --verbose             print detailed stack traces and detailed operations"
        echo "            --nocolor             do not use a colored output"
    echo
        echo " Full list of parameters:"
    echo "   --help                 show Snakemake help (or snakemake -h)"
        echo
    exit 0
}

# shellcheck disable=SC2166
if [ "$1" == "" -o "$1" == "-h" -o \( "$1" != "--configfile" -a "$1" != "--help" \) ]; then
  display_help
  exit 0
fi

#####################################################################################################################

snakemake -s "$(dirname "$0")/workflow/Snakefile" "--use-singularity" "--use-conda" "$@"

