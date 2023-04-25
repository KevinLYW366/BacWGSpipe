# Directory of conda environment tmhmm which could be created using workflow/env/tmhmm.yaml
# e.g. /home/usr/anaconda3/envs/tmhmm
envDIR="/lvdata/lyw/anaconda3/envs/tmhmm"
# Directory of tmhmm executable files
# Visit https://services.healthtech.dtu.dk/services/TMHMM-2.0/ to download
# e.g. /home/usr/software/TMHMM/tmhmm-2.0c/bin
binDIR="/lvdata/lyw/software/TMHMM/tmhmm-2.0c/bin"

# Add path to $PATH
export PATH=${envDIR}/bin:$PATH
export PATH=${binDIR}:$PATH
