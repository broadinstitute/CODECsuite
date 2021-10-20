#!/bin/bash
# properties = {properties}
. /broad/software/free/Linux/redhat_6_x86_64/pkgs/anaconda3_5.0.1/etc/profile.d/conda.sh;
conda deactivate && conda activate py36;
source /broad/software/scripts/useuse;
reuse -q .fastp-0.20.0;
reuse -q GCC-5.2;
reuse -q R-4.1;
which python;
{exec_job}
