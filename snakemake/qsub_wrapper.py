#!/usr/bin/env python

import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

defaults = {'mem': 8, 'runtime': 24, 'ncores': 1}
for k, v in defaults.items():
    if k not in job_properties['resources']:
        job_properties['resources'][k] = v

params = job_properties['resources']

qsub_cmd = (f'qsub -l h_vmem={params["mem"]}G '
            f'-pe smp {params["ncores"]} -binding linear:{params["ncores"]} '
            f'-l h_rt={params["runtime"]}:00:00 '
            f'-o logs/{job_properties["rule"]}/ '
            f'-e logs/{job_properties["rule"]}/ '
            f'-N {job_properties["rule"]} '
            f'-cwd -V -j y -sync y '
            f'{jobscript}')

qsub_cmd += ' | tail -2 | cut -d " " -f 3'
os.system(qsub_cmd)
