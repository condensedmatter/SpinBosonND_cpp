#!/bin/bash
qsub -cwd -V -N PJ -l h_data=1024M,h_rt=24:00:00 -M $HOME -m bea -t 1501-2000:1 jobarray.sh

