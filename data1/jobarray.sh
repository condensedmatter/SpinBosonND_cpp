#!/bin/bash
echo $SGE_TASK_ID
./main_proj_2 $SGE_TASK_ID
