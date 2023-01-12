#!/bin/bash
#SBATCH -p G1Part_sce   
#SBATCH -N 1            
#SBATCH -n 1            
#SBATCH -c 56           
ulimit -u 10240
export PATH=/es01/paratera/sce2335/Documents/bin:/es01/paratera/sce2335/Documents/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/es01/paratera/sce2335/.local/bin:/es01/paratera/sce2335/bin  
matlab  < path_int_fsd_script.m