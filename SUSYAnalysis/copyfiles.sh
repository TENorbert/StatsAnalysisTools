#!/bin/bash
#Simple script to copy files from folder A to folder B
for i in `ls -1 /uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/`
do 
cp $i SusyFiles $i  
#.`date +%m%d%Y`
done
