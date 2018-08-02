#!/bin/bash
# R script running table 1 nohup ./bash1.sh > nohup1.out&
# remember: chmod u+x bash1.sh
for i in {1..6}
{
echo "Code $i running"
time R CMD BATCH --slave --no-save --no-restore diffq_exam$i.R diffq_exam$i.Rout >/dev/null 2>&1
}

