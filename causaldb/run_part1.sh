#!/usr/bin/env bash
set -euo pipefail

QTL_TYPES=(e)
#QTL_INDEPS=(PrimQTL IndepQTL)
QTL_INDEPS=(IndepQTL)
#METHODS=(abf finemap paintor caviarbf susie polyfun_finemap polyfun_susie)
METHODS=(polyfun_susie) 
#CAUSALDB_INDEPS=(PrimGwas IndepGwas)
CAUSALDB_INDEPS=(IndepGwas)
#VAR_TYPES=(snp all)
VAR_TYPES=(all)


for QTL_TYPE in ${QTL_TYPES[@]}; do
  for QTL_INDEP in ${QTL_INDEPS[@]}; do
    for METHOD in ${METHODS[@]}; do
      for VAR_TYPE in ${VAR_TYPES[@]}; do
        for CAUSALDB_INDEP in ${CAUSALDB_INDEPS[@]}; do
	  echo $QTL_TYPE
	  echo $QTL_INDEP
	  echo $METHOD
	  echo $VAR_TYPE
	  echo $CAUSALDB_INDEP
	  echo logs/$QTL_TYPE.$QTL_INDEP.$METHOD.$VAR_TYPE.$CAUSALDB_INDEP.log
          bash pipeline/part1/01-05.sh \
            --qtl_type $QTL_TYPE \
            --qtl_indep $QTL_INDEP \
            --method $METHOD \
            --var_type $VAR_TYPE \
            --causaldb_indep $CAUSALDB_INDEP \
     	    > logs/$QTL_TYPE.$QTL_INDEP.$METHOD.$VAR_TYPE.$CAUSALDB_INDEP.log 2>&1
        done
      done
    done
  done
done
