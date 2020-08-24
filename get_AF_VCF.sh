#!/bin/bash

#Works with the following code:

#find . -name '*.vcf.gz' -execdir zcat {} \; | /path/to/get_AF_VCF.sh > /path/to/output_file.tsv

cat $vcf_file | awk 'BEGIN{FS=" "; print "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAF_normal\tAF_cancer\tDONOR"} {if($1 ~ /INDIVIDUAL/) {split($1,samplename,",");split(samplename[1],samplename2,"=")} else if($7=="PASS"){split($10,normal,":");split($11,cancer,":"); {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,normal[3],cancer[3],samplename2[3])}}}' 

