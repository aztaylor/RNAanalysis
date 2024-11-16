#! /bin/bash
[ -d "~/organized_fastaqs" ] && echo "Directory: organized_fastaqs exists" || mkdir ~/organized_fastaqs 
FASTAQ_FP=/Users/alec/Documents/Grad_school/Yeung_Lab/RNAanalysis/organized_fastaqs/*
for f in $FASTAQ_FP; do
		base="$(basename $f)"
		echo "Proccessing: $base"
		java -jar /Users/alec/Documents/Grad_school/Yeung_Lab/RNAanalysis/Software/Trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar $f ~/Desktop/Trimmomatic_output ILLUMINACLIP:TrueSeq3-SE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
done
