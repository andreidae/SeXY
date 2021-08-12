## $1 Autosomes.bed
## $2 Satsumaname
## $3 X_trim.bed
## $4 Output_folder
## $5 Bamfile
## $6 Output_prefix

samtools depth -a -q 25 -Q 25 -b $1 $5  > $4/$6_Autosomes_depth.txt
samtools depth -a -q 25 -Q 25 -b $3 $5  > $4/$6_X_trim_depth.txt
for num in {1..10}
do
shuf -n10000000 $4/$6_Autosomes_depth.txt | awk '{sum=sum+$3}END{print sum/10000000}' >> $4/$6_Autosomes_bootstrap.txt
shuf -n10000000 $4/$6_X_trim_depth.txt | awk '{sum=sum+$3}END{print sum/10000000}' >> $4/$6_X_trim_bootstrap.txt
done
paste $4/$6_Autosomes_bootstrap.txt $4/$6_X_trim_bootstrap.txt | awk  '{print $2/$1}' > $4/$6_ratios.txt
rm $4/$6_Autosomes_depth.txt $4/$6_X_trim_depth.txt
