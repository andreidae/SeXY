## $1 Autosomes.bed
## $2 X_trim.bed
## $3 Output_folder
## $4 Bamfile
## $5 Output_prefix

samtools depth -a -q 25 -Q 25 -b $1 $4  > $3/$5_Autosomes_depth.txt
samtools depth -a -q 25 -Q 25 -b $2 $4  > $3/$5_X_trim_depth.txt
for num in {1..10}
do
shuf -n10000000 $3/$5_Autosomes_depth.txt | awk '{sum=sum+$3}END{print sum/10000000}' >> $3/$5_Autosomes_bootstrap.txt
shuf -n10000000 $3/$5_X_trim_depth.txt | awk '{sum=sum+$3}END{print sum/10000000}' >> $3/$5_X_trim_bootstrap.txt
done
paste $4/$6_Autosomes_bootstrap.txt $3/$5_X_trim_bootstrap.txt | awk  '{print $2/$1}' > $3/$5_ratios.txt
rm $3/$5_Autosomes_depth.txt $3/$5_X_trim_depth.txt
