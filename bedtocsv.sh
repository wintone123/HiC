input="/mnt/c/hic/test3/test.bed"
output="/mnt/c/hic/test3/test2_fil.csv"
score=40

# bed to csv
awk -v S=$score 'BEGIN{print"chr1""\t""start1""\t""end1""\t""strand1""\t""chr2""\t""start2""\t""end2""\t""strand2"}NR%2==1{A=$1"\t"$2"\t"$3;B=$5;C=$6};NR%2==0{if(B>=S && $5>=S) print A"\t"C"\t"$1"\t"$2"\t"$3"\t"$6}' $input > $output

# csv split
