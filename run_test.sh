
if [ -d .gnmgit/ ]
then
   rm -rf .gnmgit/
fi
./genomegit init

./genomegit add --threads=8 /home/joe/dataSets/fruitfly/fruit_fly_75_to_1215.vcf
./genomegit add --threads=8 /home/joe/dataSets/fruitfly/GCA_002310755.1_ASM231075v1_genomic.fa
/usr/bin/time -v ./genomegit add --t=8 --a=1 --s=5000  /home/joe/dataSets/fruitfly/GCA_002310775.1_ASM231077v1_genomic.fa
./genomegit get --dataset=Variants
./genomegit discarded fruit_fly_75_to_1215.vcf
