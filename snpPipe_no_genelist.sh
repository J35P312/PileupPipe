#first argument - input directory
#second argument - output directory
#WARNING: do not set input and output directory to the same directory!!!

for file in $(ls $1/*bam)
do
    filename=$(basename "$file")
    extension="${filename##*.}"
    filename="${filename%.*}"
    sbatch SubmitGeneList.sh $file $2 config.conf
done
