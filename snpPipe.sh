#first argument - input directory
#second argument - gene list file, a text file containing one gene symbol per line
#third argument - output directory
#argument four - config
#WARNING: do not set input and output directory to the same directory!!!

for file in $(ls $1/*bam)
do
    filename=$(basename "$file")
    extension="${filename##*.}"
    filename="${filename%.*}"
    #echo "sbatch delly.sh $file $3 $2" 
    sbatch SubmitGeneList.sh $file $3 $2 $4
done
