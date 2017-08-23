#first argument - input directory
#second argument - gene list file, a text file containing one gene symbol per line
#third argument - output directory
#WARNING: do not set input and output directory to the same directory!!!

for file in $(ls $1/*bam)
do
    filename=$(basename "$file")
    extension="${filename##*.}"
    filename="${filename%.*}"
    sbatch delly.sh $1 $3 $2 
done
