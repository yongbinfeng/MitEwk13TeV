#!/bin/sh

# check if all directories contains all the output root files
dirName=ResultsNew2
counterAll=0
counter=0

for i in $(ls $dirName); do
    if [ -e $dirName/${i}/pdfsU1.root ]; then
        if [ -e $dirName/$i/pdfsU2.root ]; then
            if [ -e $dirName/$i/histos.root ]; then
                counter=$((counter + 1))
            else
                echo "Directory $dirName/${i}/histos.root does not exist"
            fi
        else
            echo "Directory $dirName/${i}/pdfsU2.root does not exist"
        fi
    else
        echo "Directory $dirName/${i}/pdfsU1.root does not exist"
    fi
    counterAll=$((counterAll + 1))
done

echo "Number of directories: $counterAll"
echo "Number of directories with all the output files: $counter"
