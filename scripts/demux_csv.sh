#!/usr/bin/bash


SampleName=$1
IndexBarcode1=$2
IndexBarcode2=$3

echo "SampleName,IndexBarcode1,IndexBarcode2" >> demux.csv
echo "${1},${2},${3}" >> demux.csv
