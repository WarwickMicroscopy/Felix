#!/bin/bash
# This script will iteratively convert every .bin output image in a directory to a
# a given image format
#
# Be in the images directory and specify image output file extension
# e.g. '../../../scripts/convert.sh jpg'
# '.jpg','jpg','gif' etc. will all work (dot can be included or not)
#
# it reads the pixel size from the name assuming filenames formatted like the following:
# GaAs_085nm_070x070_-2-2-4.bin
# where the 3rd term 070x070 is the pixel siz

img_extension=$1
if ! [[ ${img_extension:0:1} =~ "." ]]
then img_extension=".$img_extension"; fi
for file in $(ls)
do
  if [[ ${file: -4} == ".bin" ]]
  then
    IFS='_' read -r -a filearray <<< ${file}
    echo converting $file to "${file%.bin}""$img_extension"
    convert -size ${filearray[2]} -depth 64 -negate \
    -define quantum:format=floating-point \
    -define quantum:scale=65535.0 \
    -endian lsb GRAY:$file "${file%.bin}""$img_extension"
  fi
done

