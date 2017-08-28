#!/bin/bash
# This script will iteratively convert every .bin output image in a directory
#
# Be in the images directory 
#
# it reads the pixel size from the name assuming filenames formatted like the following:
# GaAs_085nm_070x070_-2-2-4.bin
# where the 3rd term 070x070 is the pixel siz

img_extension=.jpg
greyscale=0

# (output) img_extension can have '.' or not i.e. '.jpg' or 'jpg'
if ! [[ ${img_extension:0:1} =~ "." ]]
then img_extension=".$img_extension"; fi

for infile in $(ls)
do
  if [[ ${infile: -4} == ".bin" ]]
  then

    # fill array with filename parts delimted by underscores 
    IFS='_' read -r -a filename_parts <<< ${infile}
    outfile="${infile%.bin}""$img_extension"

    echo converting $infile to $outfile
    convert -size ${filename_parts[2]} -depth 64 \
    -define quantum:format=floating-point \
    -define quantum:scale=65535.0 \
    -endian lsb GRAY:$infile $outfile

    if [[ ${greyscale} == 0 ]]
    then
      convert $outfile \
      -colorspace gray \
      \( xc:black xc:red xc:orange xc:yellow xc:green1 xc:cyan xc:blue xc:blueviolet xc:black \
      +append -filter Cubic -resize 600x30! -flop \) \
      -clut $outfile
    fi

  fi
done
