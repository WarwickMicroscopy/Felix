#!/bin/bash
# This script will iteratively convert every felix .bin in a directory to a
# a given output image format
#
# Be in the images directory and specify image output file extension
# e.g. '../../../scripts/convert.sh bmp'
# '.bmp','bmp','gif' etc. will all work

img_extension=$1
if ! [[ ${img_extension:0:1} =~ "." ]]
then img_extension=".$img_extension"; fi
for ofile in $(ls)
do
  if [[ ${ofile: -4} == ".bin" ]]
  then
  echo converting $ofile to "${ofile%.bin}""$img_extension"
  convert -size 128x128 -depth 64 -negate \
  -define quantum:format=floating-point \
  -define quantum:scale=65535.0 \
  -endian lsb GRAY:$ofile "${ofile%.bin}""$img_extension"
  fi
done

