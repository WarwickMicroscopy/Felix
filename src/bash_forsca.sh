#!/bin/bash
rm felix_temp_Doyle.sca

i=1
while [ $i -lt 103 ];
do

echo -n "DATA DoyleAndTurner($i,8)/" >> felix_temp_Doyle.sca
sed "$i!d" felix.sca | cut -c 336- | sed -r 's/\s+//g' | fold -w15 | paste -sd, - >> felix_temp_Doyle.sca
let i=i+1
done
