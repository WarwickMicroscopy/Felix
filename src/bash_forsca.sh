#!/bin/bash
rm felix_temp.sca


i=1
while [ $i -lt 103 ];
do

echo -n "DATA Kirkland($i,13)/" >> felix_temp.sca
sed "$i!d" felix.sca | cut -c -208 | sed -r 's/\s+//g' | fold -w15 | paste -sd, - >> felix_temp.sca
let i=i+1
done
