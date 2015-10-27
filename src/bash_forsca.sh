#!/bin/bash
rm Lobato_initial.txt
rm felix_lobato.sca

awk 'NR % 3 !=0' Lobato_scatt_firstdraft.txt > Lobato_initial.txt
awk 'NR % 2 ==0' Lobato_initial.txt > Lobato_second_row.txt
awk 'NR % 2' Lobato_initial.txt | cut -c -89 > Lobato_first_row.txt
awk '{print $0", &"}' Lobato_first_row.txt > Lobato_first_row2.txt

i=1
let j=$i*2

while [ $i -lt 104 ]
do
echo -n "DATA Lobato($i,10)/" >> felix_lobato.sca
sed "$i!d" Lobato_first_row2.txt >> felix_lobato.sca
echo -n "" >> felix_lobato.sca
sed "$i!d" Lobato_second_row.txt >> felix_lobato.sca

let i=i+1
done

tr -s '\t' <felix_lobato.sca | tr '\t' ',' >felix_lobato2.sca 
