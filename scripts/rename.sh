#/bin/bash
# convert .raw files with ???s into .raw files without them

if [[ $# -lt 3 ]]
  then
    echo "Usage: rename .filetype OldExpression ReplacementExpression"
    exit 1
fi

for file in *$1;
do
  #echo $file
  echo $file "->" `basename $file $1`$1;
  cp $file `basename $file $1`.tmp;
  sed "s/$2/$3/g" `basename $file $1`.tmp >`basename $file $type.tmp`$type;
  #rm -f raw2dat.tmp;
done
