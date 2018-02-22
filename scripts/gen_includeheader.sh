DIR=$1
DESCR=$2
HEADER=$3

#Capture all top level headers in subdirectories

exec> $HEADER

echo "#ifndef _CPSFIT_"$DESCR"_H_"
echo "#define _CPSFIT_"$DESCR"_H_"

for i in $(ls $DIR/*.h); do
    echo "#include<"$i">"
done

echo "#endif"
