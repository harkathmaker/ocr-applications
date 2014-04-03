# Trim the raw format to make it python-friendly
grep PROFILING $1.log|sed 's/(//g'|sed 's/)//g'|cut -d\  -f2,4,6,10,12,14,18,22 | sed 's/0x//g' > $1.stage1

# Extract all the EDTs
cut -d\  -f1 $1.stage1|sort|uniq > edts

# Extract func names
# nm $1.exe |grep -i " t "| cut -d\  -f1,3 > $1.funcs
nm $1.exe > $1.funcs

# Output a sed script to convert EDTs to Func names
for i in `cat edts`; do NAME=`grep $i $1.funcs|cut -d\  -f3`; echo "s/$i/$NAME/"; done > edt_funcs.sed

# Run it through python to condense DBs into 1 line per EDT
./condense.py $1.stage1 $1.stage2

sort $1.stage2|uniq > unique_edts
IFS=$'\n' read -d '' -r -a lines < unique_edts
LINES=`wc -l unique_edts|cut -d\  -f1`
LINES=`echo $LINES-1|bc`
for i in `seq 0 $LINES`
do
   TMP=`grep -c "${lines[$i]}" $1.stage2|cut -d\: -f2`
   echo ${lines[$i]} $TMP >> tmp
done

echo 'EDT,Total Size (bytes),Instruction Count,Floating Point Ops,Reads(bytes),Writes(bytes),Number of invocations'
cat tmp|sed -f edt_funcs.sed|sed -e 's/ /,/g'

# Cleanup

rm $1.stage* $1.funcs
rm edts unique_edts tmp edt_funcs.sed
