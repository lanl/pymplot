
for i in $(ls *.py)
do
	#2to3 -w $i
    yapf -i --style='{based_on_style: pep8, column_limit: 110}' $i
    echo $i
done 
