What's New - September 28 2021

### What's new

[- New tool: Minimap2 now available]

---

## Awk tutorial

x echo "a b c d" | awk '{ print $3 }'
x awk '{ print $3 }' orders.tsv | head
x awk '{ print $0 }' orders.tsv | head  # Note that columns are not zero indexed because 0 refers to the whole line! (also works with just print)
* awk '{ print $1"-->"$2"<--"$3 }' orders.tsv | head  # Concatenating fields with custom delimiters
x awk '{ print $1,$2,$3 }' orders.tsv | head  # Use "," to use the default delimiter

x awk '{ if($3 != "Chicken") print }' orders.tsv | head  # If statements
* awk '{ if($3 != "Chicken" || NR == 1) print }' orders.tsv | head  # Still want header (NR=number of records [i.e. lines])
x EXERCISE: find large orders with 10 or more times the same item: awk '{ if($2 >= 10) print }' orders.tsv > orders_large.tsv

x awk '{ n+=$2 } END { print n }' orders.tsv  # Count total number of items
x awk 'BEGIN{ n=0; }{ n+=$2 } END { print n }' orders.tsv  # BEGIN statement is implied, but can be useful in other applications where don't want it to initialize to 0
x awk 'BEGIN{ print(5/7) }'  # awk doesn't need a file to work
x EXERCISE: find the largest value of column 2 using just awk: awk 'BEGIN{ largest=0; }{ if(NR > 1 && $2>largest)largest=$2 } END { print largest }' orders.tsv


x awk -F "\t" '{ if($3 ~ /Burrito/) myArray[$3]++ } END { for(x in myArray) print x, myArray[x] }' orders.tsv | grep "Burrito"
x awk '{ myArray[$1]++ } END { for(x in myArray) print x }' orders.tsv | head
x awk '{ if(NR > 1) myArray[$1]++ } END { print length(myArray) }' orders.tsv | head
x awk '{ if(NR > 1) myArray[$1]++ } END { print "index", "count"; for (i = 1; i <= length(myArray); i++) { print i, myArray[i] } }' orders.tsv | head
x awk '{ myArray[$1]+=$2 } END { for(x in myArray) if(myArray[x] > 10) print x"\t"myArray[x] }' orders.tsv


* echo "BEGIN{ n=0; }{ n+=\$2 } END { print n }" > script.awk
  awk -f script.awk orders.tsv    # Execute script saved in a file!
* awk '{ if($3 != "Chicken") print > "orders_not_chicken.tsv" }' orders.tsv   # Write to a file within an awk script
  head orders_not_chicken.tsv
* Functions!
awk 'function hello(a, b){ return a + b; }{ print $1, $2, hello($1, $2) }' orders.tsv | head

  awk ' { print($2) }' orders.tsv | head
