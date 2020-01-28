cd $1
cat */ranbow.single.eval > result.sing
#cat */ranbow.single.eval | grep "<--"  | cut -f3 | tr -d -c [:digit:]'\n'[:space:] | tr ' ' \\t > result.sing
#cat */ranbow.pair.eval | grep "<--"  | cut -f3 | tr -d -c [:digit:]'\n'[:space:] | tr ' ' \\t > result.pair
