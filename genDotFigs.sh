#!/bin/bash -e
for i in out/*.dot; do
  name=`basename $i .dot`
  echo Running $i
  sed -i 's/\[ label = "Bus/\[ fillcolor = "#80c0ff", label = "Bus/g' $i $i
  sed -i 's/\[ label = "Gen/\[ fillcolor = "#90ee90", label = "Gen/g' $i $i
  sed -i 's/\[ label = "Load/\[ fillcolor = "#ffa07a", label = "Load/g' $i $i
  #dot -Gdpi=800 -Tpng $i > out/$name.png &
  dot -Tsvg $i > out/$name.svg &
  dot -Tpdf $i > out/$name.pdf &
  wait
done
echo Done!

