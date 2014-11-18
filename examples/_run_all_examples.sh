#!/bin/bash

#amber/
#gui/
#parallel_pele/

basedirs="
basinhopping/
basinhopping_no_system_class/
connecting_minima/
disconnectivity_graph/
frozen_degrees_of_freedom/
heisenberg_model/
mindist/
new_potential/
using_the_system_class/
xymodel/
"

homedir=$PWD
bad_files=""

for d in $basedirs; do
  cd $homedir
  echo ""
  echo ""
  echo "going into directory $d"
  cd $d
  if [ $? -ne 0 ]; then
    echo "skipping bad dir $d"
    bad_files="$bad_files $d"
    continue
  fi
  files=*py
  for f in $files; do
    echo "running file: $f"
    python $f > /dev/null
    if [ $? -ne 0 ]; then
      echo "file $f returned an error"
      bad_files="$bad_files $d/$f"
    fi
  done
done

echo "The following files had errors"
for f in $bad_files; do
  echo $f
done
