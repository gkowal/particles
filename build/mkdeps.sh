#!/bin/sh
#
# A script to generate compile dependencies for Fortran files.
#
# Usage:
#
#   mkdep.sh (FORTRAN_SOURCE_FILE_PATH) (OBJECT_FILE_PATH)
#
# Example:
#
#   mkdep.sh ../src ./obj
#
# Requirements: awk, sort, uniq
#
files=`ls $1/*.F90`
for src in $files; do
  deps=`awk '/^\s*[Uu][Ss][Ee] / {gsub ( "[:,]","" ) ; print $2}' $src | sort | uniq`
  fname=`basename $src .F90`
  output="$2/${fname}.o: $src"
  for dep in $deps; do
    [ -z "${files##*$dep.F90*}" ] && output="${output} $2/$dep.o"
  done
  echo $output
done
