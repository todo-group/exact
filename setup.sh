#!/bin/sh

wget -O - https://github.com/todo-group/standards/archive/develop.tar.gz | tar zxf -

DIRS="config integral lattice lse"
for d in $DIRS; do
  mkdir -p $d
  mv standards-develop/$d/* $d/
done
rm -rf develop.tar.gz standards-develop
