#!/bin/sh

wget -O - https://github.com/todo-group/standards/archive/develop.tar.gz | tar zxf -

mkdir -p config
mv standards-develop/config/* config/

DIRS="integral lattice lse"
for d in $DIRS; do
  mkdir -p exact/$d
  mv standards-develop/$d/* exact/$d/
done
rm -rf develop.tar.gz standards-develop
