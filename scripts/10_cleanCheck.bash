#!/bin/bash

echo "move to app dir"
cd ../apps

echo "remove .Rcheck directories"
rm -r *.Rcheck

echo "remove R packages"
rm *.tar.gz

echo "lsit apps directory"
ls -l

echo "done"
