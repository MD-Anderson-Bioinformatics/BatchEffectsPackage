#!/bin/bash

./00_clear.bash

echo "move to app dir"
cd ../apps

echo "compile BEVIndex"
cd BEVIndex
ant -f build.xml
cd ..

echo "compile DscJava"
cd DscJava
ant -f build.xml
cd ..

echo "compile LegendJava"
cd LegendJava
ant -f build.xml
cd ..

echo "compile ReadRJava"
cd ReadRJava
ant -f build.xml
cd ..

echo "list dist files"
ls -l */dist

echo "done"
