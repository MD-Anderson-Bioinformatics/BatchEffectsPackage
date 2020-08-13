#!/bin/bash

./00_clear.bash

echo "move to app dir"
cd ../apps

echo "compile BEVIndex"
cd BEVIndex
ant clean jar
cd ..

echo "compile DscJava"
cd DscJava
ant clean jar
cd ..

echo "compile LegendJava"
cd LegendJava
ant clean jar
cd ..

echo "compile ReadRJava"
cd ReadRJava
ant clean jar
cd ..

echo "list dist files"
ls -l */dist

echo "done"
