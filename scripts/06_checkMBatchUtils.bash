#!/bin/bash

echo "move to app dir"
cd ../apps

echo "check MBatchUtils"
R CMD check MBatchUtils



echo "done"
