#!/bin/bash

echo "move to app dir"
cd ../apps

echo "build MBatchUtils"
R CMD build MBatchUtils

echo "done"
