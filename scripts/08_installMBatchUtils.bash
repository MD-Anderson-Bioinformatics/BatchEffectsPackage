#!/bin/bash

echo "move to app dir"
cd ../apps

echo "INSTALL MBatchUtils"
R CMD INSTALL MBatchUtils

echo "done"
