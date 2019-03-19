#!/bin/bash

echo "move to app dir"
cd ../apps

echo "INSTALL MBatch"
R CMD INSTALL MBatch

echo "done"
