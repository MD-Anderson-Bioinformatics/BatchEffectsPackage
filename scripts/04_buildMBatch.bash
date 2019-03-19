#!/bin/bash

echo "move to app dir"
cd ../apps

echo "build MBatch"
R CMD build MBatch

echo "done"
