#!/bin/bash

echo "move to app dir"
cd ../apps

echo "copy MBatch vignettes PDFs"
cp MBatch.Rcheck/MBatch/doc/*.pdf ../docs/MBatch/.

#echo "copy MBatchUtils vignettes PDFs"
#cp MBatchUtils.Rcheck/MBatchUtils/doc/*.pdf ../docs/MBatchUtils/.

echo "list MBatch docs"
ls -l ../docs/MBatch

echo "list MBatchUtils docs"
ls -l ../docs/MBatchUtils

echo "done"
