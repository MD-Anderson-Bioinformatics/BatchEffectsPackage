#!/bin/bash

# not done during CI
echo "START 09_copyDocs"
set -e
BASE_DIR=$1

echo "BASE_DIR=${BASE_DIR}"

echo "copy MBatch vignettes as PDFs"

cd ${BASE_DIR}/apps/MBatch.Rcheck/00_pkg_src/MBatch/inst/doc
rm -rf ${BASE_DIR}/apps/MBatch/doc
mkdir -p ${BASE_DIR}/apps/MBatch/doc
pwd
wordFiles=(*.html)
echo "${wordFiles[@]}"
for wordFile in "${wordFiles[@]}"; do
	echo "$wordFile"
	baseFile=$(echo "$wordFile" | cut -f 1 -d '.')
	echo "$baseFile"
	pandoc -o "${BASE_DIR}/docs/MBatch/${baseFile}.pdf" --self-contained "${BASE_DIR}/apps/MBatch.Rcheck/00_pkg_src/MBatch/inst/doc/${wordFile}"
	pandoc -o "${BASE_DIR}/apps/MBatch/doc/${baseFile}.pdf" --self-contained "${BASE_DIR}/apps/MBatch.Rcheck/00_pkg_src/MBatch/inst/doc/${wordFile}"
done

echo "copy MBatchUtils vignettes as PDFs"

cd ${BASE_DIR}/apps/MBatchUtils.Rcheck/00_pkg_src/MBatchUtils/inst/doc
rm -rf ${BASE_DIR}/apps/MBatchUtils/doc
mkdir -p ${BASE_DIR}/apps/MBatchUtils/doc
pwd
wordFiles=(*.html)
echo "${wordFiles[@]}"
for wordFile in "${wordFiles[@]}"; do
	echo "$wordFile"
	baseFile=$(echo "$wordFile" | cut -f 1 -d '.')
	echo "$baseFile"
	pandoc -o "${BASE_DIR}/docs/MBatchUtils/${baseFile}.pdf" --self-contained "${BASE_DIR}/apps/MBatchUtils.Rcheck/00_pkg_src/MBatchUtils/inst/doc/${wordFile}"
	pandoc -o "${BASE_DIR}/apps/MBatchUtils/doc/${baseFile}.pdf" --self-contained "${BASE_DIR}/apps/MBatchUtils.Rcheck/00_pkg_src/MBatchUtils/inst/doc/${wordFile}"
done

echo "list MBatch docs"
ls -lh ${BASE_DIR}/docs/MBatch

echo "list MBatchUtils docs"
ls -lh ${BASE_DIR}/docs/MBatchUtils

echo "FINISHED 09_copyDocs"

