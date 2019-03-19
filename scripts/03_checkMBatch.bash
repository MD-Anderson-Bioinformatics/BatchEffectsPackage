#!/bin/bash

echo "move to app dir"
cd ../apps

echo "check MBatch"
env _R_CHECK_FORCE_SUGGESTS_=0 R CMD check MBatch



echo "done"
