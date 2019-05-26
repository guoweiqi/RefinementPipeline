#!/bin/sh
tar -czvf pdbFiles.tar.gz pdbFiles/
scp pdbFiles.tar.gz gqi1@argon.hpc.uiowa.edu:/Dedicated/schnieders/gqi
rm pdbFiles.tar.gz
