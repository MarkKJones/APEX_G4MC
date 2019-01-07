#!/bin/bash

export DISPLAY=localhost:0.0
mkdir job_0003
cd job_0003
cp ../../files_in/* .
./G4MC -m 1 beam.mac
ls -ltr
