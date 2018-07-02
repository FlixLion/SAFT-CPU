#!/bin/sh

mkdir data_ref
cp -r ~/code/ipe/test/SAFT-CPU/asm_src/addsig2vol_3.0_MT+64_SSE1_\(Double\)/data/* ./data_ref 

cmp --silent ./data/ ./data_ref/

echo "COPY DONE"
