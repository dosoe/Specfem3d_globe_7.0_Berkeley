#!/bin/bash

# compute total number of nodes needed
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -c 34- `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -c 34- `
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -c 34- `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))
echo "Num of processes: $numnodes"