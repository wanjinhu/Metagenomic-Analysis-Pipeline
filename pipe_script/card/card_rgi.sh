#! /bin/bash

rgi main --input_sequence $1 --input_type protein --alignment_tool DIAMOND --output_file $2

# demo
# rgi main --input_sequence prot_nonerude.faa --input_type protein --alignment_tool DIAMOND --output_file test