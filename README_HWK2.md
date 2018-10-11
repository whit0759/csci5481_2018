# CSCI 5481: Homework 2
## Prof: Knights

## Christopher White
## October 11, 2018

## Introduction
This file processes two sequences for alignment using the Needleman-Wunsch 
algorithm. The algorithm has fix values for matches (1), mismatches (-3), and 
gaps (-2). The algorithm assumes that the head and tail of the sequences should 
be penalized for gaps.

In the case that sequence anchors are given, the match file will allow the 
algorithm to run over the set of sequences between the matches.

The output file reports the aggregate score for the alignment and the full 
aligned sequences for each species. If anchors are used, the head and tail 
portions of the sequences are appended but not aligned to each other.

## Program Flags
This program is designed to run on the command line with command flags.
    -q, --query:    Assumed to refer to the 'human' sequence files
    -r, --ref:      Assumed to refer to the 'fly' sequence files
    -o, --output:   Optional output file
    -m, --match:    Optional match file with sequence anchors
    
## Usage

```bash
python nwalign.py -q <HUMAN> -r <FLY> [-o <OUTPUT FILE>]  [-m <MATCH FILE>]
```

This file can also be imported to use the functions for processing of sequences in other use cases.

