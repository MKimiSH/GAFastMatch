# GAFastMatch
Implemetation of BMVC15 paper: [Fast Affine Template Matching over Galois Field](http://www.bmva.org/bmvc/2015/papers/paper121/index.html)

Based on [FastMatch](http://www.eng.tau.ac.il/~simonk/FastMatch/)

## Introduction
This template matching algorithm uses genetic algorithm and discretization of search field to reduce redundant computation while keeping the good configurations within the search list. 

*Note*: It is reported in the paper that the algorithm achieves comparable or better accuracy than the FastMatch algorithm, but I did not manage to reproduce the result even after checking with the author couple of times. I am sure there is something tricky about this type of heuristics. 

## Usage
run `GA_FastMatch_demo.m`

