# Bivex

Bivex is a stochastic simulator for bivalent chromatin mark patterns (H3K4me3 and H3K27me3) and gene expression. A string of nucleosomes ("chromtape", mimicking promoter region of a gene) undergoes iterations where the chromatin mark pattern gets updated by histone modifiers ("rules"). At each iteration, the chromatin mark pattern is considered in calling the gene expression status (ON/OFF).

The basic unit of operation in Bivex simulates the process of gene regulation by chromatin marks in a single cell. By running it over many cells and by looking at the overall result, Bivex also simulates population level histone mark patterns and gene expression outcome (i.e. typical data obtained from ChIP-seq and RNA-seq).

## Requirement

Bivex is built on Clojure 1.8.0., Leiningen 2.7.1, and Java 1.8.
Tested Mac OS X and Ubuntu.

## Usage

In a single-run
Command line argument processing not yet implemented;

For a quick look, `lein run`

## License
Not yet. Eventually under an open source 

Copyright © 2018 Ah-Jung Jeon


