# pairedIntersections

## Features

* Simple `nextflow` and `R` script to calculate the *intersections* (paired joint occurrences) of `gff` files
* Intersections are determined pairwise using `bedtools closest`
    * This allows to set detailed constraints for the pairs, e.g. a feature has to occur upstream of another
* A use case is the analysis of the transcriptional structure
    * For example the joint occurrence of predicted promoters, annotated genes and terminators
 * A dataframe of all joined features is exported as an `R` object for further analysis
    * Moreover, the intersections are plotted as a Venn diagram
* If the `gffs` are dense and multiple joint calls are possible, consider the usage of `bedtools` *k* parameter

## Requirements

* `nextflow >= 21.04.3.5560`
* `bedtools >= 2.27.1`
* `R >= 4.1.2`
* `tidyverse >= 1.3.1`
* `ggVennDiagram >= 1.3.1`

## How to

Run the script:

```
nextflow main.nf
```


## Notes

 * Files have to be in the `gff` format and are set in the script 
 * Parameters for `bedtools closest` and maximal distances are passed in a csv config file
    * Default parameters can be set within the script
    * Only non-default lines have to be set in the config
    * Provide at least one entry in the config
    * Use the simplified file names (`myfile.gff -> myfile`) as identifiers for the combinatorial table
 * Symmetry is not validated
    * The *k* closest function is not commutative, therefore it is necessary to calculate all permutations
 * The visualization as Venn diagrams can be counterintuitive to interpret
    * The sets (the original and non-intersected) reflect single intervals, e.g. *A={a_1, ..., a_n}* and *B={b_1, ..., b_n}*
    * The *intersections* reflect paired joined occurrences of both sets: *intersect(A, B)= {(a_1, b_1), (b_1, a_1), (a_1, b_2)}*
    * The cardinality of the *intersected* sets consists of the original set, e.g. *|A|* and a *product* part of the joined features (*|joined(AxB)|*)
 