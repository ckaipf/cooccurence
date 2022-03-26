# pairedIntersections

## Features

* Simple `nextflow` and `R` script to calculate the intersections of `gff` files
* Intersections are determined pairwise using `bedtools closest`
    * This allows to set detailed constraints for the pairs, e.g. a feature has to occure upstream of another
* A use case is the analysis of the transcriptional structure
    * For example the joint occurence of predicted promoters, annotated genes and terminators
 * A dataframe of all joined features is exported as an `R` object for further analysis
    * Moreover, the intersections are plotted as a Venn diagram

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
 * A dataframe of all joined features is exported as an `R` object
 * Symmetry is not validated
