#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow {
  main:
  cooccurrence(params.files, params.config)
}

workflow plot {
  main:
    cooccurrence(params.files, params.config)
    cooccurrence.out.order_csv | plotOrder
}


workflow cooccurrence {
  take: 
    files
    config
  main:
  files_ch = Channel.fromPath(files) | \
          sortGff | \
          map { it -> [it.getSimpleName(), it] }
  
  config = Channel.fromPath(config) | \
      splitCsv() | \
      map { row -> (row[0] == params.tag) ? row : null } |  \
      map { row -> (files.any { it.contains(row[1] + ".gff")} && files.any { it.contains(row[2] + ".gff")}) ? row : null} 

  parameters = config | \
    map { it -> [ it[1], it[2], it[3] ] }

  distances = config | \
    map { it -> [ it[1], it[2], it[4] ] }

  files_ch.combine(files_ch) | \
    filter(it -> it[0] != it[2]) | \
    map {it -> [it[0], it[2], it[1], it[3]]} |  \
    join(parameters, by: [0,1], remainder: true) |  \
    map { it -> (it[4] == null) ? it[0..3] + params.bedtools_default_parameters : it } |  \
    closestBed | \
    join(distances, by: [0,1], remainder: true) | \
    map { it -> (it[3] == null) ? it[0..2] + params.default_distance : it } | \
    monitorParams

    monitorParams.out.data | \
    rearrange | \
    catFiles | \
    collectFile(name: "collected_file.txt")  | \
    buildCompleteGraphs

    monitorParams.out.warnings.toList().forEach { if(it) log.warn(it.trim()) }
    buildCompleteGraphs.out.warnings.view()

    emit: 
      set_csv = buildCompleteGraphs.out.set_csv
      order_csv = buildCompleteGraphs.out.order_csv
}

workflow collectRuns {
  main:
  Channel.fromPath( params.dir + "/**_set.csv.gz") | \
    zcatFilesAppendFileName | \
    collect | \
    appendFiles | \
    collectedBarPlot
}

/*
    k is not large enough if
      exists max(observed_distances) < parameter(distance)
*/
process monitorParams {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'
  errorStrategy "ignore"
  input:
  tuple val(i), val(j), path(a), val(parameters)
  output:
  tuple val(i), val(j), path(a), val(parameters), emit: data
  stdout emit: warnings
  
  """
  sort -k1,1 -k4,4n -k5,5n -k7,7 ${a} | \
  awk 'function abs(v) {return v < 0 ? -v : v}
  { 
    for(i=1;i<=NF-1;i++) printf \$i"\t"; print abs(\$19)
  }' | \
  bedtools groupby -g 1,4,5,7 -c 19 -o max,count | \
  awk '{
    if(\$5<${parameters} && \$6==${params.bedtools_k})
      {print "Pairs may have been skipped for \\033[4m${i}\\033[24m and \\033[4m${j}\\033[24m. Consider to increase parameter k.";exit 0}
  }' 
  """
}

process catFiles {
  container 'debian:stable'
 // publishDir ".", mode: "copy"
  input:
  tuple val(i), val(j), path(csv), val(parameters)
  output:
  path "*.gz"
 """
 awk 'BEGIN{OFS=","}{print "${i}","${j}","${parameters}",\$0}' ${csv} | \
 gzip \
 > ${csv}.gz
 """
}

process sortGff {
  container 'debian:stable'

  input:
  path(a)
  output:
  path("${a}.sorted")

  """
  sort -k1,1 -k4,4n -k5,5n ${a} > ${a}.sorted
  """
}

process closestBed {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'
  input:
  tuple val(i), val(j), path(a), path(b), val(parameters)
  output:
  tuple val(i), val(j), path("*.closest")

  """
  bedtools closest -D a ${parameters} -k ${params.bedtools_k}  ${params.bedtools_global} -a ${a} -b ${b} > ${i}${j}.closest
  """
}

process rearrange {
  container 'debian:stable-slim'
  input:
  tuple val(i), val(j), path(file), val(parameters)
  output:
  tuple val(i), val(j), path("*.csv"), val(parameters)

"""
awk  'BEGIN{OFS=","} {print "${i}_"\$4"_"\$5"_"\$7, "${j}_"\$13"_"\$14"_"\$16, \$19}' ${file} > ${file}.csv
"""
}


process buildCompleteGraphs {
  container 'python:slim'
  publishDir "${params.dir}/${params.tag}", mode: "copy"
  input:
  path csv
  output:
  path "*_set.csv.gz", emit: set_csv
  path "*_order.csv.gz", emit: order_csv
    stdout emit: warnings

  script:
  """
  build_complete_graphs.py ${csv} ${params.tag}
  """
}

process plotOrder {
  container 'rocker/tidyverse:latest'
  publishDir "${params.dir}/${params.tag}", mode: "copy"
  input:
  path csv
  output:
  path "*.png"
 
  """
  plot_order.R ${csv} ${params.tag}
  """
}

process zcatFilesAppendFileName {
//  publishDir "${params.dir}", mode: "copy"
  input:
  path f
  output:
  path "*t.csv"
 """
 zcat ${f} | \
 awk -F ',' -v OFS=',' '{if(NR==1) {print "tag",\$0} else {print "${f.simpleName}",\$0}}' \
 > ${f.simpleName}.t.csv
 """
}

process appendFiles {
  publishDir "${params.dir}", mode: "copy"
  input:
  path f
  output:
  path "collected.csv"
 """
 head -1 ${f[0]} > collected.csv
 tail -q -n +2 ${f} >> collected.csv
 """
}

process collectedBarPlot {
  container 'rocker/tidyverse:latest'
  publishDir "${params.dir}/", mode: "copy"
  input:
  path csv
  output:
  path "*.png"

  """
  collected_bar_plot.R ${csv}
  """
}