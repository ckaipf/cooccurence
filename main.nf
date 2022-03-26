#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.files = [
"1.gff",
"2.gff",
"3.gff"
]

params.default_bedtools_parameters = ""
params.default_distance = 50
params.config = "example.config"

workflow {
  main:
  pairedIntersections(params.files, params.config)
}

workflow pairedIntersections {
  take: filePaths
        config        
  main:
  files = Channel.fromPath( filePaths ) | \
  map( it -> [it.getSimpleName(), it]) 


  parameters = Channel.fromPath(config).splitCsv().map { row -> [ "${row[1]}", "${row[2]}", "${row[3]}"] }
  distances = Channel.fromPath(config).splitCsv().map { row -> [ "${row[1]}", "${row[2]}", "${row[4]}"] }
  
  acc = []
  new File(config).splitEachLine(",") {f -> acc += [f[0]]}
  assert acc.every {it == acc[0]} : "Run identifiers have to be equal for all entries in the config file"
  params.tag = acc.get(0)

  
  betool_parameters = parameters
  max_distance = distances


  files.combine(files) | \
  filter(it -> it[0] != it[2]) | \
  map {it -> [it[0], it[2], it[1], it[3]]} |  \
  join(betool_parameters, by: [0,1], remainder: true) |  \
  map{it -> (it[4] == null) ? it[0..3] + params.default_bedtools_parameters : it} | \
  closestBed | \
  join(max_distance, by: [0,1], remainder: true) |  \
  map{it -> (it[3] == null) ? it[0..2] + params.default_distance : it} | \
  map{it -> ['"'+it[0]+'"', '"'+it[1]+'"', '"'+it[3]+'"', '"'+it[2]+'"']} | \
  toList | \
  fullJoin | \
  plotVenn
}



process closestBed {
  input:
  tuple val(i), val(j), path(a), path(b), val(parameters)
  output:
  tuple val(i), val(j), path("*.closest")

"""
sort -k4,4n ${a} > ${a}.sorted
sort -k4,4n ${b} > ${b}.sorted
bedtools closest -D a ${parameters} -a ${a}.sorted -b ${b}.sorted > ${i}${j}.closest
"""
}

process fullJoin {
  publishDir ".", mode: "copy"
  input:
  val closestOutput
  output:
  path "*.RData"
  
  script:
  //  Groovy list to R list string
  input_list = closestOutput
      .collect{f -> f.toString().replace("[", "c(").replace("]", ")")}
      .toString()
      .minus("[")
      .minus("]")

  """
  #!/usr/bin/Rscript
  library(tidyverse)

  x <- list(${input_list}) 
  permutations <- map(x, ~ paste0(.[1],.[2])) %>% unlist
  sets <- map(x, ~ .[1]) %>% unlist(use.names = F) %>% unique
  names(x) <- permutations
  names(permutations) <- permutations

  pairs <- x %>%
  map(function(y) {
    read_delim(file = y[[4]], delim = "\t", col_names = F, col_types = "ccccccccccccccccccc") %>%
    mutate(id_a = paste0(y[[1]], X4, X5, X7)) %>%
    mutate(id_b = paste0(y[[2]], X13, X14, X16)) %>%
    rename_with(.fn = function(x) y[[1]], .cols =  c("id_a")) %>%
    rename_with(.fn = function(x) y[[2]], .cols =  c("id_b")) %>%
    mutate(distance = abs(as.numeric(X19))) %>%
    select(!starts_with("X")) %>%
    distinct(across(everything()))
    }
  )

  pairs_filtered <- 
    map2(.x = pairs, .y = x, .f = ~ filter(.x, distance < as.numeric(.y[[3]]))) %>%
    map(~ select(., -distance))

  sets_total <- sets %>%
    map(~ permutations[str_starts(permutations, .)]) %>%
    map(~ first(.)) %>%
    map(~ map(., .f = function(x) pairs[[x]])) %>% 
    flatten %>%
    map(~ .[1]) 

  joined_filtered <- sets %>%
    map(~ permutations[str_starts(permutations, .)]) %>%
    map(~ map(., .f = function(x) pairs_filtered[[x]])) %>%
    map(~ reduce(., function(acc, y) full_join(acc, y, na_matches = "na"))) %>%
    reduce(function(acc, y) full_join(acc, y, na_matches = "na"))

  venn_complete <- sets_total %>%
    reduce(.f = function(acc, y) full_join(acc, y, na_matches = "na"), .init = joined_filtered) %>%
    rowid_to_column()
    
  save(venn_complete, sets, file = "${params.tag}.RData")
  """
}

process plotVenn{
  publishDir ".", mode: "copy"

  input:
  path Rdata
  
  output:
  path "*.png"


  """
  #!/usr/bin/Rscript

  library(tidyverse)
  load("${Rdata}")
  
  names(sets) <- sets
  ids <- sets %>%
    map(~ filter(venn_complete, !is.na(.data[[.]]))) %>%
    map(~ select(., rowid) %>% unlist %>% unique) 

  p <- ggVennDiagram::ggVennDiagram(ids)
  ggsave("${params.tag}_venn.png", plot = p, device = "png")
  """
}