library(testthat)
library(alakazam)
library(igraph)

ExampleDb <- file.path("..","data-tests","ExampleDb.gz")
ExampleDb <- readChangeoDb(ExampleDb)

test_check("alakazam")
