bib <- readLines('../../../articles.bib')
bib <- bib[!grepl('eprint', bib)]
url <- grepl("url\\s*=", bib)
rproj <- grepl("R-project\\.org", bib)
bib <- bib[!url | rproj]
writeLines(bib, 'articles.bib')
