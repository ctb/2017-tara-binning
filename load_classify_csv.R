# load the output of 'classify-genome-sigs.py'
library('metacoder')
library('dplyr')

load_classify_csv <- function(filename) {
  csvdata = read.csv(filename, header=TRUE)

  # eliminate empty lineages
  csvdata <- filter(csvdata, lineage != "")

  # create a mashup of name and lineage to satisfy extract_taxonomy
  csvnames <- paste(csvdata$name, csvdata$lineage, sep='\t')

  # heck if I know why I have to do this :), but set list names to list values
  names(csvnames) <- csvnames

  taxdata <- extract_taxonomy(csvnames, regex = "^(.*)\\\t(.*)",
                              key = c(id = "obs_info", "class"), class_sep=';')
  
  taxdata
}

