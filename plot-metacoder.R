library('metacoder')
library('dplyr')

# load the output of 'classify-genome-sigs.py'
csvdata = read.csv('delmont-genome-sigs.taxonomy.csv', header=TRUE)

# eliminate empty lineages
csvdata <- filter(csvdata, lineage != "")

# createa mashup of name and lineage to satify extract_taxonomy
csvnames <- paste(csvdata$name, csvdata$lineage, sep='\t')

# heck if I know why I have to do this :), but set list names to list values
names(csvnames) <- csvnames

taxdata <- extract_taxonomy(csvnames, regex = "^(.*)\\\t(.*)",
                         key = c(id = "obs_info", "class"), class_sep=';')

taxon_data(taxdata)
heat_tree(taxdata, node_size = n_obs, node_label = name, node_color = n_obs)