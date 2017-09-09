
args <- commandArgs(TRUE)

source('load_classify_csv.R')

# load the specified csv file
taxdata <- load_classify_csv(args[1])

# open a PDF file -> output graph
pdf(args[2], onefile=FALSE)

# plot a heat tree
heat_tree(taxdata, node_size = n_obs, node_label = name, node_color = n_obs)
dev.off()

