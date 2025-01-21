# mfuzz: using noise-robust soft clustering of omics time-series data
# in our case, this package allowed us to visualize how different cell 
# surface markers change in their expression over time; this is really
# helpful since we have very distinct timepoints.

# we chose to use this package to expand on our previous k-means clustering. 
# this is useful to illustrate exactly what differentiates each cluster from each other, and 
# helps identify future patterns for us to build on

# col1 -> avg of myo data, col2 -> avg of d10, col3 -> avg of d15, col4 -< avg of d20

# library loading! woohoo!
library(Mfuzz) # soft clustering of omics time-series data
library(dplyr) # data manipulation
library(tidyr) # tidying data + more manipulation

# loading in data
mfuzz_mat <- read.csv("/Users/bellapfeiffer/Dropbox (Harvard University)/Circle Tx/R Scripts/LegendScreen Analysis BP/timepoint analysis.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
# converting to expressionset data type
eset <- ExpressionSet(assayData = as.matrix(mfuzz_mat))

# conduct mfuzz analysis
c <- mfuzz(eset, c = 8, m = 1.25)  # c = number of clusters, m = fuzziness parameter
# after more research: no need to do hard k-means; does not assist in comparison

# output which markers are in which cluster as .txt file (for cross reference with PCA)
membership_values <- c$membership #only include membership values
max_membership <- apply(membership_values, 1, which.max) 
write.table(max_membership, file = "mfuzz.cluster.membership.v1.txt") # save which 

# visualize as clusters - set ylim to be 100
mfuzz.plot2(eset, c, ylim=c(0,100), mfrow = c(2, 4), new.window = FALSE)


