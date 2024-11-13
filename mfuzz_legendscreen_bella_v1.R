# library loading! woohoo!
library("Mfuzz")

# i sorted this data in excel, the columns need to be in order of "timepoint" for mfuzz
# col1 -> avg of myo data, col2 -> avg of d10, col3 -> avg of d15, col4 -< avg of d20
mfuzz_mat <- read.csv("/Users/bellapfeiffer/Dropbox (Harvard University)/Circle Tx/R Scripts/LegendScreen Analysis BP/timepoint analysis.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
eset <- ExpressionSet(assayData = as.matrix(mfuzz_mat))

# conduct analysis
c <- mfuzz(eset, c = 8, m = 1.25)  # c = number of clusters, m = fuzziness parameter
# after more research: no need to do hard k-means; does not assist in comparison
membership_values <- c$membership 
max_membership <- apply(membership_values, 1, which.max)
write.table(max_membership, file = "mfuzz.cluster.membership.v1.txt")

# visualize as clusters - set ylim to be 100
mfuzz.plot2(eset, c, ylim=c(0,100), mfrow = c(2, 4), new.window = FALSE)