# the clustering + visualization seen below helps illustrate which genes are expressed 
# similarly between different timepoints - in our research, this code helped us to identify targets
# that turned on and off at timepoints of interest for further screening.
# we were also able to identify pathways that may be relevant based on k-means clustering.

# for data privacy purposes, i removed the target names from the list and added id numbers

# an outline for the first few portions was given to me by KH in the Rubin Lab
# i am now comfortable performing this type of data manipulation after this project

# loading libraries
library(dplyr) # data manipulation
library(tidyr) # tidying data + more manipulation
library(tibble) # more tibble manipulation
library(ggplot2) # visualization
library(factoextra) # basic PCA + clustering visualization
library(FactoMineR) # multivariate analysis 
library(openxlsx) # input / output of excel files for ease of use by biologists
library(RColorBrewer) # pretty colors :)

# read data as csv
data <- read.csv("/Users/bellap/legend_screen/inputs/LS_241113.csv", header = TRUE, stringsAsFactors = FALSE)

# data manipulation ----
# replace any empty passage field with the string "SkMO"
data$Passage <- gsub("^$", "SkMO", data$Passage)

# converting all character columns to NA if empty string exists
# then drop rows where 'Type' column = NA
# then create new 'ID' column combining, type, line, passage, and day columns
data <- data%>% 
  mutate(across(where(is.character), ~ na_if(.,""))) %>% 
  drop_na(Type) %>% 
  mutate(ID=paste0(Type, "_", Line, "_", Passage, "_", Day))

# create matrix with target columns as row names + values for ID and pct_pop
mat <- data %>% 
  select(Target, ID, pct_pop) %>% # only including these columns
  pivot_wider(names_from=ID, values_from=pct_pop,values_fn = mean) %>% # reshaping + taking averages for each marker
  column_to_rownames("Target") %>% # set the target column -> new rownames
  as.matrix() # converting to matrix

COL<-data.frame(COl2=colnames(mat)) # creating coldata 
hab<-COL %>% separate(COl2, c("Type", "Line", "Passage", "Day")) # separate COL2 column into different columns

# perform PCA on re-transposed matrix 
res.pca <- PCA(t(mat),  graph = FALSE)

# adjust the hab for Line/Type/Passage and give a color vector according to contribution level
g_pca_ind<-fviz_pca_ind(res.pca,
             geom="point",
             col.ind="black",
             label = "all", # hide individual labels
             habillage = as.factor(hab$Type), # color by groups
             repel=TRUE,
             palette = c("yellow", "green4"),
             addEllipses = FALSE # Concentration ellipses
             )     

# get the eigenvalues of the PCA, helps explain variance
get_eig(res.pca)

# visualize using a screeplot, showing the contributions of the top 10 variables for PC1
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

#visualize the pca 
fviz_pca_var(res.pca, geom = "point", col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #repel = TRUE # Avoid text overlapping
)

fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

# store pca in var
var <- get_pca_var(res.pca)
# get pca individual's data and store
ind <- get_pca_ind(res.pca)

# 1. k-means clustering - PC1 vs PC2
new_mat <- na.omit(mat)
df <- scale(new_mat)
fviz_nbclust(new_mat, kmeans, method = "gap_stat") #find optimal number of clusters (in this case, 8)

# 2. compute k-means + cluster
set.seed(126)
mat.res <- kmeans(scale(new_mat), 8, nstart = 25) #nstart = random initialization
write.xlsx(mat.res$cluster, "matrix.v1.xlsx") # saving output as xslx file for biologists

# 3. visualization outputs
# at the time i made this, i created multiple outputs; i have reduced it to 1 for simplicity

# creating a beautiful plot
plot <- fviz_cluster(mat.res, data = df,
             pointsize = 4,
             labelsize = 9,
             palette = brewer.pal(8, "Dark2"),
             ggtheme = theme_minimal(base_size = 20),
             main = "PCA and K-means Clustering of Cell Surface Marker Expression\nin Human SkMOs Across 20-Day Timepoints",
             alpha = 0.7, 
             ellipse.type = "convex", 
             ellipse.alpha = 0.2
) + 
  theme(
    legend.position = "right", 
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5), 
    axis.title = element_text(size = 10), 
    axis.text = element_text(size = 10) 
  )

print(plot)