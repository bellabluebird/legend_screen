#loading libraries
library(dplyr) # data manipulation
library(tidyr) # tidying data + more manipulation
library(ggplot2) # visualization
library(factoextra) # basic PCA + clustering visualization
library(FactoMineR) # multivariate analysis 
library(openxlsx) # input / output of excel files for ease of use by biologists

data <- read.csv("/Users/bellapfeiffer/Dropbox (Harvard University)/Circle Tx/Raw Data/Aim 1/Aim 1 Cell Characterization/LegendScreen/Data Analysis/LS_240624.csv", header = TRUE, stringsAsFactors = FALSE)

data$Passage<-gsub("^$", "SkMO", data$Passage)
data<-data %>% mutate(across(where(is.character), ~ na_if(.,"")))
data <- data %>% drop_na(Type)
data<-data %>% mutate(ID=paste0(Type, "_", Line, "_", Passage, "_", Day))

mat<- data %>% 
select(Target, ID, pct_pop) %>% 
pivot_wider(names_from=ID, values_from=pct_pop,values_fn = mean) %>%
column_to_rownames("Target") %>%
as.matrix()

COL<-data.frame(COl2=colnames(mat))
hab<-COL %>% separate(COl2, c("Type", "Line", "Passage", "Day"))

res.pca <- PCA(t(mat),  graph = FALSE)
fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("yellow", "green4"),
             repel = TRUE
             )

#modify the hab for Line/Type/Passage and give a color vector 
g_pca_ind<-fviz_pca_ind(res.pca,
             geom="point",
             #col.ind="black",
             label = "all", # hide individual labels
             habillage = as.factor(hab$Type), # color by groups
             repel=TRUE,
             palette = c("yellow", "green4"),
             addEllipses = FALSE # Concentration ellipses
             )     

ggsave("fviz_pca_ind_type.v1.pdf", g_pca_ind)        
get_eig(res.pca)
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
var <- get_pca_var(res.pca)
fviz_pca_var(res.pca, col.var = "black")
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #repel = TRUE # Avoid text overlapping
             )

fviz_pca_var(res.pca, geom = "point", col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #repel = TRUE # Avoid text overlapping
)

fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
ind <- get_pca_ind(res.pca)

#k-means clustering - PC1 vs PC2
new_mat <- na.omit(mat)
df <- scale(new_mat)
fviz_nbclust(new_mat, kmeans, method = "gap_stat")

# 2. compute k-means + cluster
set.seed(126)
mat.res <- kmeans(scale(new_mat), 8, nstart = 25)
write.xlsx(mat.res$cluster, "matrix.v1.xlsx")
print(x)

# 3. visualize
# as requested: no labels, larger dots
fviz_cluster(mat.res, data = df, pointsize = 1, labelsize = 0,
             palette = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF", "#FFA500", "#800080"),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot v1"
)

# as requested: 
fviz_cluster(mat.res, data = df, pointsize = 1, labelsize = 7,
             palette = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF", "#FFA500", "#800080"),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot v2"
)
