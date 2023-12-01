## Load/Install packeages -----
library(tidyverse)
install.packages("readxl")
library(readxl)
install.packages('edgeR', 'matrixStats', 'cowplot')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)
library(matrixStats)

## Read in excel file and tidy -----
df <- read_excel("/Users/evanvallenari/Proteomics/proteinGroups_Olav_LMD_698-707-original.xlsx")
df.select <- df %>% 
  select(`Protein names`, `Gene names`, `Fasta headers`, contains('iBAQ'))
#df.select.unite <- unite(df.select, 'iBAQ 1500', 'iBAQ 1507', 'iBAQ 1508', 'iBAQ 1509', 'iBAQ 1510', 'iBAQ 1511', 'iBAQ 1512', 'iBAQ 1513')
df.select.pivot <- df.select %>% 
  pivot_longer(
    cols = contains('iBAQ'),
    names_to = 'sampleID',
    values_to = 'iBAQ'
  )
df.select.pivot.grouped <- df.select.pivot %>% 
  group_by(sampleID) %>% 
  arrange(.by_group = TRUE) %>% 
  select(c(4,1,2,3,5))
?write_csv()
write_csv(df.select.pivot.grouped, 'ProtData_tidy.csv')

df.select.pivot.grouped

## PCA Analysis -----
install.packages("DT", "plotly", "gt")
library(tidyverse)
library(DT)
library(plotly)
library(gt)

df.select
df.select.matrix <- df.select %>% #Create data matrix from data frame
  select(contains('iBAQ 1')) %>% #Select columns containing iBAQ values
  data.matrix() #As data matrix
  row.names(df.select.matrix) <- df.select$`Protein names` #w/ Protein names as row names

#Read in target data -- identify variables of interest
targets <- read.delim("targets.txt")
group <- factor(targets$group)

#Hierarchical clustering
distance <- dist(t(df.select.matrix), method = "euclidean")
clusters <- hclust(distance, method = "single")
plot(clusters)
?dist

#PCA
df.select.NAomit.matrix <- na.omit(df.select.matrix) #prcomp cannot run with NA variables. na.omit omits all NAs
pca.res <- prcomp(t(df.select.NAomit.matrix), scale. = F, retx=T)
summary(pca.res) #Prints variance summary 
screeplot(pca.res) #Screeplot is standard way to view eigenvalues for each PCA, not great for presenting
pc.var <- pca.res$sdev^2 #Captures eigenvalues from PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) #Uses eigenvalues to calculate percentage variance accounted for by each PC
pc.per
pca.res.df <- as_tibble(pca.res$x)

sampleLabels <- as.character(targets$sample)
#Plotting with ggplot2 -- 2D XY plot
ggplot(pca.res.df) + #Identify dataframe
  aes(x=PC1, y=PC2, label=sampleLabels, color=group) + #Aesthetics, X, Y, labels, color by group
  geom_point(size=4) +
  #geom_label(nudge_x = 20, nudge_y = 20) + #Haven't figured this out
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + #Create x axis label by pasting together without spaces
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

#Small multiples PCA plot
pca.res.df4 <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group) 
  pca.res.df4.rename <- rename(pca.res.df4, "PC1 (70.1%)" = "PC1", "PC2 (17.1%)" = "PC2", "PC3 (7.9%)" = "PC3", "PC4 (4.3%)" = "PC4") # Rename columns 'new = old' format
?rename
pca.pivot <- pivot_longer(pca.res.df4.rename, #Identify datafram to be pivoted
                          cols = 'PC1 (70.1%)':'PC4 (4.3%)', #Column names to be stored as single variable
                          names_to = "PC", #Column name of that variable
                          values_to = "loadings") #Name of new column storing data values
#labelsSM.PCA <- as.character(c(paste0("PC1 (",pc.per[1],"%",")"), paste0("PC2 (",pc.per[2],"%",")"), paste0("PC3 (",pc.per[3],"%",")"), paste0("PC4 (",pc.per[4],"%",")")))


ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) +
  geom_bar(stat='identity') +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()


