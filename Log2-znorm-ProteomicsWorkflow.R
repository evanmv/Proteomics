##Load packages -----
library(tidyverse)
install.packages("BBmisc")
library(BBmisc)
library(DT)
library(plotly)
library(gt)
library(readxl)
library(limma)
library(edgeR)

##Pre-processing -----

#Filter out peptides only IDd by site, reverse, and potential contaminants
df.filtered.tv <- filter(df, Only.identified.by.site != "+" & 
                           Reverse != "+" &
                           Potential.contaminant != "+")
df.filtered.select <- df.filtered.tv %>%
  select(Protein.names, contains('LFQ'))

#Transform data (rows as columns) and keep rows numeric. Apply prot names as column names and remove first row (was column names, not numerical, so we have to remove)
df.t <- t(sapply(df.filtered.select, as.numeric))

colnames(df.t) <- df.filtered.select$Protein.names
df.t <- df.t[-1,]

df.t2 <- log2(df.t)

#Uses function scale to Z-score normalize within protein column, confirm by taking mean/sd
df.t.scale <- scale(df.t2) %>%
  as.data.frame()
mean(df.t.scale$V14) #Confirm z-score normalization - mean = 0, st.dev = 1
sd(df.t.scale$V14)

#Transform (Each column is a sample, each row is a protein)
df.revert <- t(df.t.scale)

#Omit rows with NA 
is.na(df.t.scale) <- sapply(df.t.scale, is.infinite)
df.scale.NAomit <- na.omit(df.revert)


##PCA -----
#Heirarchical clustering
distance <- dist(t(df.scale.NAomit), method = "euclidean")
clusters <- hclust(distance, method = "single")
plot(clusters)

#PCA
pca.res <- prcomp(t(df.scale.NAomit), scale. = F, retx=T)
summary(pca.res) #Prints variance summary 
screeplot(pca.res) #Screeplot is standard way to view eigenvalues for each PCA, not great for presenting
pc.var <- pca.res$sdev^2 #Captures eigenvalues from PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) #Uses eigenvalues to calculate percentage variance accounted for by each PC
pc.per
pca.res.df <- as_tibble(pca.res$x)


sampleLabels <- as.character(targets$sample) #Creates character string of sample labels
#Plotting with ggplot2 -- 2D XY plot
ggplot(pca.res.df) + #Identify dataframe
  aes(x=PC1, y=PC2, label=sampleLabels, color=group) + #Aesthetics, X, Y, labels, color by group
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + #Create x axis label by pasting together without spaces
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggsave("Scaled_PCA_1v2.png")


Small multiples PCA plot
pca.res.df2 <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group) 
pca.res.df2.rename <- rename(pca.res.df2, "PC1 (57.2%)" = "PC1", "PC2 (13.4%)" = "PC2", "PC3 (10.2%)" = "PC3", "PC4 (8.4%)" = "PC4") # Rename columns 'new = old' format
?rename
pca.pivot <- pivot_longer(pca.res.df2.rename, #Identify dataframe to be pivoted
                          cols = 'PC1 (57.2%)':'PC4 (8.4%)', #Column names to be stored as single variable
                          names_to = "PC", #Column name of that variable
                          values_to = "loadings") #Name of new column storing data values

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) +
  geom_bar(stat='identity') +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

ggsave("Scaled_PCA_SmallMultiples.png")


## Calculate fold change -----

mydata.scaled.df <- df.scale.NAomit %>% #Create mydata.df with new columns for averages and fold change
  as_tibble() %>%
  mutate(Protein.names = row.names(df.scale.NAomit)) %>%
  mutate(cancer.far.AVG = LFQ.intensity.1500 + LFQ.intensity.1512 / 2,
         cancer.near.AVG = LFQ.intensity.1507,
         muscle.far.AVG = LFQ.intensity.1508 + LFQ.intensity.1510 / 2,
         muscle.near.AVG = LFQ.intensity.1509 + LFQ.intensity.1511 + LFQ.intensity.1513 / 3,
         LogFC.cancer = cancer.far.AVG - cancer.near.AVG,
         LogFC.muscle = muscle.far.AVG - muscle.near.AVG) %>%
    relocate(Protein.names) %>%
  mutate_if(is.numeric, round, 2) 

#Separate out cancer and muscle
mydata.df.cancer <- mydata.scaled.df %>%
  dplyr::select(Protein.names, cancer.far.AVG, cancer.near.AVG, LogFC.cancer) %>% #Pull out two columns
  dplyr::arrange(desc(LogFC.cancer)) #Arrange based on descending order of fold change

mydata.df.muscle <- mydata.scaled.df %>%
  dplyr::select(Protein.names, muscle.far.AVG, muscle.near.AVG, LogFC.muscle) %>%
  dplyr::arrange(desc(LogFC.muscle))

#Static table
gt(mydata.df.cancer)

#Create interactive table
datatable(mydata.df.cancer,
          extensions = c('KeyTable', 'FixedHeader'),
          filter = 'top',
          options = list(keys = TRUE,
                         searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))

datatable(mydata.df.muscle,
          extensions = c('KeyTable', 'FixedHeader'),
          filter = 'top',
          options = list(keys = TRUE,
                         searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))

#Create scatter plot of average values
canplot <- ggplot(mydata.df.cancer) +
  aes(x=cancer.far.AVG, y=cancer.near.AVG,
      text = paste(mydata.df.cancer$Protein.names)) + #Adds protein names for mouseover tooltip in ggplotly
  geom_point(shape=16, size=1) +
  ggtitle("near vs. far") +
  theme_bw()

ggplotly(canplot)

muscleplot <- ggplot(mydata.df.muscle) +
aes(x=muscle.far.AVG , y=muscle.near.AVG,
    text = paste(mydata.df.muscle$Protein.names)) + #Adds protein names for mouseover tooltip in ggplotly
  geom_point(shape=16, size=1) +
  ggtitle("near vs. far") +
  theme_bw()

ggplotly(muscleplot)

##Volcano plot -----

group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.cancer <- voom(mydata.df.cancer, design, plot = TRUE)
