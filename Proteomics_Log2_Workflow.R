##Install and load packages -----
install.packages("DT", "plotly", "gt")
library(tidyverse)
library(DT)
library(plotly)
library(gt)
library(readxl)

##Read in and prep data -----
#Read in excel file and select columns
df <- read_excel("/Users/evanvallenari/Proteomics/proteinGroups_Olav_LMD_698-707-original.xlsx")
df.select <- df %>% 
  select(`Protein.names`, `Gene.names`, `Fasta.headers`, contains('LFQ'))

df.select2 <- df %>%
  select(Protein.names, contains('LFQ')) 

as_tibble(df.select2)
  
#Convert to matrix
df.select.matrix <- df.select %>% #Create data matrix from data frame
  select(contains('LFQ.')) %>% #Select columns containing LFQ values
  data.matrix() #As data matrix
row.names(df.select.matrix) <- df.select$Protein.names #w/ Protein names as row names

#Read in target data -- identify variables of interest
targets <- read.delim("targets.txt")
group <- factor(targets$group)

#Log2 transformation of data matrix containing LFQ intensity
log2.data.matrix <- log2(df.select.matrix)

#log2.df <- df.select %>% #
  select(contains('LFQ.')) %>% #Select columns containing LFQ values
  log2()
  .rowNamesDF(log2.df, make.names = TRUE) <- df.select$Protein.names
?row.names.data.frame

log2.df <- df.select %>%
  select(contains('LFQ.')) %>%
  log2() %>%
  mutate(Protein.names = df.select$Protein.names) %>%
  as_tibble() %>%
  relocate(Protein.names)
#Remove infinite values (0.0 LFQ)
log2.data.matrix[is.infinite(log2.data.matrix)] <- NA #Convert infinite values to NA
log2.data.matrix.NAomit <- na.omit(log2.data.matrix) #Omit NA values

#Remove infinite values from dataframe
is.na(log2.df) <- sapply(log2.df, is.infinite)
log2.df.NAomit <- na.omit(log2.df)

##PCA -----
# Hierarchical clustering w/ log2 normalized data
distance <- dist(t(log2.data.matrix.NAomit), method = "euclidean")
clusters <- hclust(distance, method = "single")
plot(clusters)

#PCA
pca.res <- prcomp(t(log2.data.matrix.NAomit), scale. = F, retx=T)
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

ggsave("Log2_PCA_1v2.png")

#Small multiples PCA plot
pca.res.df4 <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group) 
pca.res.df4.rename <- rename(pca.res.df4, "PC1 (79.0%)" = "PC1", "PC2 (7.9%)" = "PC2", "PC3 (5.1%)" = "PC3", "PC4 (2.7%)" = "PC4") # Rename columns 'new = old' format
?rename
pca.pivot <- pivot_longer(pca.res.df4.rename, #Identify datafram to be pivoted
                          cols = 'PC1 (79.0%)':'PC4 (2.7%)', #Column names to be stored as single variable
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

ggsave("Log2_PCA_SmallMultiples.png")

## Calculate fold change -----
mydata.df <- log2.df.NAomit %>% #Create mydata.df with new columns for averages and fold change
  mutate(cancer.far.AVG = LFQ.intensity.1500 + LFQ.intensity.1512 / 2,
         cancer.near.AVG = LFQ.intensity.1507,
         muscle.far.AVG = LFQ.intensity.1508 + LFQ.intensity.1510 / 2,
         muscle.near.AVG = LFQ.intensity.1509 + LFQ.intensity.1511 + LFQ.intensity.1513 / 3,
         LogFC.cancer = cancer.far.AVG - cancer.near.AVG,
         LogFC.muscle = muscle.far.AVG - muscle.near.AVG) %>%
  mutate_if(is.numeric, round, 2)

#Separate out cancer and muscle
mydata.df.cancer <- mydata.df %>%
  dplyr::select(Protein.names, cancer.far.AVG, cancer.near.AVG, LogFC.cancer) %>% #Pull out two columns
  dplyr::arrange(desc(LogFC.cancer)) #Arrange based on descending order of fold change

mydata.df.muscle <- mydata.df %>%
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

#Make plot interactive