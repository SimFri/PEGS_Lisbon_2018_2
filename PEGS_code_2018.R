
####################################################
## 2018 PEGS Lisbon 
## Immune repertoires with R
## Simon Friedensohn, Alexander Yermanos, Sai T Reddy
## ETH Zurich 
## simon.friedensohn@bsse.ethz.ch, ayermano@ethz.ch
## sai.reddy@ethz.ch
####################################################

## This file walks you through an introduction in R.
## and analysis of NGS of immune repertoires.

###################  PART I  #######################
####################################################
##### Quick start in R and immune repertoires ######
####################################################

# Script and Console ------------------------------------------------------

## Console, Environment, Code

## short cut to run code from script
## > indicates that code starts immediately after
## incomplete code shown by + in console
## errors/warnings shown in red  e.g. > print()


# Basic arithmetic operations
1+1

4-5

2*8

8/2

3^3

## use hashtag to comment (code will not run)

# 3*3


# Variables ---------------------------------------------------------------

x <- 1

y <- 2

x+y

x*y

x <- 5

x+y 

# Data types

# numeric 

5

5*4

# character

"CARWYGGGAWY"

my_cdr3 <- "CARWYGGGAWY"

# boolean / logical

sky_blue <- TRUE

5 >= 4 # greater or equal to

1 > 2 # greater than

4 == 4 #equal to

1 != 4  # is not equal

5>3 & 10>6 # and : both conditions must be met 

1>3 | 10>6 # or : one condition must be met

#??# Fill in the ? in the following code to return TRUE

(7-5)*3 >= ?

#??# Fill in the ? in the following code to return FALSE
"CARYW" == "?"

#??# Fill in the ? in the following code to return FALSE

(T & F) | (T | (F & '?'))   ## gift from Simon


# Data structures - Array ---------------------------------------------------------


## Vector

# numeric vector with the c function, 
# "combine" or "concatenate"

cdr3_lengths <- c(11, 15, 22, 12, 17, 20, 19, 12, 21, 19)

### Character vector of 10 random sequences
cdr3_sequences <- c('CARIMTERSLKW',
                    'CARFEPHVW',
                    'CARRWHQYELCSFANW',
                    'CARHLMKDPWVCIFRGTESYAW',
                    'CARRQKYIEACTGHNW',
                    'CARGKYNFIMCQEWTDLVSPRHW',
                    'CARRIGPEMCSQDAW',
                    'CARVRWIHPKQFDLGCAW',
                    'CARWNVPGHDAKCFSW',
                    'CARLYEPISDQRHAWGMVTFNCKW')

logical_vector <- c(TRUE,TRUE,FALSE,TRUE)

# Select a subset of a vector:

cdr3_lengths[1]

cdr3_lengths[1:4] 

#??# Access the 3rd CDR3 in "cdr3_sequences"

cdr3_lengths[c(1,4,6)] ## can select items using a vector 

logical_vector[-1] ## can remove specified items in vector

cdr3_lengths==22

cdr3_lengths[cdr3_lengths==22] # all values excluding 22


## combining vectors

c(cdr3_lengths,10) ## append element to last spot

c(cdr3_lengths,logical_vector)


#??# What happens when you combine logical_vector (boolean) with cdr3_sequences (character)
c(logical_vector,cdr3_sequences)

#??# What happens when you combine cdr3_lengths (numeric) with cdr3_sequences (character)
c(cdr3_lengths,cdr3_sequences)


# Logical statements can also be vectorized:

x <- c(4, 5)
y <- c(4, 7)

x == y



## ----Data_structures - Matrix----------------------------------------------

mixed_vector <- 1:12

matrix_example <- matrix(mixed_vector, nrow = 4, ncol = 3)
matrix_example

# To order the filling of the matrix by row, set byrow = TRUE
matrix_by_row <- matrix(mixed_vector, nrow = 3, ncol = 4,byrow = T)
matrix_by_row

## can also transpose
t(matrix_example)


# Check the dimension of your matrix
dim(matrix_example)
dim(matrix_by_row)

# Exercise: Try to select only the first row
# Try to select the second column

## ----Data_structures - Data frame----------------------------------------------

# Transform matrix into dataframe
dataframe_example <- data.frame(1:4,c(5,3,10,4))
dataframe_example

# Set row and column names
colnames(dataframe_example) <- c("CDR Index", "Abundance")
rownames(dataframe_example) <- c("Zurich","Geneva","Basel","Lausanne")
dataframe_example

# Add new column to data frame
dataframe_example$CDR3 <- cdr3_sequences[1:4]
dataframe_example$Vgene <- c("V1-2", "V1-69", "V3-46", "V5-1")
dataframe_example

# Use the dollar-sign operator to extract specific columns
dataframe_example$Vgene

#??# Display the first CDR3 sequence from the data frame


## ----Data_structures - List----------------------------------------------

example_list <- list()
example_list[[1]] <- 1
example_list[[2]] <- c("hat","cat","sat")
example_list[[3]] <- dataframe_example

example_list[[2]]


## Can add names to the list elements
names(example_list) <- c('Number', 'Words', 'Dataframe')

# can index the list by names using the $ 
example_list$...


#??# add a numeric vector to the list



# Functions ---------------------------------------------------------------


# functions have parentheses ()

help()

help("matrix")
help("dim")
#Optional usage: ?mean
?dim


## ----help_in_R

apropos("mean") # Returns the names of all objects in the search

example("mean") # Examples part of R's online help

# Optional usage: ?mean

# Documentation on a topic with name (typically, an R object or a data set) 
# can be displayed by either help("name") or ?name.

?mean # Access the documentation on a topic with name (e.g. "mean")
?plot # Access the documentation on a topic (e.g. "plot")

??mean # Search the Help System 

#??? Which functions have "help" in them?



### Built in functions

class(cdr3_lengths) # data class (numeric)

mean(cdr3_lengths)
min(cdr3_lengths)
max(cdr3_lengths)

unique(cdr3_lengths)

summary(cdr3_lengths) # summary statistics


t.test(cdr3_lengths, cdr3_lengths)

sd(cdr3_lengths)

str(dataframe_example) # Try ?str # Internal structure of an R object


#Subsample a dataset
set.seed(1) # Fix the outcome of the pseudo-random number generator

#How to sample
sample(x = cdr3_sequences, size = 4)

# Repeat an element n times
rep("A",3)

print(cdr3_sequences)

c("hat",0)
# Remove NA 
mean(as.numeric(mixed_vector), na.rm = T)

# "intersect" performs set intersection

set.seed(42)
sample1 <- sample(cdr3_sequences, 5)
sample2 <- sample(cdr3_sequences, 5)

which_lengths_overlap <- intersect(sample1, sample2) 

print(which_lengths_overlap) 

## ----data_structures_matrix----------------------------------------------

matrix_example <- matrix(mixed_vector, nrow = 4, ncol = 3)
matrix_example

# To order the filling of the matrix by row, set byrow = TRUE
matrix_by_row <- matrix(mixed_vector, nrow = 3, ncol = 4,
                        byrow = T)

matrix_by_row

# Check the dimension of your matrix
dim(matrix_example)
dim(matrix_by_row)

# Exercise: Try to select only the first row
# Try to select the second column

## ----data_structures_dataframes------------------------------------------

# Transform matrix into dataframe
dataframe_example <- as.data.frame(matrix_example)

# Set row and column names
colnames(dataframe_example) <- c("Sample", "CDR3", "Abundance")
rownames(dataframe_example) <- c(1:nrow(dataframe_example))

dataframe_example

# Add new column
dataframe_example$Vgene <- c("V1-2", "V1-69", "V3-46", "V5-1")
dataframe_example

# Use the dollar-sign operator to extract specific columns
dataframe_example$...

## ----data_structures_arrays----------------------------------------------
array_example <- array(0, dim = c(2, 3, 2))

## ----data_structures_list------------------------------------------------

list_example <- list(cdr3_lengths, cdr3_sequences, dataframe_example,
                     matrix_example)

list_example

length(list_example)

list_example[[3]]


# Let's add one more element

list_example[[5]] <- array_example

# Lists can also be named

names(list_example) <- c('cdr3_lengths', 'cdr3_sequences', 'dataframe_example',
                         'matrix_example', 'array')


#Exercise: Select an element by its name



# Functions ---------------------------------------------------------------

# functions have parentheses ()

help()

help("mean")


## ----help_in_R

apropos("mean") # Returns the names of all objects in the search

example("mean") # Examples part of R's online help

help("mean") # Require help regarding the function "mean", equivalent to "?mean". 
# Optional usage: ?mean

# Documentation on a topic with name (typically, an R object or a data set) 
# can be displayed by either help("name") or ?name.

?mean # Access the documentation on a topic with name (e.g. "mean")
?plot # Access the documentation on a topic (e.g. "plot")

??mean # Search the Help System 

#??? Which functions have "help" in them?


## ----accessing_data------------------------------------------------------

cdr3_lengths

### Algebric manipulation
# Algebric operations are vectorized
cdr3_lengths + 1 # Result is each element of the vector +1:

# More operations - ,  /,  %*%


# Elements in column "CDR3", equivalent to dataframe_example[,"CDR3"]

dataframe_example$CDR3 
colnames(dataframe_example) # Returns column names
names(dataframe_example) # Returns column names

### Extract unique values
unique(cdr3_lengths)

### Check data: class, length, summary, structure and head
class(cdr3_lengths) # data class (numeric)

summary(cdr3_lengths) # summary statistics

mean(cdr3_lengths)
min(cdr3_lengths)
max(cdr3_lengths)

# Try finding the standard deviation of 'cdr3_lengths' - Hint: ?sd

t.test(cdr3_lengths, cdr3_lengths)

str(dataframe_example) # Try ?str # Internal structure of an R object

### apply family
apply(matrix_example, 2, nchar) # Count the characters by row (1)

sapply(1:10, function(x) x+1) # see apply, lapply, tapply

## Write your custom function

our_square_function <- function(x){
  x <- x * x
  return(x)
}

our_square_function(3)

#??# Write your own function that adds 1 to every element in a numeric vector
#??# amd perform it on this vector c(3, 6, 2, 6, 2)    

# Control flow statements

### If statement
# if(logical boolean){function to perform if the boolean is TRUE}
x <- 100000000

if(x > 10000){
  "x is a really big number"
}

if(x < 10000){
  "x isn't that big"
} ## Nothing happens here because x < 10000 evaluates to false

### For loop

for(i in 1:10){
  print(i)
} 

### While loop

x <- 0
while(x <= 5){
  print(x)
  x <- x + 1
}
  
#??#  Add an amino acid to every element in cdr3_sequences using a for loop
# Example output:
# [1] "CARIMTERSLKWW"             "CARFEPHVWW"                "CARRWHQYELCSFANWW"
# [4] "CARHLMKDPWVCIFRGTESYAWW"   "CARRQKYIEACTGHNWW"        
# [6] "CARGKYNFIMCQEWTDLVSPRHWW"  "CARRIGPEMCSQDAWW"          "CARVRWIHPKQFDLGCAWW"
# [9] "CARWNVPGHDAKCFSWW"         "CARLYEPISDQRHAWGMVTFNCKWW"


## Note: Many operations in R are vectorized -> use loops only when needed

paste(cdr3_sequences, 'W', sep ="")

# Install R packages ------------------------------------------------------

#install.packages(c("ggplot2", "stringdist" ),
#                 dependencies = TRUE,
#                 repos = "http://cran.us.r-project.org")


## Load R packages

library("ggplot2")
library("stringdist")

stringdist('CARDGYSSGYAMDY', 'CARDGYSLGYAMDY')

# You can also use function without loading the whole package

stringdist::stringdist('CARDGYSSGYAMDY', 'CARDGYSLGYAMDY')

#Set up bioconductor (Additional repository for biological sciences)

source("http://www.bioconductor.org/biocLite.R")
biocLite("seqLogo")  # Package to generate PWM/SeqLogos

## ----plot_data

plot(cdr3_lengths)

# Save as pdf:
pdf(file = 'cdr3_length_plot.pdf', width = 3, height = 3)
plot(cdr3_lengths)
dev.off()

### Line plot of CDR3 lengths

# Store CDR3 length in seperate object
nchar_smpl1 <- nchar(sample1)
nchar_smpl2 <- nchar(sample2)

# Graph CDR3 in function of their length in amino acids
plot(nchar_smpl1, type="o", col="blue", xlab="Antibody clone", 
     ylab="CDR3 length (a.a.)", ylim = c(0, 25))
# Graph lengths of the CDR3 sequences simulated with red dashed line and square points
lines(nchar_smpl2, type="o", pch=22, lty=2, col="red")

# Create a title with a red, bold/italic font
title(main="CDR3 length of antibodies in two samples", col.main="black", font.main=1)


## ----visualization_barplot_boxplot
table(cdr3_lengths)

barplot(table(cdr3_lengths))

boxplot(nchar_smpl1, nchar_smpl2)


## ----visualization_vgene

set.seed(7)

vgene <- dataframe_example$Abundance
names(vgene) <- dataframe_example$Vgene

# Note: prop.table generates frequency table

barplot(prop.table(vgene)*100, ylab = "Frequency (%)")



## ----visualization_venn

### Produce figure of overlapping CDR3 lengths in the two repertoires

# install.packages("VennDiagram") # Uncomment to install package "VennDiagram"

library(VennDiagram) # Load package "VennDiagram"


#Exercise save intersect of sample1 and sample 2 as object

dev.off() # Close the plotting device

cross_area <- intersect(sample1, sample2)

venn_r1_r2 <- draw.pairwise.venn(area1 = length(sample1),
                                 area2 = length(sample2),
                                 cross.area = length(cross_area))

grid.draw(venn_r1_r2)

# Use ?venn.diagram for a detailed explanation of the arguments used


## ----visualization_phylogenies

install.packages('ape')

library(ape)
library(stringdist)

cdr3_cluster <- nj(stringdistmatrix(as.character(cdr3_sequences),
                                    as.character(cdr3_sequences),
                                    method = "lv"))

cdr3_cluster$tip.label = as.character(cdr3_sequences)

plot(cdr3_cluster, cex = 0.5)




# Hands on repertorie -----------------------------------------------------

## Module 3 - loading into R or python, and filtering,
mixcr_output <- read.table("~/Data/mixcr_output_seqs.txt",
           sep=",", header = TRUE, stringsAsFactors=FALSE, fill=TRUE)

## Mixcr output is read in as a data frame with the columns corresponding to 
# annotations regarding CDR sequences, germline gene usage, and clonal frequency

View(mixcr_output)

colnames(mixcr_output)

## Where are the CDR3 sequences stored?
mixcr_output$AA..Seq.CDR3

## How many clones are in the file?
## Each row is considered a clone from mixcr 
nrow(mixcr_output)
## or
length(mixcr_output$AA..Seq.CDR3)

# Singletons are clones with a clonal count (cloneCount) equal to 1
# Sometimes removed to account for PCR/sequencing errors

## How many singletons are there in the file?

singleton_index <- which(mixcr_output$Clone.count==1)
not_singletons <- which(mixcr_output$Clone.count>1)
head(singleton_index)
singleton_number <- length(which(mixcr_output$Clone.count==1))

## What is the percentage of singletons in this data set?
100*singleton_number/nrow(mixcr_output)

## Let's look at the breakdown of all observed clone counts
table(mixcr_output$Clone.count)

prop.table(table(mixcr_output$Clone.count))

## Filter out the singletons and save it to a new object
## use comma to select all columns in a data frame
mixcr_no_singletons <- mixcr_output[not_singletons,]
dim(mixcr_no_singletons)
# for example, to only extract the first two columns of non-singletons
first_twocolumns <- mixcr_output[not_singletons,1:2]
dim(first_twocolumns)


## Now extract just the CDR3 and the clone cpunt from those sequences that ARE singletons
cdrs_only <- mixcr_output[singleton_index,c("Clone.count","AA..Seq.CDR3")]
dim(cdrs_only)
cdrs_only$AA..Seq.CDR3

## _ and * in the CDR3s are out of frame or stop codons respectively
## We will just work on the data without singletons
grepl(pattern="*","CARDYVYWY")
grepl(pattern="\\*","CARDYVYWY")
grepl(pattern="\\*","CAR*DDWY")

non_stopcodon_index <- which(grepl(pattern = "\\*",x = mixcr_no_singletons$AA..Seq.CDR3)==FALSE)
nosingletons_nostopcodon <- mixcr_no_singletons[non_stopcodon_index,]


## Lets remove CDR3s less than 4. 

# what is the minimum CDR3 length in the data set?
min(nchar(nosingletons_nostopcodon$AA..Seq.CDR3))

## how can we find out the index of this element?
which.min(nchar(nosingletons_nostopcodon$AA..Seq.CDR3))

# Module 4 - Summary statistics on filtered data frame

## Clone count distribution
table(nosingletons_nostopcodon$Clone.count)

## Update clone fraction after dropping singletons and stop codons
nosingletons_nostopcodon$new_cloneFraction <- nosingletons_nostopcodon$Clone.count/sum(nosingletons_nostopcodon$Clone.count)
summary(nosingletons_nostopcodon$Clone.count)

## V gene distribution
table(nosingletons_nostopcodon$Best.V.gene)

## VJ gene distribution
## Let's first create a new column in our data object containing the VJ combinations
nosingletons_nostopcodon$VJgenes <- paste(nosingletons_nostopcodon$Best.V.gene,
                                          nosingletons_nostopcodon$Best.J.gene,
                                          sep="")

table(nosingletons_nostopcodon$VJgenes)


## C gene
table(nosingletons_nostopcodon$Best.C.gene)


## CDR3 length distribution
nchar(nosingletons_nostopcodon$AA..Seq.CDR3)
table(nchar(nosingletons_nostopcodon$AA..Seq.CDR3))
summary(nchar(nosingletons_nostopcodon$AA..Seq.CDR3))

## Mutations
###summary(nosingletons_nostopcodon$) ## add still

## How many Serines are there in average per CDR3?
library(stringr)
serines_in_CDRs <- str_count(string = nosingletons_nostopcodon$AA..Seq.CDR3,pattern = "S")
mean(serines_in_CDRs)
summary(serines_in_CDRs)
which.max(serines_in_CDRs)


## now lets normalize by CDR3 length
normalized_serines <- serines_in_CDRs/nchar(nosingletons_nostopcodon$AA..Seq.CDR3)
summary(normalized_serines)

## Amino acid distribution for a given length
source("https://bioconductor.org/biocLite.R")
biocLite("motifStack")
library(motifStack)

desired_length <- 13
cdr3_length13 <- nosingletons_nostopcodon$AA..Seq.CDR3[nchar(nosingletons_nostopcodon$AA..Seq.CDR3)==desired_length]
consensusMatrix(cdr3_length13)


# Module 5 - Visualization

## Clone count distribution
plot(table(nosingletons_nostopcodon$Clone.count))
hist(nosingletons_nostopcodon$Clone.count)
## Update clone fraction after dropping singletons and stop codons
nosingletons_nostopcodon$new_cloneFraction <- nosingletons_nostopcodon$Clone.count/sum(nosingletons_nostopcodon$Clone.count)
summary(nosingletons_nostopcodon$Clone.count)

library(reshape)
library(ggplot2)
clone_fraction_dataframe <- data.frame(Count=nosingletons_nostopcodon$Clone.count,Name=1:nrow(nosingletons_nostopcodon))
clone_fraction_plot <- ggplot(data=clone_fraction_dataframe,aes(x=Name, y=Count)) + geom_point() + theme_bw() + labs(x="Clone ID")
clone_fraction_plot

## V gene distribution
table(nosingletons_nostopcodon$Best.V.gene)
vgene_dataframe <- data.frame(Vgenes=names(table(nosingletons_nostopcodon$Best.V.gene)),
                              Value=as.numeric(table(nosingletons_nostopcodon$Best.V.gene)))
vgene_ggplot <- ggplot(data=vgene_dataframe,aes(x=Vgenes,y=Value,fill=Vgenes)) 
vgene_ggplot <- vgene_ggplot + geom_bar(stat = "identity", position="dodge")
vgene_ggplot <- vgene_ggplot + theme(axis.text.x = element_text(angle = 90))
vgene_ggplot

pdf('~/Figures/vgene_ggplot.pdf')

vgene_ggplot

dev.off()

## VJ gene distribution
## Let's first create a new column in our data object containing the VJ combinations

nosingletons_nostopcodon$VJgenes <- paste(nosingletons_nostopcodon$Best.V.gene,
                                          nosingletons_nostopcodon$Best.J.gene,
                                          sep="-")

vjgene_dataframe <- data.frame(VJgenes=names(table(nosingletons_nostopcodon$VJgenes)),
                              Value=as.numeric(table(nosingletons_nostopcodon$VJgenes)))

vjgene_dataframe <- vjgene_dataframe[order(vjgene_dataframe$Value, decreasing = T)[1:nrow(vjgene_dataframe)], ]

vjgene_ggplot <- ggplot(data=vjgene_dataframe,aes(x=VJgenes,y=Value,fill=VJgenes)) 
vjgene_ggplot <- vjgene_ggplot + geom_bar(stat = "identity", position="dodge")
vjgene_ggplot <- vjgene_ggplot + theme(axis.text.x = element_text(angle = 90))
vjgene_ggplot

## C gene
cgene_dataframe <- data.frame(CGenes = names(table(nosingletons_nostopcodon$Best.C.gene)),
                              Value = as.numeric(table(nosingletons_nostopcodon$Best.C.gene)))

cgene_ggplot <- ggplot(data=cgene_dataframe,aes(x=CGenes,y=Value,fill=CGenes)) 
cgene_ggplot <- cgene_ggplot + geom_bar(stat = "identity", position="dodge")
cgene_ggplot <- cgene_ggplot + theme_bw() + scale_y_continuous(trans='log10')
cgene_ggplot <- cgene_ggplot + theme(axis.text.x = element_text(angle = 90, hjust = 1))
cgene_ggplot

## CDR3 length distribution

nchar(nosingletons_nostopcodon$AA..Seq.CDR3)
table(nchar(nosingletons_nostopcodon$AA..Seq.CDR3))
summary(nchar(nosingletons_nostopcodon$AA..Seq.CDR3))

cdr3_length_dataframe <- data.frame(Length = nchar(nosingletons_nostopcodon$AA..Seq.CDR3))

cdr3_length_plot <- ggplot(data=cdr3_length_dataframe, aes(y = Length, x = 1))
cdr3_length_plot <- cdr3_length_plot + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
cdr3_length_plot <- cdr3_length_plot + theme_bw()
cdr3_length_plot

## How many Serines are there in average per CDR3?
library(stringr)

normalized_aa_list <- list()
aa <- c('S', 'Y', 'W', 'D', 'C')

for(i in 1:5){
  
  normalized_aa_list[[i]] <- str_count(string = nosingletons_nostopcodon$AA..Seq.CDR3,
                                       pattern = aa[i])
}


aa_dataframe <- data.frame('S' = normalized_aa_list[[1]],
                           'Y' = normalized_aa_list[[2]],
                           'W' = normalized_aa_list[[3]],
                           'D' = normalized_aa_list[[4]],
                           'C' = normalized_aa_list[[5]])

aa_dataframe_melt <- melt(aa_dataframe)

aa_plot <- ggplot(data=aa_dataframe_melt, aes(y = value, x =""))
aa_plot <- aa_plot + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
aa_plot <- aa_plot + theme_bw() + facet_grid(.~variable)
aa_plot <- aa_plot + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
aa_plot

## Amino acid distribution for a given length
source("https://bioconductor.org/biocLite.R")
biocLite("motifStack")
library(motifStack)

desired_length <- 13
cdr3_length13 <- nosingletons_nostopcodon$AA..Seq.CDR3[nchar(nosingletons_nostopcodon$AA..Seq.CDR3)==desired_length]

consensusMatrix(cdr3_length13)
logo_variable <- pcm2pfm(consensusMatrix(cdr3_length12))
plot(new("pfm",name="Temp",mat=logo_variable,color=colorset(alphabet = "AA",colorScheme = "chemistry")),ic.scale=F)
plot(new("pfm",name="Temp",mat=logo_variable,color=colorset(alphabet = "AA",colorScheme = "chemistry")),ic.scale=T)

cdr3_length13

# Module 6 - clonotyping

## Basic idea of clonotyping: Group clones together that derive from the same
## precursor cell -> Find lineages of clones

# One basic idea is to cluster (i.e. group) clones by sequence similarity
# Here we use an hierarchical clustering scheme based on edit distance

# Sometimes people also add additional constraints, e.g. fix CDR3 length or
# V- and J-Gene

# Let us first select all the clones with a CDR3 that is 15 AA long

# Find index & extrtact clones

idx_15AA <- nchar(mixcr_output$AA..Seq.CDR3) == 15
clones_15AA <- mixcr_output[idx_15AA, ]

# Extract all clones with V-Gene = "IGHV5-17" an J-Gene = "IGHJ4"

idx_vj <- clones_15AA$Best.V.gene == "IGHV5-17" & clones_15AA$Best.J.gene == "IGHJ4"
clones_to_cluster <- clones_15AA[idx_vj, ]

# Extract CDR3

cdr3s_to_cluster <- clones_to_cluster$AA..Seq.CDR3

# Calculate the edit distance between all extracted sequences

cdr3s_distmat <- stringdistmatrix(cdr3s_to_cluster, cdr3s_to_cluster, 'lv') 

# Normalize by length -> so that we can cluster based on %-similarity

cdr3s_distmat <- cdr3s_distmat/15

# Turn into distance matrix (i.e. make sure that matrix is symmetric)

cdr3s_distmat <- as.dist(cdr3s_distmat)

# Calculate clustering based on cdr3s_distmat

cdr3s_clustering <- hclust(cdr3s_distmat, method = 'single')

# Let's have a look at the clustering object

summary(cdr3s_clustering)

# Let's combine sequences with at least 80% similarity

cdr3s_final_clusters <- cutree(cdr3s_clustering, h = 0.2)

# Inspect final results:

cdr3s_final_clusters

table(cdr3s_final_clusters)

#??# We find one large cluster -> what is its label and how many clones does it contain?

# Let's select all sequences from this cluster

idxs_cluster <- cdr3s_final_clusters == 3

cdr3s_to_cluster[idxs_cluster] 


#??# Obviously it would be nice to cluster a complete dataset based on these
#??# parameters (Same V/J-Gene, same length, 80% similarity)
#??# This is implemented in the code below: Try to understand what it does

# Define variables that are needed

Combined_Datasets <- list()
Clusters <- list()

# Calculate the length for all sequences

mixcr_output$Len <- nchar(mixcr_output$AA..Seq.CDR3)

# Helper function to select the relevant subsets

select_columns <- function(x){
  x_new <- dplyr::select(x, c('AA..Seq.CDR3', 'Best.V.Gene', 'Best.J.Gene'))
  return(x_new)
}

# Helper function to calculate distance matrix

dist_calc <- function(x){
  
  dist_mat <- as.dist(stringdistmatrix(x, x, method = 'hamming')/nchar(x[1]))
  return(dist_mat)
  
}

# Function to compute clonotypes

clonotyping <- function(data){
  
  meta_list <- split(data, 
                     list(data$Best.V.gene, 
                          data$Best.J.gene, 
                          data$Len))
  
  idxs <- which(sapply(meta_list, function(x) length(x$Clone.ID)) > 0)
  
  meta_list <- meta_list[idxs]
  
  dist_mat <- lapply(meta_list, function(x) dist_calc(x$AA..Seq.CDR3))
  
  # Hierarchical clustering step, complete linkage, cut tree at 80% similarity
  
  clusts <- lapply(dist_mat, function(x) {
    if(length(x) > 0){
      return(cutree(hclust(x, method = 'single'), h = 0.2))
    } else {
      return(1)
    }
  }
  )
  
  # Needed to increase the clonotype numbering correctly
  add_nr <- 0
  
  # Renumber clonotypes 
  for(i in 1:length(clusts)){
    clusts[[i]] <- clusts[[i]] + add_nr
    add_nr <- max(clusts[[i]])
  }
  
  meta_list <- do.call(rbind, meta_list)
  meta_list$Clonotype <- unlist(clusts)
  
  
  return(unique(meta_list))
}
  
# Calculate complete clustering

mixcr_output_clonotyped <- clonotyping(mixcr_output)

# Module 7 - finding clones from known sequence,
# e.g. how many diff v gene, exact same cdr3, v gene/cdr3 length.

# Let's assume a previous experiment identified the clone with the
# CDR3 'CARDYYDYGYAMDYW' as a binding antibody -> How can we find more variants
# of that clone?

# First let's see which V-Gene and J-Gene this clone uses:

idx_clone <- mixcr_output_clonotyped$AA..Seq.CDR3 == 'CARRRPYAMDYW'

mixcr_output_clonotyped[idx_clone,]

#??# What can you say about this clone?

#Find further variants from the clonotyping:

relevant_cl <- mixcr_output_clonotyped$Clonotype[idx_clone]
idx_cl <- mixcr_output_clonotyped$Clonotype == relevant_cl

potential_variants <- mixcr_output_clonotyped[idx_cl, ]

potential_variants$AA..Seq.CDR3

 
# Module 8 - Phylogenetics and Networks -----------------------------------

###  Both networks and phylogenetics can start from a distance matrix,
# although this is not always the case. 


library(stringdist)
lv_distance_cdr3s <- stringdistmatrix(mixcr_output$AA..Seq.CDR3[1:500],mixcr_output$AA..Seq.CDR3[1:500],method="lv")
## lv_distance_cdr3s is a matrix containing the distances from all-against-all comparison

library(pheatmap) ### Heatmaps can be very useful to visualize large matrices 
## the more red, the higher number of mutations needed to change one string to another
pheatmap(lv_distance_cdr3s,cluster_rows = F,cluster_cols  = F)

## Networks are made of vertices with edges connecting them
## often for CDR3-based networks, an edge is drawn between two CDR3 nodes if it is more similar than a set threshold

# first we will save the output from stringdistmatrix as a new object (our future adjacency matrix)
adjacency_matrix_1 <- lv_distance_cdr3s

# then we will draw edges between points who are separated by one mutation (e.g. levenshtein distance is <1)
adjacency_matrix_1[adjacency_matrix_1<2] <- 1
# similarly, we will exclude edges between points who are separated by multiple mutations (e.g. levenshtein distance is >1)
adjacency_matrix_1[adjacency_matrix_1>=2] <- 0
# then we can visualize now how the distances look
pheatmap(adjacency_matrix_1,cluster_rows = F,cluster_cols = F)
# then we will set the diagonal values to NA, as we do not want self-connections
diag(adjacency_matrix_1) <- NA

## To compare, we will also make a network using a larger threshold 
adjacency_matrix_3 <- lv_distance_cdr3s
adjacency_matrix_3[adjacency_matrix_3<4] <- 1
adjacency_matrix_3[adjacency_matrix_3>=4] <- 0
pheatmap(adjacency_matrix_3,cluster_rows = F,cluster_cols = F)
diag(adjacency_matrix_3) <- NA

library(igraph)
graph_lv_1 <- graph_from_adjacency_matrix(adjacency_matrix_1, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
graph_lv_3 <- graph_from_adjacency_matrix(adjacency_matrix_3, mode = c("undirected"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)

par(mfrow=c(1,2))
plot(graph_lv_1,vertex.size=1,vertex.label=NA,main="Levenshtein distance 1")
plot(graph_lv_3,vertex.size=1,vertex.label=NA, main="Levenshtein distance 3")

## Graph statistics 
## Number of edges - i.e. how many lines are there connecting each point
number_of_edges_lv1 <- length(E(graph_lv_1))
number_of_edges_lv3 <- length(E(graph_lv_3))

## Density of a graph is the ratio of the number of edges and the number of possible edges
density_lv1 <- graph.density(graph_lv_1)
density_lv3 <- graph.density(graph_lv_3)

## Calculates the connected components of a graph  - the number of subgraphs is given by $no
cluster_number_lv1 <- clusters(graph_lv_1)$no
cluster_number_lv3 <- clusters(graph_lv_3)$no

### Transitivity measures the probability that the adjacent vertices are connected. 
# also known as clustering coefficient
transitive_lv1 <- transitivity(graph_lv_1,type="global")
transitive_lv3 <- transitivity(graph_lv_3,type="global")


######## Phylogenetics

## Usually phylogenetic trees are inferred using full length (entire VDJ region) sequences that arise from a single VDJ recombination event

## Typical workflow involves 
# 1. Assigning sequences to germline root 
# 2. Phylogenetic inference method
# 3. Tree visualization and analysis 

# 1. While some specialized software exists for this (Clonify, Partis, SONAR), many people simply
# group sequences by germline usage, sometimes adding restrictions on CDR3 length or sequence homology

# in our example we can work with the most used V gene, which is IGHV5-17

## first want to get all CDR3 clones using this V gene
ighv5_17 <- mixcr_output[mixcr_output$Best.V.gene=="IGHV5-17",]
ighv5_17$cdr3_length <- nchar(ighv5_17$AA..Seq.CDR3)

## now can split by cdr3 length and j gene

distance_matrix_tree1 <- stringdistmatrix(ighv5_17$AA..Seq.CDR3, ighv5_17$AA..Seq.CDR3,method = "lv")
library(ape)
my_tree <- nj(distance_matrix_tree1)
class(my_tree)
attributes(my_tree)
my_tree$tip.label <- ighv5_17$AA..Seq.CDR3
plot(my_tree,cex=.3)
dev.off()


save.image(file="~/PHD/pegs/pegs_temp.RData")
load("~/Downloads/pegs_temp.RData")
