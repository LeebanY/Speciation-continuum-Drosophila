library(tidyverse)
library(gt)
library(plotly)
library(scales)
library(ggpubr)
library(ape)
library(data.table)
library(viridis)
library(stringr)

Drosophila_sampling <- read.csv('UpdatedDrosophilomicsProjectMasterSpreadsheet.csv')
#Look at dataset
glimpse(Drosophila_sampling)

######################### PREPROCESS DATA TO GET IMPORTANT INFO ########################


#Get accurate boolean identification of whether samples are pooled male and female, just female, or just male.
#Code asking if any element in any column contains sex information, then add TRUE/FALSE classification to new column.
Drosophila_sampling$MaleFemale <- !!rowSums(sapply(Drosophila_sampling[1:96], grepl, pattern = "pooledmandf"))
Drosophila_sampling$Female <- !!rowSums(sapply(Drosophila_sampling[1:96], grepl, pattern = "\\bfemale\\b"))
Drosophila_sampling$Male <- !!rowSums(sapply(Drosophila_sampling[1:96], grepl, pattern = "\\bmale\\b"))


#Convert true and false to integer representation
Drosophila_sampling$MaleFemale<-as.integer(as.logical(Drosophila_sampling$MaleFemale))
Drosophila_sampling$Female<-as.integer(as.logical(Drosophila_sampling$Female))
Drosophila_sampling$Male<-as.integer(as.logical(Drosophila_sampling$Male))
glimpse(Drosophila_sampling)

#Now for the case where this is no sex label, we want to assign that to the Unknown column as 1 (TRUE).
Drosophila_sampling<-Drosophila_sampling %>% 
  mutate(Unknown= case_when(Male == 0 & Female == 0 & MaleFemale == 0 ~ 1, #if no males/females/pooled, then we categorise as unknown
                            Male == 1 & Female == 0 & MaleFemale == 0 ~ 0, #if male, then unknown is false
                            Male == 0 & Female == 1 & MaleFemale == 0 ~ 0, #if female, then unknown is false
                            Male == 0 & Female == 0 & MaleFemale == 1 ~ 0)) #if pooled male and female, unknown is false
                            


###############     TIDY, FILTER AND OUTPUT FINAL TABLE
#First do some preliminary filtering of data to get rid of anything with particularly
#problematic features, for example, anything in 'crosses' is not suitable, or pooled individuals.




#First do some preliminary filtering of data to get rid of anything with particularly
#problematic features, for example, anything in 'crosses' is not suitable, or pooled individuals.


Randomly_sampled_dataset<-Drosophila_sampling %>%
  select(Organism,Run, Assay.Type, Bases, BioSample, BioProject, #Select variables that we want to keep in our dataset
         Bytes, Center.Name, Consent, DATASTORE.filetype, Experiment,
         Instrument, LibraryLayout, LibrarySource, Platform, ReleaseDate,
         Sample.Name, Male, Female, MaleFemale, Unknown) %>%
  filter_all(any_vars(!str_detect(., pattern = "cross"))) %>% #remove any samples containing crosses -- we do not any samples to have been manipulated.
  filter(Assay.Type =='WGS', Bases > 1000000000, Bases < 50,000,000,000) %>% #Preferably we want samples with large enough, but not too large a file size. We preferably don't want pooled male and female samples.
  group_by(Organism) %>% #Group samples by organism.
  slice_sample() #randomly sample one sample per organism. 


#Quick check to see if we've not removed any species with filtering. 
length(unique(Randomly_sampled_dataset$Organism))
length(unique(Drosophila_sampling$Organism))


#Now output the code in tsv and csv formats. 
write.csv(Randomly_sampled_dataset, 'FinalRandomlySampledDrosophilomicsDataset.csv')



 ##############   Redundant code #############
#Check what species are missing
# #Read in table with all pairs
# Drosophila_pairs <- read.csv('UpdatedDrosophilomicsProjectMasterSpreadsheetLeebanPairsNonIndependent.csv')
# #adjust names so that they're equivalent
# Drosophila_pairs$Species_A <- paste("Drosophila", Drosophila_pairs$Species_A)
# Drosophila_pairs$Species_B <- paste("Drosophila", Drosophila_pairs$Species_B)
# Drosophila_pairs$Species_A<-gsub("yacuba", "yakuba", Drosophila_pairs$Species_A)
# Drosophila_pairs$Species_A<-gsub("sulfurigaster neonasuta", "neonasuta", Drosophila_pairs$Species_A)
# Drosophila_pairs$Species_A<-gsub("kepulauna", "kepulauana", Drosophila_pairs$Species_A)
# Drosophila_pairs$Species_A<-gsub("loweii", "lowei", Drosophila_pairs$Species_A)
# Drosophila_pairs$Species_A<-gsub("melanogaster ", "melanogaster", Drosophila_pairs$Species_A)
# #Then use intersect function to extract those that do match and outersect for those that do not.
# intersect(Drosophila_pairs$Species_A, Randomly_sampled_dataset$Organism)
# setdiff(Drosophila_pairs$Species_A, Randomly_sampled_dataset$Organism)


#Checks over. 