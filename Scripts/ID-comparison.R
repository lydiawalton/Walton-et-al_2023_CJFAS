dev.off() #Shut down open graphics devices
rm(list = ls(all=T)) #Clear environment

#Setting your working directory
library(here) 
#Data manipulation and analysis
library(dplyr)
library(geomorph)
library(quadcleanR)
#Data visualization
library(kableExtra)

#Import csv files
genetics <- read.csv("Data/Myct_geneticIDs.csv") #Genetic identifications
otolith <- read.csv("Data/Myct_otolithIDs.csv") #Otolith identifications


#----------Identification methods comparison

#Merge the two csv's together so you can compare the ID's
mergedID <- merge(otolith, genetics, by="Sample.ID")

#Remove the barcode non-compliant (compliant column = no) and flagged = yes
mergedID.filtered <- mergedID %>% 
  filter(Barcode.Compliant == "Yes" & Flagged.Record == "")

#Now compare the Identification column (genetics) to the Otolith.ID column (otoliths)
mergedID.subset <- mergedID.filtered %>% 
  subset(select=c(Sample.ID, Otolith.ID, Collection.Date, Identification))

#Rename the genetics column
mergedID.subset <- mergedID.subset %>% 
  rename("Genetic.ID" = "Identification")

#Compare if the answers in the Otolith.ID column match the answers in the Genetic.ID column
mergedID.comparison <- mergedID.subset %>% 
  mutate(Final_ID = if_else(Otolith.ID == Genetic.ID, "Match", "Mismatch"))

#Create a table showing the ID comparison
(mergedID.table <- mergedID.comparison %>% 
  group_by(Final_ID) %>% 
  summarise("Otolith ID" = length(Final_ID)))

#We have 179 Otolith IDs that match the genetics and 37 that do not


#Make a table showing which species were mismatched
(Mismatch.table <- mergedID.comparison %>% 
    group_by(Otolith.ID, Genetic.ID, Final_ID) %>% 
    filter(Final_ID == "Mismatch") %>%   
    summarise("Number of IDs" = length(Final_ID)))










