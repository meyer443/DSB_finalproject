# Header ------------------------------------------------------------------

# Author name: Justin Meyer
# Date: 2022-10-21
# Notes: script to process Michigan DNR formatted Saginaw Bay trawl data
# and calculate species and feeding functional group abundances, then create a
# plot of abundances through time
# catch per unit effort = cpue = relative abundance

# to be updated later with 1970-2008 data

# file setup --------------------------------------------------------------

rm(list = ls())
setwd("G:/SagBayComm")

# load required packages --------------------------------------------------

library(tidyverse)
library(readxl)

# load data ---------------------------------------------------------------

# MDNR refers to their tables as "cards"
card1 <- read_xlsx("sagbay.xlsx", sheet = 1)
card2 <- read_xlsx("sagbay.xlsx", sheet = 2)
card3 <- read_xlsx("sagbay.xlsx", sheet = 3)
card4lengths <- read_xlsx("sagbay.xlsx", sheet = 4)
card4weights <- read_xlsx("sagbay.xlsx", sheet = 5)
card5 <- read_xlsx("sagbay.xlsx", sheet = 6)
# feeding functional group information
functional <- read_xlsx("Func_final.xlsx")

# prepare data for joins --------------------------------------------------

# trawl info, select identifying columns
card1cpue <- 
  card1 %>% 
  select(SYear, Idnum, SDate, Towtime)

# non-forage species sheet 1/2
card2cpue <- 
  card2 %>% 
  select(SYear, Idnum, Species, Totnum) %>% 
  group_by(SYear, Idnum, Species) %>%
  distinct() %>% # rows are duplicated for being split by length group
  ungroup()

# non-forage sheet 2/2
card3cpue <- 
  card3 %>% 
  select(SYear, Idnum, Species) %>% 
  group_by(SYear, Idnum, Species) %>% 
  count(Species, name = "Totnum") %>% # count fish and match columns from card2
  ungroup()

# forage sheet 1/3
card4lengthscpue <- 
  card4lengths %>% 
  select(SYear, Idnum, Species, N) %>% # N is this sheet's Totnum
  group_by(SYear, Idnum, Species) %>% 
  summarise(Totnum = sum(N)) %>% # combine species together
  ungroup()

# forage sheet 2/3
card4weightscpue <- 
  card4weights %>% 
  select(SYear, Idnum, Species, Number) %>% # Number = Totnum
  group_by(SYear, Idnum, Species) %>% 
  summarise(Totnum = sum(Number)) %>% # tally species
  ungroup()

# forage sheet 3/3
card5cpue <- 
  card5 %>% 
  select(SYear, Idnum, Species, Totnum) %>% 
  group_by(SYear, Idnum, Species) %>%
  summarise(Totnum = sum(Totnum)) %>% # tally species
  ungroup()

# functional groups
func <- 
  functional %>% 
  group_by(Func_ID) %>% 
  filter(Func_ID != 6) %>% # parasites ommitted from this analysis
  
# bind cards 2-5 ----------------------------------------------------------

# columns match from data processing, bind together
catch <- 
  card2cpue %>% 
  bind_rows(card3cpue) %>% 
  bind_rows(card4lengthscpue) %>%
  bind_rows(card4weightscpue) %>% 
  bind_rows(card5cpue) %>% 
  group_by(SYear, Idnum, Species) %>% 
  summarise(Totnum = sum(Totnum))

# execute joins -----------------------------------------------------------

# join using year-specific trawl ID
cpue <- 
  card1cpue %>% 
  left_join(catch) %>% 
  group_by(SYear, Idnum, Species) %>% 
  summarise(cpue = Totnum/Towtime) %>% # calculate cpue for each species per trawl
  ungroup() %>% 
  group_by(SYear, Species) %>% 
  summarise(cpue = mean(cpue)) # mean cpue for species per year

# write_csv(cpue, "speciescpue.csv")

# join to feeding functional group relational table using species codes
cpuedata <- 
  cpue %>% 
  left_join(func) %>%
  na.omit() # unidentified species omitted

# write_csv(cpuedata, "cpue_species_func.csv")

# calculate functional group cpue
meancpue <-
  cpuedata %>%
  group_by(SYear, Func_ID) %>%
  summarise(meancpue = mean(cpue)) # yearly mean cpue for functional groups

# factorise for plot
meancpue$Func_ID <- 
  factor(meancpue$Func_ID, 
            levels = c(1, 2, 3, 4, 5),
            labels = c("Benthic invertivore", "Planktivore", 
                       "General invertivore", "Piscivore", "Omnivore"))

# plot --------------------------------------------------------------------

ggplot(meancpue) +
  aes(x = SYear, y = meancpue, group=Func_ID, color=Func_ID) +
  geom_line(size = 1) +
  xlab("Year") +
  ylab("CPUE") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("Functional Group Abundance")
