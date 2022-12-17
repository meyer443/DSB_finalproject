# Header ------------------------------------------------------------------

# This is a script to take Michigan DNR trawl data and calculate yearly species biomass estimates
# Weighed fish are used to estimate the mass of unweighed fish of the same species
# Relative abundances are estimated using fish weights and counts
# Smoothing plots are used to visualize trends over time at various scales
# An NMDS plot is made to visualize trends in fish community composition through time

# Author: Justin Meyer
# Date: 2022-12-09
# Email: meyer443@purdue.edu
# github data repository: https://github.com/meyer443/DSB_finalproject

# Set directory -----------------------------------------------------------
rm(list=ls())
setwd("H:/DataScienceFinalProject")

# Load required packages --------------------------------------------------
library(tidyverse)
library(readxl)
# library(ggplot2)
library(GGally)
library(vegan)

# create function to use later
`%notin%` <- Negate(`%in%`)

# load data ---------------------------------------------------------------
card1a <- read_xlsx("./datasets/SagBay2009-2021.xlsx", sheet = 1)
card1b <- read_xlsx("./datasets/SagBay1970-2008.xlsx", sheet = 1)
card2a <- read_xlsx("./datasets/SagBay2009-2021.xlsx", sheet = 2)
card2b <- read_xlsx("./datasets/SagBay1970-2008.xlsx", sheet = 2)
card3a <- read_xlsx("./datasets/SagBay2009-2021.xlsx", sheet = 3)
card3b <- read_xlsx("./datasets/SagBay1970-2008.xlsx", sheet = 3) # warnings not about data we care about
card4a <- read_xlsx("./datasets/SagBay2009-2021.xlsx", sheet = 5)
card4b <- read_xlsx("./datasets/SagBay1970-2008.xlsx", sheet = 5)
card5a <- read_xlsx("./datasets/SagBay2009-2021.xlsx", sheet = 6)
card5b <- read_xlsx("./datasets/SagBay1970-2008.xlsx", sheet = 6)
sppcodes <- read_xlsx("./datasets/SpeciesCodes.xlsx")

# append cards
card1 <- card1a %>% bind_rows(card1b)
card2 <- card2a %>% bind_rows(card2b)
card3 <- card3a %>% bind_rows(card3b)
card4 <- card4a %>% bind_rows(card4b)
card5 <- card5a %>% bind_rows(card5b)

# Card 5 ------------------------------------------------------------------

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                 Summarize Biomass and Unweighed Fish                # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# card 5
# isolate age classes, differentiating between weighed and unweighed fish

# convert all unknown ages (0+) to code 5
card5 <- 
  card5 %>%
  mutate(Age = replace(Age, Age == 0, 5)) %>% 
  mutate(Age = replace(Age, Age == 7, 5)) %>% 
  mutate(Age = replace(Age, Age %in% NA, 5))

# Age 0 fish --------------------------------------------------------------

# for 'AGE == 1' age 0 fish with weights
bulkdf0 <- 
  card5 %>% 
  group_by(SYear, Species) %>% 
  filter(Age == 1) %>% 
  filter(Totwt != 'NA' &
           Totwt != 0) %>% 
  summarise(Age = mean(Age),
            biomass = sum(Totwt),
            numAge = sum(Totnum)) %>% 
  filter(numAge %notin% NA)

# age 0 fish without weights
bulkdf0.noweight <- 
  card5 %>% 
  group_by(SYear, Species) %>% 
  filter(Age == 1) %>% 
  filter(Totwt == 'NA' |
           Totwt == 0 | 
           Totwt %in% NA) %>% 
  summarise(Age = mean(Age),
            numAge = sum(Totnum)) %>% 
  filter(numAge %notin% NA)

# Adult fish --------------------------------------------------------------

# for 'AGE = 2-4', for age 1, 1+, and 2+ fish with weights
bulkdfAdult <- 
  card5 %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Age == 2 |
           Age == 3 |
           Age == 4) %>% 
  filter(Totwt != 'NA' &
           Totwt != 0) %>% 
  summarise(biomass = sum(Totwt),
            numAge = sum(Totnum))

# for 'AGE = 2-4' weighout weights
bulkdfAdult.noweight <- 
  card5 %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Age == 2 |
           Age == 3 |
           Age == 4) %>% 
  filter(Totwt == 'NA' |
           Totwt == 0) %>% 
  summarise(numAge = sum(Totnum))

# Unaged fish -------------------------------------------------------------

# for AGE == 5, age 0+ fish with weights
bulkdf5 <- 
  card5 %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Age == 5) %>%
  filter(Totwt != 'NA' &
           Totwt != 0) %>% 
  summarise(biomass = sum(Totwt),
            numAge = sum(Totnum))

# for AGE ==  5 without weights
bulkdf5.noweight <- 
  card5 %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Age == 5) %>% 
  filter(Totwt == 'NA' |
           Totwt == 0) %>% 
  summarise(numAge = sum(Totnum)) %>% 
  filter(numAge %notin% NA)

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                           Biomass Estimates                         # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# Estimate mass of Age 0 fish ---------------------------------------------

# average weight of age 0 fish
avg.weight5 <-
  bulkdf0 %>% 
  group_by(Species, Age) %>% 
  filter(numAge != 'NA' & 
           numAge != 0) %>%
  summarise(avgweight = sum(biomass)/sum(numAge))

# join no weight list with average weight based on species id, mutate biomass label by avgweight*totnum
noweight5 <- 
  bulkdf0.noweight %>% 
  left_join(avg.weight5) %>%
  mutate(biomass = numAge*avgweight) %>% 
  select(-avgweight)

# age 0 all
card5.age0 <- 
  bulkdf0 %>%
  bind_rows(noweight5)

# Estimate mass of adult fish ---------------------------------------------

# avg weight of adult ages for each species
avg.weight.adult5 <-
  bulkdfAdult %>% 
  group_by(Species, Age) %>% 
  summarise(avgweight = sum(biomass)/sum(numAge))

# join to no weight list
# remove remaining NAs, not enough data to fill in (fish not caught any other time)
bulkadult.noweight <- 
  bulkdfAdult.noweight %>% 
  left_join(avg.weight.adult5) %>% 
  filter(avgweight != 'NA') %>%
  mutate(biomass = numAge*avgweight) %>%
  select(-avgweight)

# all adults
card5.adult <- 
  bulkdfAdult %>% 
  bind_rows(bulkadult.noweight)

# all weights
allweights5 <- 
  avg.weight5 %>% 
  bind_rows(avg.weight.adult5)

# Estimate mass of unaged fish --------------------------------------------

# unknown ages
age.proportion5 <- 
  card5.age0 %>%
  bind_rows(card5.adult) %>%
  group_by(Species, Age) %>%
  summarise(count = sum(numAge)) %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  mutate(prop = count/sum(count))

# join with all weights
propweight5 <- 
  age.proportion5 %>% 
  left_join(allweights5) %>% 
  group_by(Species) %>%
  summarise(weight = sum(prop*avgweight))

# apply age proportions to totnum counts
addweights5 <- 
  bulkdf5.noweight %>% 
  left_join(propweight5) %>% 
  group_by(SYear, Species, Age) %>%
  summarise(numAge = numAge,
            biomass = numAge*weight)

# all unknown ages
card5.unknown <- 
  bulkdf5 %>% 
  bind_rows(addweights5)

# compile all ages back into one tibble for this card
card5.weights <- 
  card5.age0 %>% 
  bind_rows(card5.adult) %>% 
  bind_rows(card5.unknown) %>% 
  filter(biomass %notin% NA)

# Card 4 ------------------------------------------------------------------

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                 Summarize Biomass and Unweighed Fish                # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# card 4
# isolate age classes, differentiating between weighed and unweighed fish

# convert all unknown age codes to 5
card4 <- 
  card4 %>%
  mutate(AGE = replace(AGE, AGE == 0, 5))

# there is 1 row where the weight and count = 0, so I am removing it since I can't be sure the fish was actually caught
card4 <- 
  card4 %>%
  filter(NUMBER != 0)

# Age 0 fish --------------------------------------------------------------

# for 'AGE == 1' age 0 fish with weights
df0 <- 
  card4 %>% 
  group_by(SYear, SPECIES) %>% 
  filter(AGE == 1) %>% 
  filter(Weight != 'NA' &
           Weight != 0) %>% 
  summarise(AGE = mean(AGE),
            biomass = sum(Weight),
            numAge = sum(NUMBER))

# age 0 fish without weights
df0.noweight <- 
  card4 %>% 
  group_by(SYear, SPECIES) %>% 
  filter(AGE == 1) %>% 
  filter(Weight == 'NA' | 
           Weight == 0) %>% 
  summarise(AGE = mean(AGE),
            numAge = sum(NUMBER))

# Adult fish --------------------------------------------------------------

# for 'AGE = 2-4', for age 1, 1+, and 2+ fish with weights
dfAdult <- 
  card4 %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  filter(AGE == 2 |
           AGE == 3 |
           AGE == 4) %>% 
  filter(Weight != 'NA' &
           Weight != 0) %>% 
  summarise(biomass = sum(Weight),
            numAge = sum(NUMBER))

# for 'AGE = 2-4' weighout weights
dfAdult.noweight <- 
  card4 %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  filter(AGE == 2 |
           AGE == 3 |
           AGE == 4) %>% 
  filter(Weight == 'NA' |
           Weight == 0) %>% 
  summarise(numAge = sum(NUMBER))

# Unaged fish -------------------------------------------------------------

# for 'AGE == 5', age 0+ fish with weights
df5 <- 
  card4 %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  filter(AGE == 5) %>%
  filter(Weight != 'NA' &
           Weight != 0) %>% 
  summarise(biomass = sum(Weight),
            numAge = sum(NUMBER))

# for 'AGE ==  5' without weights
df5.noweight <-
  card4 %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  filter(AGE == 5) %>% 
  filter(Weight == 'NA' |
           Weight == 0) %>% 
  summarise(totnum = sum(NUMBER))

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                           Biomass Estimates                         # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# Estimate mass of Age 0 fish -----------------------------------------

# avg weight of age 0 for each species
avg.weight <- 
  df0 %>% 
  group_by(SPECIES, AGE) %>% 
  summarise(avgweight = sum(biomass)/sum(numAge))

# join no weight list with average weight based on species id, mutate biomass label by avgweight*totnum
noweight <-
  df0.noweight %>% 
  left_join(avg.weight) %>% 
  mutate(biomass = numAge*avgweight) %>% 
  select(-avgweight)

# age 0 all
card4.age0 <- 
  df0 %>% 
  bind_rows(noweight)

# Estimate mass of adult fish ---------------------------------------------

# avg weight of adult ages for each species
avg.weight.adult <-
  dfAdult %>% 
  group_by(SPECIES, AGE) %>% 
  summarise(avgweight = sum(biomass)/sum(numAge))

# join to no weight list
# remove remaining NAs, not enough species weight data to fill in
adult.noweight <- 
  dfAdult.noweight %>% 
  left_join(avg.weight.adult) %>% 
  filter(avgweight != 'NA') %>% 
  mutate(biomass = numAge*avgweight) %>%
  select(-avgweight)

# all adults
card4.adult <- 
  dfAdult %>% 
  bind_rows(adult.noweight)

# all weights
allweights <- 
  avg.weight %>% 
  bind_rows(avg.weight.adult)

# Estimate mass of unaged fish (age 0+) -----------------------------------

# calculate age structure from sheet
age.proportion <- 
  card4.age0 %>% 
  bind_rows(card4.adult) %>%
  group_by(SPECIES, AGE) %>% 
  summarise(count = sum(numAge)) %>% 
  ungroup() %>% 
  group_by(SPECIES) %>% 
  mutate(prop = count/sum(count))

# join proportions with allweights
propweight <- 
  age.proportion %>% 
  left_join(allweights) %>% 
  group_by(SPECIES) %>% 
  summarise(avgweight = sum(prop*avgweight))

# apply age proportions to totnum counts
addweights <- 
  df5.noweight %>% 
  left_join(propweight) %>% 
  filter(avgweight != 'NA') %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  summarise(numAge = totnum, 
            biomass = totnum*avgweight)

# all unknown ages
card4.unknown <- 
  df5 %>% 
  bind_rows(addweights)

# compile all ages
card4.weights <- 
  card4.age0 %>% 
  bind_rows(card4.adult) %>%
  bind_rows(card4.unknown)

# Card 2 ------------------------------------------------------------------

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                 Summarize Biomass and Unweighed Fish                # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# card 2
# isolate age classes, differentiating between weighed and unweighed fish
# all aged fish have weights

# edit ages so 0, 7, and NA are 5 (age 0+/unknown)
card2 <- 
  card2 %>% 
  mutate(Age = replace(Age, Age == 0, 5)) %>%
  mutate(Age = replace(Age, Age == 7, 5)) %>%
  mutate(Age = replace(Age, Age %in% NA, 5))


# age 0 with weights
foragedf0 <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 1) %>% 
  group_by(SYear, Species) %>% 
  filter(Weight != 0) %>% 
  summarise(Age = mean(Age),
            biomass = sum(Weight),
            numAge = sum(Totnum))

# age 0 fish without weights
foragedf0.noweight <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 1) %>% 
  group_by(SYear, Species) %>% 
  filter(Weight == 0 |
           Weight == 'NA') %>% 
  summarise(Age = mean(Age),
            numAge = sum(Totnum))

# Adult fish --------------------------------------------------------------

# age 1, 2, 3 with weights
foragedfAdult <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 2 |
           Age == 3 | 
           Age == 4) %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Weight != 0) %>% 
  summarise(biomass = sum(Weight),
            numAge = sum(Totnum))

# adult fish of unknown weight
foragedfAdult.noweight <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 2 |
           Age == 3 | 
           Age == 4) %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Weight == 0 |
           Weight == 'NA') %>% 
  summarise(biomass = sum(Weight),
            numAge = sum(Totnum))

# Unaged fish -------------------------------------------------------------

# unknown ages with weights
foragedf5 <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 5) %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Weight != 0) %>% 
  summarise(Age = mean(Age),
            biomass = sum(Weight),
            numAge = sum(Totnum))

# unknown ages no weights
foragedf5.noweight <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 5) %>% 
  filter(Totnum != 0) %>% 
  filter(Weight == 0 |
           Weight == 'NA') %>% 
  group_by(SYear, Species, Age) %>% 
  summarise(numAge = sum(Totnum))

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                           Biomass Estimates                         # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# Estimate mass of Age 0 fish ---------------------------------------------

# average weight of age 0 fish
avg.weight2 <- 
  foragedf0 %>% 
  group_by(Species, Age) %>% 
  filter(numAge != 'NA' & 
           numAge != 0) %>%
  summarise(avgweight = sum(biomass)/sum(numAge))

# Adult fish --------------------------------------------------------------

# avg weight of adult ages for each species
avg.weight.adult2 <- 
  foragedfAdult %>%
  group_by(Species, Age) %>%
  summarise(avgweight = sum(biomass)/sum(numAge))

# join average weights
allweights2 <- 
  avg.weight2 %>% 
  bind_rows(avg.weight.adult2)

# Unaged fish -------------------------------------------------------------

# calculate age structure
age.proportion2 <- 
  foragedf0 %>%
  bind_rows(foragedfAdult) %>% 
  group_by(Species, Age) %>%
  summarise(count = sum(numAge)) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(prop = count/sum(count))

# join with allweights
propweight2 <- 
  age.proportion2 %>% 
  left_join(allweights2) %>% 
  group_by(Species) %>% 
  summarise(weight = sum(prop*avgweight))

# apply age proportions to totnum counts
addweights2 <- 
  foragedf5.noweight %>% 
  left_join(propweight2) %>% 
  group_by(SYear, Species, Age) %>%
  summarise(numAge = numAge, 
            biomass = numAge*weight)

# compile card together into single tibble
card2.weights <- 
  foragedf0 %>% 
  bind_rows(foragedfAdult) %>%
  bind_rows(foragedf5) %>%
  bind_rows(addweights2) %>% 
  filter(biomass %notin% NA)

# Card 3 ------------------------------------------------------------------

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                 Summarize Biomass and Unweighed Fish                # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# Calculate yearly species biomass and counts -----------------------------

# card 3
# isolate age classes, differentiating between weighed and unweighed fish

# age = actual age
# convert ages to match other cards
card3 <- 
  card3 %>% 
  mutate(Age = replace(Age, Age > 4 & Age < 99, 4)) %>%
  mutate(Age = replace(Age, Age == -9, 5)) %>% 
  mutate(Age = replace(Age, Age == 99, 5)) %>% 
  mutate(Age = replace(Age, Age == 3, 4)) %>%
  mutate(Age = replace(Age, Age == 2, 3)) %>% 
  mutate(Age = replace(Age, Age == 1, 2)) %>%
  mutate(Age = replace(Age, Age == 0, 1)) %>% 
  mutate(Age = replace(Age, Age %in% NA, 5))

# All weighed fish --------------------------------------------------------

# all fish with weights
card3biomass <- 
  card3 %>% 
  filter(Weight != 'NA' & 
           Weight != 0) %>% 
  group_by(SYear, Species, Age) %>% 
  summarise(biomass = sum(Weight),
            numAge = n())

# pull out age classes
tagdf0 <-
  card3biomass %>%
  filter(Age == 1)

tagdfAdult <- 
  card3biomass %>%
  filter(Age != 1 & Age != 5)

tagdfunknown <- 
  card3biomass %>%
  filter(Age == 5)

# Unweighed fish ----------------------------------------------------------

# all fish without weights
card3.noweight <- 
  card3 %>% 
  filter(Weight == "NA" |
           Weight == 0) %>%
  group_by(SYear, Species, Age) %>%
  summarise(biomass = sum(Weight), 
            numAge = n())

# pull out age classes
tagdf0.noweight <- 
  card3.noweight %>% 
  filter(Age == 1)

tagdfAdult.noweight <- 
  card3.noweight %>%
  filter(Age != 1 &
           Age != 5)

tagdfunknown.noweight <- 
  card3.noweight %>% 
  filter(Age == 5)

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                           Biomass Estimates                         # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# Estimate mass of Age 0 fish ---------------------------------------------

# average weight of age 0 fish
avg.weight3 <- 
  tagdf0 %>% 
  group_by(Species, Age) %>% 
  summarise(avgweight = sum(biomass)/sum(numAge))

# age 0 no weights
noweight3 <- 
  tagdf0.noweight %>%
  left_join(avg.weight3) %>%
  mutate(biomass = numAge*avgweight) %>%
  select(-avgweight)

# age 0 all
card3.age0 <- 
  tagdf0 %>%
  bind_rows(noweight3)

# Adult fish --------------------------------------------------------------

# average weight of aged fish
avg.weight.adult3 <- 
  tagdfAdult %>% 
  group_by(Species, Age) %>% 
  filter(numAge != 'NA' & numAge != 0) %>%
  summarise(avgweight = sum(biomass)/sum(numAge))

# aged fish no weights
adult.noweight3 <- 
  tagdfAdult.noweight %>% 
  left_join(avg.weight.adult3) %>%
  mutate(biomass = numAge*avgweight) %>%
  select(-avgweight)

# all adults
card3.adult <- 
  tagdfAdult %>% 
  bind_rows(adult.noweight3)

# all weights
tag.weights <- 
  avg.weight3 %>%
  bind_rows(avg.weight.adult3)

# Unaged fish -------------------------------------------------------------

# age structure
age.proportion3 <- 
  card3.age0 %>%
  bind_rows(card3.adult) %>%
  group_by(Species, Age) %>%
  summarise(count = sum(numAge)) %>% 
  ungroup() %>%
  group_by(Species) %>% 
  mutate(prop = count/sum(count))

# join with allweights
propweight3 <- 
  age.proportion3 %>%
  left_join(tag.weights) %>%
  group_by(Species) %>%
  summarise(avgweight = sum(prop*avgweight))

# fill in weights
noweight3 <- 
  card3.noweight %>%
  left_join(propweight3) %>%
  mutate(biomass = numAge*avgweight) %>% 
  select(-avgweight)

# all unaged fish
card3.unknown <- 
  tagdfunknown %>%
  bind_rows(noweight3)

# all weights all fish
card3.weights <- 
  card3.age0 %>% 
  bind_rows(card3.adult) %>% 
  bind_rows(card3.unknown)

# Plots -------------------------------------------------------------------

# change column names on card 4 so they match the other cards
card4.weights <- 
  card4.weights %>%
  rename(Age = AGE, 
         Species = SPECIES)

# catch info
allfish <- 
  card2.weights %>% 
  bind_rows(card3.weights) %>% 
  bind_rows(card4.weights) %>% 
  bind_rows(card5.weights)

# write_csv(allfish, "biomass.csv")

# find any remaining 0s
allfish %>% 
  filter(Age == 0 | Age == 'NA' | Age %in% NA |
           biomass == 0 | biomass == 'NA' | biomass %in% NA | 
           numAge == 0 | numAge == 'NA' | numAge %in% NA) %>%
  print(n = 29)
# all remaining 0/NAs are for fish counts that have weights (weighed on board the boat but not enumerated) so I am keeping them for biomass estimates
allfish <- 
  allfish %>%
  mutate(numAge = replace(numAge, numAge %in% NA, 0))

# Standardize effort ------------------------------------------------------

# remove unidentified fish, dreissenid mussels, and the invasive sea lamprey since they are not informative of commmunity productivity
# join species info from 'sppcodes' object
df.main <- 
  allfish %>% 
  left_join(sppcodes) %>%
  filter(Species != 999 &
           Species != 998 &
           Species != 920 &
           Species != 5) %>% 
  mutate(Func_group = str_replace(Func_group, "_", " "))

# convert species id code to character from double
df.main$Species <- as.character(df.main$Species)

# === need to standardize per trawl and log-transform biomass === #

# convert NA towtime to 0, these are tows that got hung up
card1 <- 
  card1 %>%
  mutate(Towtime = replace(Towtime, Towtime %in% NA, 0))

# create trawls/year tibble from card1 then join it to dataframe
trawl <- 
  card1 %>%
  group_by(SYear) %>%
  summarise(towtime = sum(Towtime))

# join
all.trawl <- 
  df.main %>%
  left_join(trawl)

# write_csv(all.trawl, "catch_trawl.csv")

# create new biomass/trawl column
fishpertrawl <-
  all.trawl %>%
  mutate(BioPerTrawl = biomass/towtime)
# fishpertrawl will be used for nmds later

# write_csv(fishpertrawl, "fishpertrawl.csv")

# Bay-wide ----------------------------------------------------------------

# get yearly bay-wide values
yearly <- 
  fishpertrawl %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

# log-transform biomass since fish growth is not linear
forgraph <- 
  yearly %>% 
  mutate(Log_Biomass = log10(biomass),
         Log_Biomass_per_Trawl = log10(bioper))

# low sampling effort in 1970 and only yellow perch reported
forgraph <- forgraph %>% filter(SYear != 1970)

# bay-wide biomass/trawl graph
ggplot(forgraph) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  xlab("Year") +
  ylab(expression(paste("Log biomass (kg " %*% " minute"^"-1"*")"))) +
  ggtitle("Fish biomass per tow minute") +
  theme_bw(base_size = 16)
# ggsave(filename = "biomassplot.jpg", device='jpg', dpi=700)

# Functional group graph ------------------------------------------

# biomass per trawl grouped by feeding functional group
feedingpertrawl <-   
  fishpertrawl %>%
  group_by(SYear, Func_group) %>%
  summarise(biomass = sum(biomass),
            bioper = sum(BioPerTrawl))
# write_csv(feedingpertrawl, "FuncEffort.csv")

# log-transform bio
# get rid of parasite data, not interested
feedgraph <- 
  feedingpertrawl %>% 
  mutate(Log_Biomass = log10(biomass),
         Log_Biomass_per_Trawl = log10(bioper)) %>% 
  filter(Func_group != 'Parasite')

# graph of feeding functional group biomass per trawl
# I added a color-blind friendly color palette after reading the report rubric
# cbbPalette and scale_colour_manual are the only changes from the figure in my class presentation

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(feedgraph) +
  aes(x = SYear, y = Log_Biomass_per_Trawl, group = Func_group) +
  # geom_point(aes(shape = Func_group, colour = Func_group), size = 2.5) +
  geom_smooth(aes(color = Func_group), se = TRUE, linewidth = 1.5) +
  xlab("Year") +
  ylab(expression(paste("Log biomass (kg " %*% " minute"^"-1"*")"))) +
  ggtitle("Feeding functional group biomass per tow minute") +
  labs(color="Functional Group") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values=cbbPalette)
# ggsave(filename = "functionalplot.jpg", device='jpg', dpi=700)

# Pairs plot --------------------------------------------------------------

# create object to compare groups to one another
func.compare <- 
  feedgraph %>% 
  group_by(SYear) %>% 
  select(-bioper, -Log_Biomass, -biomass) %>%
  pivot_wider(names_from="Func_group", values_from="Log_Biomass_per_Trawl") %>%
  filter(SYear != 1970) %>% 
  ungroup() %>% 
  rename('Ben invert' = 'Benthic invertivore', 'Gen invert' = 'General invertivore') %>% 
  select(-SYear)

# plot pairwise correlations
ggpairs(func.compare) +
  ggtitle("Functional Group Biomass Correlations") +
  theme_bw(base_size = 14)
# ggsave(filename = "pairsplot.jpg", device='jpg', dpi=700)

# Calculate catch-per-unit-effort -----------------------------------------

# cpue all fish
cpue <- 
  fishpertrawl %>%
  ungroup() %>%
  group_by(SYear) %>%
  summarise(cpue = sum(numAge)/towtime) %>%
  unique()

# log-transform
cpue <- 
  cpue %>% 
  mutate(logcpue = log10(cpue))

# plot
ggplot(cpue) +
  aes(x = SYear, y = logcpue) +
  geom_point(size = 2.5) +
  geom_smooth(linewidth = 1.5, se = TRUE) +
  xlab("Year") +
  ylab(expression(paste("Log cpue (no. fish " %*% " minute"^"-1"*")"))) +
  ggtitle("Count per tow minute") +
  theme_bw(base_size = 16)
# ggsave(filename = "cpueplot.jpg", device='jpg', dpi=700)

# NMDS --------------------------------------------------------------------

# create decade column for nmds grouping
fishpertrawl <- 
  fishpertrawl %>% 
  mutate(decade = ifelse(SYear < 1980, '70s',
                         ifelse(SYear > 1979 & SYear < 1990, '80s',
                                ifelse(SYear > 1989 & SYear < 2000, '90s',
                                       ifelse(SYear > 1999 & SYear < 2010, '00s',
                                              ifelse(SYear > 2009 & SYear < 2020, '10s',
                                                     ifelse (SYear > 2019, '20s', SYear)))))),
         Func_group = str_replace(Func_group, "_", " "))

# convert data types
fishpertrawl$Species <- as.character(fishpertrawl$Species)
fishpertrawl$SYear <- as.character(fishpertrawl$SYear)

# Species nmds ------------------------------------------------------------

# in order to help nmds converge, only keep species present in at least 15 years
# sum cpue/tow for each species-year combination and pivot wider for matrix conversion
fish <- 
  fishpertrawl %>% 
  filter(Species == 106 |
           Species == 108 |
           Species == 109 |
           Species == 119 |
           Species == 131 |
           Species == 132 | 
           Species == 133 |
           Species == 134 | 
           Species == 203 | 
           Species == 402 | 
           Species == 403 |
           Species == 405 |
           Species == 413 | 
           Species == 504 | 
           Species == 508 | 
           Species == 511 | 
           Species == 601 | 
           Species == 604 |
           Species == 706 | 
           Species == 707 |
           Species == 801 |
           Species == 803 |
           Species == 906) %>%
  filter(SYear != 1970) %>% 
  group_by(SYear, Common_name, decade) %>% 
  summarise(BioPerTrawl = sum(BioPerTrawl)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Common_name", values_from = BioPerTrawl)

# NAs are for years where certain species were not caught, change this to 0 for relative abundance estimates in nmds
# object created will be used for nmds plot creation
fish <- 
  fish %>%
  replace(is.na(.), 0)

# select only species data
com <- 
  fish %>%
  select(-SYear, -decade) %>% 
  as.matrix()

# set variable for nmds
set.seed(123)

# run Bray-Curtis nmds on matrix
nmds <- metaMDS(com, distance = 'bray')

# generate nmds plot with ellipses and labels
x <- -1.5:1.5
y <- -2:1
# tiff("ellipse_species.tiff", width = 9.85, height = 6.32, units = 'in', res = 300)
plot(x,y, type = "n", main='Species NMDS', xlab='Axis 1', ylab='Axis 2')
ordiellipse(nmds, fish$decade, draw = "polygon", label = TRUE)
orditorp(nmds,display="species",col="red",air=1,cex=1)
# dev.off()