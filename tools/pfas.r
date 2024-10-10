library(tidyverse)
rm(list = ls())

source('./Instructions,_Source_File,_and_Test_Data_for_the_LCMRL_Project/100518 MRL.LCMRL.Stats.r')


## Data file conversion
files <- list.files('./data')
str_sub(files[17], end = -5)

for(i in 1:length(files)) {
  file.path <- paste('./data', files[i], sep = '/')
  comp.name <- str_sub(files[i], end = -5)
  dat <- read.table(file = file.path, header = TRUE)
  cn <- as.character(read.table(file = file.path)[1,])
  colnames(dat) <- cn
  dat.long <- gather(dat, cn, key = 'Spike', value = 'Result')
  lrb <- data.frame(Spike = c(0, 0, 0, 0), Result = c(0, 0, 0, 0))
  dat.long <- rbind(lrb, dat.long)
  dat.long <- add_column(dat.long, Analyte = comp.name, .before = 'Spike')
  dat.long <- add_column(dat.long, Lab = 'NIES', .after = 'Analyte')
  dat.long <- add_column(dat.long, 'Dilution Factor' = 1, .after = 'Result')
  dat.long <- add_column(dat.long, Units = 'ng/l', .after = 'Dilution Factor')
  write.table(dat.long, file = 'pfas-lcmrl.csv', append = TRUE, row.names = FALSE, col.names = FALSE, sep = ',')
}

filename <- 'pfas-lcmrl.csv'
LCMRL.Values(filename, rnnr = 1)
LCMRL.Graphs(filename, rnnr = 1)

library(lcmrl4jecs)
fname <- 'pfas-lcmrl.csv'
LCMRL.Values(fname, rnnr = 1)
LCMRL.Graphs(fname, rnnr = 1)

## Extra
pfas1 <- read.csv('pfas-lcmrl.csv', header = FALSE)
data(test)
coln <- colnames(test) 
colnames(pfas1) <- coln
names(pfas1)
saveRDS(pfas1, 'pfas.rda')
readRDS('pfas.rda')



