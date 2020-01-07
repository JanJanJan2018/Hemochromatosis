
# These files are from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121620

# Hemochromatosis is a disorder that can be fatal, it is when someone's
# body has too much iron in it or their gut takes out more iron from food than normal
# it poison's the blood system.

# This study analyzes 40-65 year old males and females who have two genotypes of 
# homozygous strands and one of a heterozygous genotype for the hemochromatosis
# gene compared to healthy 40-65 year old males and females. 


# The platform for the GPL16686 family is 6.3 GB in file size
# But a lot of the genes are NA, so read in the first 53,621 observations
# the same as the length of the hemotochromatosis gene expressions file


# hemo <- read.delim('GPL16686_family.soft', sep='\t', comment.char='#',
#                    na.strings=c('',' '), skip=8087, header=FALSE, nrows=53621)

iron <- read.csv('hemotochromatosisPatientsGeneExpressionValues.csv', sep=',',
                 header=TRUE, na.strings=c('',' '), skip=3)
ironMeta <- read.csv('hemotochromatosisPatientsGeneExpressionValues.csv', sep=',',
                     header=FALSE, na.strings=c('',' '), nrows=4)

# Hemo <- hemo[2:53621, c(1,6)]
# colnames(Hemo) <- c('ID','GB_ACC')
# Hemo <- Hemo[3:53619,]
# Hemo1 <- Hemo[complete.cases(Hemo$GB_ACC),]
# row.names(Hemo1) <- NULL


# combined <- merge(Hemo1, iron, by.x='ID', by.y='ID_REF')

###
# The gene symbols can be taken from ensebl AFFY HuGene 2 0 st v1 probe in bioMart
# this file is too large for Github and is 157.4 mb in size it is split apart below

# BioMart <- read.delim('mart_export.txt', sep='\t', header=TRUE,
#                       na.strings=c('',' '))
# 
# # the following is the above file split into smaller file sizes to load in
# 
# mart1 <-bioMart[1:211142,]
# mart2 <-bioMart[211143:422284,]
# mart3 <-bioMart[422285:633426,]
# mart4 <-bioMart[633427:844569,]
# mart5 <-bioMart[844570:1055712,]
# mart6 <-bioMart[1055713:1266855,]
# mart7 <-bioMart[1266856:1477998,]
# mart8 <-bioMart[1477999:1689134,]
# 
# write.csv(mart1, 'mart1.csv', row.names=FALSE)
# write.csv(mart2, 'mart2.csv', row.names=FALSE)
# write.csv(mart3, 'mart3.csv', row.names=FALSE)
# write.csv(mart4, 'mart4.csv', row.names=FALSE)
# write.csv(mart5, 'mart5.csv', row.names=FALSE)
# write.csv(mart6, 'mart6.csv', row.names=FALSE)
# write.csv(mart7, 'mart7.csv', row.names=FALSE)
# write.csv(mart8, 'mart8.csv', row.names=FALSE)
###

mart1 <- read.csv('mart1.csv', sep=',', header=TRUE, na.strings=c('',' '))
mart2 <- read.csv('mart2.csv', sep=',', header=TRUE, na.strings=c('',' '))
mart3 <- read.csv('mart3.csv', sep=',', header=TRUE, na.strings=c('',' '))
mart4 <- read.csv('mart4.csv', sep=',', header=TRUE, na.strings=c('',' '))
mart5 <- read.csv('mart5.csv', sep=',', header=TRUE, na.strings=c('',' '))
mart6 <- read.csv('mart6.csv', sep=',', header=TRUE, na.strings=c('',' '))
mart7 <- read.csv('mart7.csv', sep=',', header=TRUE, na.strings=c('',' '))
mart8 <- read.csv('mart8.csv', sep=',', header=TRUE, na.strings=c('',' '))

BioMart <- rbind(mart1, mart2, mart3, mart4, mart5, mart6, mart7, mart8)

# The biomart file has the probe IDs used in the platform but not by accession number

#merge the bioMart file with the iron file using the HGNC Symbol

bio <- BioMart[,c(2:4)]
Bio <- bio[complete.cases(bio$HGNC.symbol),]
Bio <- Bio[!duplicated(Bio),]

Blend <- merge(iron, Bio, by.x='ID_REF', by.y='AFFY.HuGene.2.0.st.v1.probe')
# There are still duplicate HGNC names because of the probe ID identifying the gene
# where it is a pseudo gene

Blend1 <- Blend[,c(2:25,27)]

Blend-meta <- Blend[,26:27]

# group by genes and get the gene means so that there is only one unique gene in the table

library(dplyr)
samples <- colnames(Blend1[1:24]) 
samples <- as.vector(samples)
geneMeans <- Blend1 %>% group_by(HGNC.symbol) %>% summarise_at(samples, mean, na.rm = TRUE)


healthyMales <- geneMeans[,c(1,20:22,25)]
healthyFemales <- geneMeans[,c(1,23,24)]

#ignore the groups where 1 and 2 are one type of homozygous strand, 
# and 3 is another or a heterozygous strand
hemoMales <- geneMeans[,c(1,2:7,9,10,12:15,19)]
hemoFemales <- geneMeans[,c(1,8,11,16:18)]

write.csv(geneMeans, 'hemochromatosis-cleaned.csv', row.names=FALSE)
write.csv(healthyMales, 'hemoHealthyMales.csv', row.names=FALSE)
write.csv(healthyFemales, 'hemoHealthyFemales.csv', row.names=FALSE)
write.csv(hemoMales, 'hemoMales.csv', row.names=FALSE)
write.csv(hemoFemales, 'hemoFemales.csv', row.names=FALSE)


