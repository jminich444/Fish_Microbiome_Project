##mantel test for phylogenetic distance and gut microbiota distance
##subset of species that were found in timetree

library("ape")
library("ade4")
library("spaa")

####gill####

#Load phylogenetic tree
fish_phylo_gill <- read.tree(file = "C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/species_gill_sub_timetree.nwk")
fish_phylo_gill_divtime <- cophenetic.phylo(fish_phylo_gill)

fish_phylo_gill_divtime <- fish_phylo_gill_divtime[order(rownames(fish_phylo_gill_divtime)) , order(colnames(fish_phylo_gill_divtime))]
fish_phylo_gill_divtime_dist <- as.dist(fish_phylo_gill_divtime)

#read metadata
fish_meta_gill <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/fish_meta_gill_sub_phylo.csv", sep = ";", dec = ".")

#unweighted unifrac
fish_dm_uwuf_gill_nr <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/non-rarified/fish_dm_uwuf_gill_nr_phylo.tsv", sep = "\t", dec = ".", row.names = 1, header= TRUE, check.names = FALSE)

rownames(fish_dm_uwuf_gill_nr) <- fish_meta_gill$species_name_timetree[match(rownames(fish_dm_uwuf_gill_nr), fish_meta_gill$id)]
colnames(fish_dm_uwuf_gill_nr) <- fish_meta_gill$species_name_timetree[match(colnames(fish_dm_uwuf_gill_nr), fish_meta_gill$id)]

fish_dm_uwuf_gill_nr <- fish_dm_uwuf_gill_nr[order(rownames(fish_dm_uwuf_gill_nr)) , order(colnames(fish_dm_uwuf_gill_nr))]

fish_dm_uwuf_gill_nr_dist <- as.dist(fish_dm_uwuf_gill_nr)
fish_dm_uwuf_gill_nr_matrix <- as.matrix(fish_dm_uwuf_gill_nr)

#weighted unifrac
fish_dm_wuf_gill_nr <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/non-rarified/fish_dm_wuf_gill_nr_phylo.tsv", sep = "\t", dec = ".", row.names = 1, header= TRUE, check.names = FALSE)

rownames(fish_dm_wuf_gill_nr) <- fish_meta_gill$species_name_timetree[match(rownames(fish_dm_wuf_gill_nr), fish_meta_gill$id)]
colnames(fish_dm_wuf_gill_nr) <- fish_meta_gill$species_name_timetree[match(colnames(fish_dm_wuf_gill_nr), fish_meta_gill$id)]

fish_dm_wuf_gill_nr <- fish_dm_wuf_gill_nr[order(rownames(fish_dm_wuf_gill_nr)) , order(colnames(fish_dm_wuf_gill_nr))]

fish_dm_wuf_gill_nr_dist <- as.dist(fish_dm_wuf_gill_nr)
fish_dm_wuf_gill_nr_matrix <- as.matrix(fish_dm_wuf_gill_nr)

#Mantel test

#unweighted unifrac
mantel.test(fish_phylo_gill_divtime, fish_dm_uwuf_gill_nr_matrix, nperm = 999, graph = TRUE)
mantel.rtest(fish_phylo_gill_divtime_dist, fish_dm_uwuf_gill_nr_dist, nrepet = 999)
mantel(fish_phylo_gill_divtime_dist, fish_dm_uwuf_gill_nr_dist, permutations = 999, method = "spearman")
mantel(fish_phylo_gill_divtime_dist, fish_dm_uwuf_gill_nr_dist, permutations = 999, method = "pearson")

#weighted unifrac
mantel.test(fish_phylo_gill_divtime, fish_dm_wuf_gill_nr_matrix, nperm = 999, graph = TRUE)
mantel.rtest(fish_phylo_gill_divtime_dist, fish_dm_wuf_gill_nr_dist, nrepet = 999)
mantel(fish_phylo_gill_divtime_dist, fish_dm_wuf_gill_nr_dist, permutations = 999, method = "spearman")
mantel(fish_phylo_gill_divtime_dist, fish_dm_wuf_gill_nr_dist, permutations = 999, method = "pearson")

#Make data frame in long format from distance matrix
fish_gill_divtime_df <- dist2list(fish_phylo_gill_divtime_dist)
fish_gill_divtime_df$value <- fish_gill_divtime_df$value/2
fish_uwuf_gill_df <- dist2list(fish_dm_uwuf_gill_nr_dist)
fish_wuf_gill_df <- dist2list(fish_dm_wuf_gill_nr_dist)

#Remove rows with comparisons of the same species
fish_gill_divtime_df2 <- fish_gill_divtime_df[fish_gill_divtime_df$value != 0,]
fish_uwuf_gill_df2 <- fish_uwuf_gill_df[fish_uwuf_gill_df$value != 0,]
fish_wuf_gill_df2 <- fish_wuf_gill_df[fish_wuf_gill_df$value != 0,]

plot(fish_uwuf_gill_df2$value ~ fish_gill_divtime_df2$value, ylab = "gut mb distance", xlab = "divergence time", main = "unweighted unifrac-gills", pch = 16)
abline(lm(fish_uwuf_gill_df2$value ~ fish_gill_divtime_df2$value))
summary(lm(fish_uwuf_gill_df2$value ~ fish_gill_divtime_df2$value))

plot(fish_wuf_gill_df2$value ~ fish_gill_divtime_df2$value, ylab = "gut mb distance", xlab = "divergence time", main = "weighted unifrac-gills", pch = 16)
abline(lm(fish_wuf_gill_df2$value ~ fish_gill_divtime_df2$value))
summary(lm(fish_wuf_gill_df2$value ~ fish_gill_divtime_df2$value))

####hindgut####

#Load phylogenetic tree
fish_phylo_hg <- read.tree(file = "C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/species_hg_sub_timetree.nwk")
fish_phylo_hg_divtime <- cophenetic.phylo(fish_phylo_hg)

fish_phylo_hg_divtime <- fish_phylo_hg_divtime[order(rownames(fish_phylo_hg_divtime)) , order(colnames(fish_phylo_hg_divtime))]

fish_phylo_hg_divtime_dist <- as.dist(fish_phylo_hg_divtime)

#read metadata
fish_meta_hg <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/fish_meta_hg_sub_phylo.csv", sep = ";", dec = ".")

#unweighted unifrac
fish_dm_uwuf_hg_nr <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/non-rarified/fish_dm_uwuf_hg_nr_phylo.tsv", sep = "\t", dec = ".", row.names = 1, header= TRUE, check.names = FALSE)

rownames(fish_dm_uwuf_hg_nr) <- fish_meta_hg$species_name_timetree[match(rownames(fish_dm_uwuf_hg_nr), fish_meta_hg$id)]
colnames(fish_dm_uwuf_hg_nr) <- fish_meta_hg$species_name_timetree[match(colnames(fish_dm_uwuf_hg_nr), fish_meta_hg$id)]

fish_dm_uwuf_hg_nr <- fish_dm_uwuf_hg_nr[order(rownames(fish_dm_uwuf_hg_nr)) , order(colnames(fish_dm_uwuf_hg_nr))]

fish_dm_uwuf_hg_nr_dist <- as.dist(fish_dm_uwuf_hg_nr)
fish_dm_uwuf_hg_nr_matrix <- as.matrix(fish_dm_uwuf_hg_nr)

#weighted unifrac
fish_dm_wuf_hg_nr <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/non-rarified/fish_dm_wuf_hg_nr_phylo.tsv", sep = "\t", dec = ".", row.names = 1, header= TRUE, check.names = FALSE)

rownames(fish_dm_wuf_hg_nr) <- fish_meta_hg$species_name_timetree[match(rownames(fish_dm_wuf_hg_nr), fish_meta_hg$id)]
colnames(fish_dm_wuf_hg_nr) <- fish_meta_hg$species_name_timetree[match(colnames(fish_dm_wuf_hg_nr), fish_meta_hg$id)]

fish_dm_wuf_hg_nr <- fish_dm_wuf_hg_nr[order(rownames(fish_dm_wuf_hg_nr)) , order(colnames(fish_dm_wuf_hg_nr))]

fish_dm_wuf_hg_nr_dist <- as.dist(fish_dm_wuf_hg_nr)
fish_dm_wuf_hg_nr_matrix <- as.matrix(fish_dm_wuf_hg_nr)

#Mantel test

#unweighted unifrac
mantel.test(fish_phylo_hg_divtime, fish_dm_uwuf_hg_nr_matrix, nperm = 999, graph = TRUE)
mantel.rtest(fish_phylo_hg_divtime_dist, fish_dm_uwuf_hg_nr_dist, nrepet = 999)
mantel(fish_phylo_hg_divtime_dist, fish_dm_uwuf_hg_nr_dist, permutations = 999, method = "spearman")
mantel(fish_phylo_hg_divtime_dist, fish_dm_uwuf_hg_nr_dist, permutations = 999, method = "pearson")

#weighted unifrac
mantel.test(fish_phylo_hg_divtime, fish_dm_wuf_hg_nr_matrix, nperm = 999, graph = TRUE)
mantel.rtest(fish_phylo_hg_divtime_dist, fish_dm_wuf_hg_nr_dist, nrepet = 999)
mantel(fish_phylo_hg_divtime_dist, fish_dm_wuf_hg_nr_dist, permutations = 999, method = "spearman")
mantel(fish_phylo_hg_divtime_dist, fish_dm_wuf_hg_nr_dist, permutations = 999, method = "pearson")

#Make data frame in long format from distance matrix
fish_hg_divtime_df <- dist2list(fish_phylo_hg_divtime_dist)
fish_hg_divtime_df$value <- fish_hg_divtime_df$value/2
fish_uwuf_hg_df <- dist2list(fish_dm_uwuf_hg_nr_dist)
fish_wuf_hg_df <- dist2list(fish_dm_wuf_hg_nr_dist)

#Remove rows with comparisons of the same species
fish_hg_divtime_df2 <- fish_hg_divtime_df[fish_hg_divtime_df$value != 0,]
fish_uwuf_hg_df2 <- fish_uwuf_hg_df[fish_uwuf_hg_df$value != 0,]
fish_wuf_hg_df2 <- fish_wuf_hg_df[fish_wuf_hg_df$value != 0,]

plot(fish_uwuf_hg_df2$value ~ fish_hg_divtime_df2$value, ylab = "gut mb distance", xlab = "divergence time", main = "unweighted unifrac-hg", pch = 16)
abline(lm(fish_uwuf_hg_df2$value ~ fish_hg_divtime_df2$value))
summary(lm(fish_uwuf_hg_df2$value ~ fish_hg_divtime_df2$value))

plot(fish_wuf_hg_df2$value ~ fish_hg_divtime_df2$value, ylab = "gut mb distance", xlab = "divergence time", main = "weighted unifrac-hg", pch = 16)
abline(lm(fish_wuf_hg_df2$value ~ fish_hg_divtime_df2$value))
summary(lm(fish_wuf_hg_df2$value ~ fish_hg_divtime_df2$value))

####midgut####

#Load phylogenetic tree
fish_phylo_mg <- read.tree(file = "C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/species_mg_sub_timetree.nwk")
fish_phylo_mg_divtime <- cophenetic.phylo(fish_phylo_mg)

fish_phylo_mg_divtime <- fish_phylo_mg_divtime[order(rownames(fish_phylo_mg_divtime)) , order(colnames(fish_phylo_mg_divtime))]

fish_phylo_mg_divtime_dist <- as.dist(fish_phylo_mg_divtime)

#read metadata
fish_meta_mg <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/fish_meta_mg_sub_phylo.csv", sep = ";", dec = ".")

#unweighted unifrac
fish_dm_uwuf_mg_nr <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/non-rarified/fish_dm_uwuf_mg_nr_phylo.tsv", sep = "\t", dec = ".", row.names = 1, header= TRUE, check.names = FALSE)

rownames(fish_dm_uwuf_mg_nr) <- fish_meta_mg$species_name_timetree[match(rownames(fish_dm_uwuf_mg_nr), fish_meta_mg$id)]
colnames(fish_dm_uwuf_mg_nr) <- fish_meta_mg$species_name_timetree[match(colnames(fish_dm_uwuf_mg_nr), fish_meta_mg$id)]

fish_dm_uwuf_mg_nr <- fish_dm_uwuf_mg_nr[order(rownames(fish_dm_uwuf_mg_nr)) , order(colnames(fish_dm_uwuf_mg_nr))]

fish_dm_uwuf_mg_nr_dist <- as.dist(fish_dm_uwuf_mg_nr)
fish_dm_uwuf_mg_nr_matrix <- as.matrix(fish_dm_uwuf_mg_nr)

#weighted unifrac
fish_dm_wuf_mg_nr <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/non-rarified/fish_dm_wuf_mg_nr_phylo.tsv", sep = "\t", dec = ".", row.names = 1, header= TRUE, check.names = FALSE)

rownames(fish_dm_wuf_mg_nr) <- fish_meta_mg$species_name_timetree[match(rownames(fish_dm_wuf_mg_nr), fish_meta_mg$id)]
colnames(fish_dm_wuf_mg_nr) <- fish_meta_mg$species_name_timetree[match(colnames(fish_dm_wuf_mg_nr), fish_meta_mg$id)]

fish_dm_wuf_mg_nr <- fish_dm_wuf_mg_nr[order(rownames(fish_dm_wuf_mg_nr)) , order(colnames(fish_dm_wuf_mg_nr))]

fish_dm_wuf_mg_nr_dist <- as.dist(fish_dm_wuf_mg_nr)
fish_dm_wuf_mg_nr_matrix <- as.matrix(fish_dm_wuf_mg_nr)

#Mantel test

#unweighted unifrac
mantel.test(fish_phylo_mg_divtime, fish_dm_uwuf_mg_nr_matrix, nperm = 999, graph = TRUE)
mantel.rtest(fish_phylo_mg_divtime_dist, fish_dm_uwuf_mg_nr_dist, nrepet = 999)
mantel(fish_phylo_mg_divtime_dist, fish_dm_uwuf_mg_nr_dist, permutations = 999, method = "spearman")
mantel(fish_phylo_mg_divtime_dist, fish_dm_uwuf_mg_nr_dist, permutations = 999, method = "pearson")

#weighted unifrac
mantel.test(fish_phylo_mg_divtime, fish_dm_wuf_mg_nr_matrix, nperm = 999, graph = TRUE)
mantel.rtest(fish_phylo_mg_divtime_dist, fish_dm_wuf_mg_nr_dist, nrepet = 999)
mantel(fish_phylo_mg_divtime_dist, fish_dm_wuf_mg_nr_dist, permutations = 999, method = "spearman")
mantel(fish_phylo_mg_divtime_dist, fish_dm_wuf_mg_nr_dist, permutations = 999, method = "pearson")

#Make data frame in long format from distance matrix
fish_mg_divtime_df <- dist2list(fish_phylo_mg_divtime_dist)
fish_mg_divtime_df$value <- fish_mg_divtime_df$value/2
fish_uwuf_mg_df <- dist2list(fish_dm_uwuf_mg_nr_dist)
fish_wuf_mg_df <- dist2list(fish_dm_wuf_mg_nr_dist)

#Remove rows with comparisons of the same species
fish_mg_divtime_df2 <- fish_mg_divtime_df[fish_mg_divtime_df$value != 0,]
fish_uwuf_mg_df2 <- fish_uwuf_mg_df[fish_uwuf_mg_df$value != 0,]
fish_wuf_mg_df2 <- fish_wuf_mg_df[fish_wuf_mg_df$value != 0,]

plot(fish_uwuf_mg_df2$value ~ fish_mg_divtime_df2$value, ylab = "gut mb distance", xlab = "divergence time", main = "unweighted unifrac-mg", pch = 16)
abline(lm(fish_uwuf_mg_df2$value ~ fish_mg_divtime_df2$value))
summary(lm(fish_uwuf_mg_df2$value ~ fish_mg_divtime_df2$value))

plot(fish_wuf_mg_df2$value ~ fish_mg_divtime_df2$value, ylab = "gut mb distance", xlab = "divergence time", main = "weighted unifrac-mg", pch = 16)
abline(lm(fish_wuf_mg_df2$value ~ fish_mg_divtime_df2$value))
summary(lm(fish_wuf_mg_df2$value ~ fish_mg_divtime_df2$value))

####skin####

#Load phylogenetic tree
fish_phylo_skin <- read.tree(file = "C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/species_skin_sub_timetree.nwk")
fish_phylo_skin_divtime <- cophenetic.phylo(fish_phylo_skin)

fish_phylo_skin_divtime <- fish_phylo_skin_divtime[order(rownames(fish_phylo_skin_divtime)) , order(colnames(fish_phylo_skin_divtime))]

fish_phylo_skin_divtime_dist <- as.dist(fish_phylo_skin_divtime)

#read metadata
fish_meta_skin <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/fish_meta_skin_sub_phylo.csv", sep = ";", dec = ".")

#unweighted unifrac
fish_dm_uwuf_skin_nr <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/non-rarified/fish_dm_uwuf_skin_nr_phylo.tsv", sep = "\t", dec = ".", row.names = 1, header= TRUE, check.names = FALSE)

rownames(fish_dm_uwuf_skin_nr) <- fish_meta_skin$species_name_timetree[match(rownames(fish_dm_uwuf_skin_nr), fish_meta_skin$id)]
colnames(fish_dm_uwuf_skin_nr) <- fish_meta_skin$species_name_timetree[match(colnames(fish_dm_uwuf_skin_nr), fish_meta_skin$id)]

fish_dm_uwuf_skin_nr <- fish_dm_uwuf_skin_nr[order(rownames(fish_dm_uwuf_skin_nr)) , order(colnames(fish_dm_uwuf_skin_nr))]

fish_dm_uwuf_skin_nr_dist <- as.dist(fish_dm_uwuf_skin_nr)
fish_dm_uwuf_skin_nr_matrix <- as.matrix(fish_dm_uwuf_skin_nr)

#weighted unifrac
fish_dm_wuf_skin_nr <- read.csv("C:/Users/Andre/Dropbox/Work/UCSD/Lab work/Metaanalysis-Parallelism of gut mb/analysis/fish/mantel_tests/non-rarified/fish_dm_wuf_skin_nr_phylo.tsv", sep = "\t", dec = ".", row.names = 1, header= TRUE, check.names = FALSE)

rownames(fish_dm_wuf_skin_nr) <- fish_meta_skin$species_name_timetree[match(rownames(fish_dm_wuf_skin_nr), fish_meta_skin$id)]
colnames(fish_dm_wuf_skin_nr) <- fish_meta_skin$species_name_timetree[match(colnames(fish_dm_wuf_skin_nr), fish_meta_skin$id)]

fish_dm_wuf_skin_nr <- fish_dm_wuf_skin_nr[order(rownames(fish_dm_wuf_skin_nr)) , order(colnames(fish_dm_wuf_skin_nr))]

fish_dm_wuf_skin_nr_dist <- as.dist(fish_dm_wuf_skin_nr)
fish_dm_wuf_skin_nr_matrix <- as.matrix(fish_dm_wuf_skin_nr)

#Mantel test

#unweighted unifrac
mantel.test(fish_phylo_skin_divtime, fish_dm_uwuf_skin_nr_matrix, nperm = 999, graph = TRUE)
mantel.rtest(fish_phylo_skin_divtime_dist, fish_dm_uwuf_skin_nr_dist, nrepet = 999)
mantel(fish_phylo_skin_divtime_dist, fish_dm_uwuf_skin_nr_dist, permutations = 999, method = "spearman")
mantel(fish_phylo_skin_divtime_dist, fish_dm_uwuf_skin_nr_dist, permutations = 999, method = "pearson")

#weighted unifrac
mantel.test(fish_phylo_skin_divtime, fish_dm_wuf_skin_nr_matrix, nperm = 999, graph = TRUE)
mantel.rtest(fish_phylo_skin_divtime_dist, fish_dm_wuf_skin_nr_dist, nrepet = 999)
mantel(fish_phylo_skin_divtime_dist, fish_dm_wuf_skin_nr_dist, permutations = 999, method = "spearman")
mantel(fish_phylo_skin_divtime_dist, fish_dm_wuf_skin_nr_dist, permutations = 999, method = "pearson")

#Make data frame in long format from distance matrix
fish_skin_divtime_df <- dist2list(fish_phylo_skin_divtime_dist)
fish_skin_divtime_df$value <- fish_skin_divtime_df$value/2
fish_uwuf_skin_df <- dist2list(fish_dm_uwuf_skin_nr_dist)
fish_wuf_skin_df <- dist2list(fish_dm_wuf_skin_nr_dist)

#Remove rows with comparisons of the same species
fish_skin_divtime_df2 <- fish_skin_divtime_df[fish_skin_divtime_df$value != 0,]
fish_uwuf_skin_df2 <- fish_uwuf_skin_df[fish_uwuf_skin_df$value != 0,]
fish_wuf_skin_df2 <- fish_wuf_skin_df[fish_wuf_skin_df$value != 0,]

plot(fish_uwuf_skin_df2$value ~ fish_skin_divtime_df2$value, ylab = "gut mb distance", xlab = "divergence time", main = "unweighted unifrac-skins", pch = 16)
abline(lm(fish_uwuf_skin_df2$value ~ fish_skin_divtime_df2$value))
summary(lm(fish_uwuf_skin_df2$value ~ fish_skin_divtime_df2$value))

plot(fish_wuf_skin_df2$value ~ fish_skin_divtime_df2$value, ylab = "gut mb distance", xlab = "divergence time", main = "weighted unifrac-skins", pch = 16)
abline(lm(fish_wuf_skin_df2$value ~ fish_skin_divtime_df2$value))
summary(lm(fish_wuf_skin_df2$value ~ fish_skin_divtime_df2$value))
