##R notebook for alpha diversity analysis of FMP

library(tidyverse)
library(ggplot2)
library(gridExtra)
#I. loading in alpha diversity

getwd()
setwd("Desktop/FMP_biom/alpha/")
alpha<-read_csv("MF.alpha_2021-09-20.csv")
head(alpha)

#subset by rows...create new df for each gill, skin, and digesta
gill.df<-filter(alpha,sample_type=="fish gill")
skin.df<-filter(alpha,sample_type=="fish skin")
mg.df<-filter(alpha,sample_type=="fish midgut digesta")
hg.df<-filter(alpha,sample_type=="fish hindgut digesta")

#geom boxplot reference
#http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization

hab1_g<-ggplot(gill.df, aes(habitat_depth_level1, chao1))+
  geom_boxplot()
hab1_g
kruskal.test(chao1 ~ habitat_depth_level1, data = gill.df)

hab2_g<-ggplot(gill.df, aes(habitat_depth_level2, chao1))+
  geom_boxplot()
hab2_g
kruskal.test(chao1 ~ habitat_depth_level2, data = gill.df)

clim_g<-ggplot(gill.df, aes(climate, chao1))+
  geom_boxplot()
clim_g
kruskal.test(chao1 ~ climate, data = gill.df)

collsub_g<-ggplot(gill.df, aes(collection_substrate, chao1))+
  geom_boxplot()
collsub_g
kruskal.test(chao1 ~ collection_substrate, data = gill.df)

subgroup_g<-ggplot(gill.df, aes(substrata_group, chao1))+
  geom_boxplot()
subgroup_g
kruskal.test(chao1~substrata_group, data = gill.df)

subgroup_s<-ggplot(skin.df, aes(substrata_group, chao1))+
  geom_boxplot()
subgroup_s

subgroup_mg<-ggplot(mg.df, aes(substrata_group, chao1))+
  geom_boxplot()
subgroup_mg

subgroup_hg<-ggplot(hg.df, aes(substrata_group, chao1))+
  geom_boxplot()
subgroup_hg

salinity_g<-ggplot(gill.df, aes(salinity_tolerance, chao1))+
  geom_boxplot()
salinity_g
kruskal.test(chao1~salinity_tolerance, data = gill.df)

swimperf_g<-ggplot(gill.df, aes(swim_performance, chao1))+
  geom_boxplot()
swimperf_g
kruskal.test(chao1~swim_performance, data = gill.df)

swimmode_g<-ggplot(gill.df, aes(swim_mode, chao1))+
  geom_boxplot()
swimmode_g
kruskal.test(chao1~swim_mode, data = gill.df)

## SKIN

hab1_s<-ggplot(skin.df, aes(habitat_depth_level1, chao1))+
  geom_boxplot()
hab1_g

hab2_g<-ggplot(skin.df, aes(habitat_depth_level2, chao1))+
  geom_boxplot()
hab2_g

clim_g<-ggplot(skin.df, aes(climate, chao1))+
  geom_boxplot()
clim_g

collsub_g<-ggplot(skin.df, aes(collection_substrate, chao1))+
  geom_boxplot()
collsub_g

subgroup_g<-ggplot(skin.df, aes(substrata_group, chao1))+
  geom_boxplot()
subgroup_g

salinity_g<-ggplot(skin.df, aes(salinity_tolerance, chao1))+
  geom_boxplot()
salinity_g

swimperf_g<-ggplot(skin.df, aes(swim_performance, chao1))+
  geom_boxplot()
swimperf_g

swimmode_g<-ggplot(skin.df, aes(swim_mode, chao1))+
  geom_boxplot()
swimmode_g



#Stat test of categorical chao1
# Gill
kruskal.test(chao1 ~ habitat_depth_level1, data = gill.df)
kruskal.test(chao1 ~ habitat_depth_level2, data = gill.df)
kruskal.test(chao1 ~ climate, data = gill.df)
kruskal.test(chao1 ~ collection_substrate, data = gill.df)
kruskal.test(chao1~substrata_group, data = gill.df)
kruskal.test(chao1~type_habitat_nonbay_bol, data = gill.df)
kruskal.test(chao1~salinity_tolerance, data = gill.df)
kruskal.test(chao1~swim_performance, data = gill.df)
kruskal.test(chao1~swim_mode, data = gill.df)
kruskal.test(chao1~p_class, data = gill.df)

#Spearman correlation
cor.test(gill.df$trophic_likely, gill.df$chao1, method='spearman')
cor.test(gill.df$biomass_even, gill.df$chao1, method='spearman')
cor.test(gill.df$biomass_quartile, gill.df$chao1, method='spearman')
cor.test(gill.df$swim_acceleration, gill.df$chao1, method='spearman')
cor.test(gill.df$swim_endurance, gill.df$chao1, method='spearman')
cor.test(gill.df$ratio_gi_to_tl, gill.df$chao1, method='spearman')
cor.test(gill.df$ratio_dorsal_to_tl, gill.df$chao1, method='spearman')
cor.test(gill.df$ratio_gape_to_tl, gill.df$chao1, method='spearman')
cor.test(gill.df$ph_cells_per_g_log, gill.df$chao1, method='spearman')

# SKIN
kruskal.test(chao1 ~ habitat_depth_level1, data = skin.df)
kruskal.test(chao1 ~ habitat_depth_level2, data = skin.df)
kruskal.test(chao1 ~ climate, data = skin.df)
kruskal.test(chao1 ~ collection_substrate, data = skin.df)
kruskal.test(chao1~substrata_group, data = skin.df)
kruskal.test(chao1~type_habitat_nonbay_bol, data = skin.df)
kruskal.test(chao1~salinity_tolerance, data = skin.df)
kruskal.test(chao1~swim_performance, data = skin.df)
kruskal.test(chao1~swim_mode, data = skin.df)
kruskal.test(chao1~p_class, data = skin.df)

#chao
cor.test(skin.df$trophic_likely, skin.df$chao1, method='spearman')
cor.test(skin.df$biomass_even, skin.df$chao1, method='spearman')
cor.test(skin.df$biomass_quartile, skin.df$chao1, method='spearman')
cor.test(skin.df$swim_acceleration, skin.df$chao1, method='spearman')
cor.test(skin.df$swim_endurance, skin.df$chao1, method='spearman')
cor.test(skin.df$ratio_gi_to_tl, skin.df$chao1, method='spearman')
cor.test(skin.df$ratio_dorsal_to_tl, skin.df$chao1, method='spearman')
cor.test(skin.df$ratio_gape_to_tl, skin.df$chao1, method='spearman')
cor.test(skin.df$ph_cells_per_g_log, skin.df$chao1, method='spearman')

#MG
kruskal.test(chao1 ~ habitat_depth_level1, data = mg.df)
kruskal.test(chao1 ~ habitat_depth_level2, data = mg.df)
kruskal.test(chao1 ~ climate, data = mg.df)
kruskal.test(chao1 ~ collection_substrate, data = mg.df)
kruskal.test(chao1~substrata_group, data = mg.df)
kruskal.test(chao1~type_habitat_nonbay_bol, data = mg.df)
kruskal.test(chao1~salinity_tolerance, data = mg.df)
kruskal.test(chao1~swim_performance, data = mg.df)
kruskal.test(chao1~swim_mode, data = mg.df)
kruskal.test(chao1~p_class, data = mg.df)

#chao
cor.test(mg.df$trophic_likely, mg.df$chao1, method='spearman')
cor.test(mg.df$biomass_even, mg.df$chao1, method='spearman')
cor.test(mg.df$biomass_quartile, mg.df$chao1, method='spearman')
cor.test(mg.df$swim_acceleration, mg.df$chao1, method='spearman')
cor.test(mg.df$swim_endurance, mg.df$chao1, method='spearman')
cor.test(mg.df$ratio_gi_to_tl, mg.df$chao1, method='spearman')
cor.test(mg.df$ratio_dorsal_to_tl, mg.df$chao1, method='spearman')
cor.test(mg.df$ratio_gape_to_tl, mg.df$chao1, method='spearman')
cor.test(mg.df$ph_cells_per_g_log, mg.df$chao1, method='spearman')

#HG
kruskal.test(chao1 ~ habitat_depth_level1, data = hg.df)
kruskal.test(chao1 ~ habitat_depth_level2, data = hg.df)
kruskal.test(chao1 ~ climate, data = hg.df)
kruskal.test(chao1 ~ collection_substrate, data = hg.df)
kruskal.test(chao1~substrata_group, data = hg.df)
kruskal.test(chao1~type_habitat_nonbay_bol, data = hg.df)
kruskal.test(chao1~salinity_tolerance, data = hg.df)
kruskal.test(chao1~swim_performance, data = hg.df)
kruskal.test(chao1~swim_mode, data = hg.df)
kruskal.test(chao1~p_class, data = hg.df)

#chao
cor.test(hg.df$trophic_likely, hg.df$chao1, method='spearman')
cor.test(hg.df$biomass_even, hg.df$chao1, method='spearman')
cor.test(hg.df$biomass_quartile, hg.df$chao1, method='spearman')
cor.test(hg.df$swim_acceleration, hg.df$chao1, method='spearman')
cor.test(hg.df$swim_endurance, hg.df$chao1, method='spearman')
cor.test(hg.df$ratio_gi_to_tl, hg.df$chao1, method='spearman')
cor.test(hg.df$ratio_dorsal_to_tl, hg.df$chao1, method='spearman')
cor.test(hg.df$ratio_gape_to_tl, hg.df$chao1, method='spearman')
cor.test(hg.df$ph_cells_per_g_log, hg.df$chao1, method='spearman')

##Stat test of categorical Faith PD
# Gill
kruskal.test(PD_whole_tree ~ habitat_depth_level1, data = gill.df)
kruskal.test(PD_whole_tree ~ habitat_depth_level2, data = gill.df)
kruskal.test(PD_whole_tree ~ climate, data = gill.df)
kruskal.test(PD_whole_tree ~ collection_substrate, data = gill.df)
kruskal.test(PD_whole_tree~substrata_group, data = gill.df)
kruskal.test(PD_whole_tree~type_habitat_nonbay_bol, data = gill.df)
kruskal.test(PD_whole_tree~salinity_tolerance, data = gill.df)
kruskal.test(PD_whole_tree~swim_performance, data = gill.df)
kruskal.test(PD_whole_tree~swim_mode, data = gill.df)
kruskal.test(PD_whole_tree~p_class, data = gill.df)

#Spearman correlation
cor.test(gill.df$trophic_likely, gill.df$PD_whole_tree, method='spearman')
cor.test(gill.df$biomass_even, gill.df$PD_whole_tree, method='spearman')
cor.test(gill.df$biomass_quartile, gill.df$PD_whole_tree, method='spearman')
cor.test(gill.df$swim_acceleration, gill.df$PD_whole_tree, method='spearman')
cor.test(gill.df$swim_endurance, gill.df$PD_whole_tree, method='spearman')
cor.test(gill.df$ratio_gi_to_tl, gill.df$PD_whole_tree, method='spearman')
cor.test(gill.df$ratio_dorsal_to_tl, gill.df$PD_whole_tree, method='spearman')
cor.test(gill.df$ratio_gape_to_tl, gill.df$PD_whole_tree, method='spearman')
cor.test(gill.df$ph_cells_per_g_log, gill.df$PD_whole_tree, method='spearman')

# SKIN
kruskal.test(PD_whole_tree ~ habitat_depth_level1, data = skin.df)
kruskal.test(PD_whole_tree ~ habitat_depth_level2, data = skin.df)
kruskal.test(PD_whole_tree ~ climate, data = skin.df)
kruskal.test(PD_whole_tree ~ collection_substrate, data = skin.df)
kruskal.test(PD_whole_tree~substrata_group, data = skin.df)
kruskal.test(PD_whole_tree~type_habitat_nonbay_bol, data = skin.df)
kruskal.test(PD_whole_tree~salinity_tolerance, data = skin.df)
kruskal.test(PD_whole_tree~swim_performance, data = skin.df)
kruskal.test(PD_whole_tree~swim_mode, data = skin.df)
kruskal.test(PD_whole_tree~p_class, data = skin.df)

#faiths
cor.test(skin.df$trophic_likely, skin.df$PD_whole_tree, method='spearman')
cor.test(skin.df$biomass_even, skin.df$PD_whole_tree, method='spearman')
cor.test(skin.df$biomass_quartile, skin.df$PD_whole_tree, method='spearman')
cor.test(skin.df$swim_acceleration, skin.df$PD_whole_tree, method='spearman')
cor.test(skin.df$swim_endurance, skin.df$PD_whole_tree, method='spearman')
cor.test(skin.df$ratio_gi_to_tl, skin.df$PD_whole_tree, method='spearman')
cor.test(skin.df$ratio_dorsal_to_tl, skin.df$PD_whole_tree, method='spearman')
cor.test(skin.df$ratio_gape_to_tl, skin.df$PD_whole_tree, method='spearman')
cor.test(skin.df$ph_cells_per_g_log, skin.df$PD_whole_tree, method='spearman')

#MG
kruskal.test(PD_whole_tree ~ habitat_depth_level1, data = mg.df)
kruskal.test(PD_whole_tree ~ habitat_depth_level2, data = mg.df)
kruskal.test(PD_whole_tree ~ climate, data = mg.df)
kruskal.test(PD_whole_tree ~ collection_substrate, data = mg.df)
kruskal.test(PD_whole_tree~substrata_group, data = mg.df)
kruskal.test(PD_whole_tree~type_habitat_nonbay_bol, data = mg.df)
kruskal.test(PD_whole_tree~salinity_tolerance, data = mg.df)
kruskal.test(PD_whole_tree~swim_performance, data = mg.df)
kruskal.test(PD_whole_tree~swim_mode, data = mg.df)
kruskal.test(PD_whole_tree~p_class, data = mg.df)

#faiths
cor.test(mg.df$trophic_likely, mg.df$PD_whole_tree, method='spearman')
cor.test(mg.df$biomass_even, mg.df$PD_whole_tree, method='spearman')
cor.test(mg.df$biomass_quartile, mg.df$PD_whole_tree, method='spearman')
cor.test(mg.df$swim_acceleration, mg.df$PD_whole_tree, method='spearman')
cor.test(mg.df$swim_endurance, mg.df$PD_whole_tree, method='spearman')
cor.test(mg.df$ratio_gi_to_tl, mg.df$PD_whole_tree, method='spearman')
cor.test(mg.df$ratio_dorsal_to_tl, mg.df$PD_whole_tree, method='spearman')
cor.test(mg.df$ratio_gape_to_tl, mg.df$PD_whole_tree, method='spearman')
cor.test(mg.df$ph_cells_per_g_log, mg.df$PD_whole_tree, method='spearman')

#HG
kruskal.test(PD_whole_tree ~ habitat_depth_level1, data = hg.df)
kruskal.test(PD_whole_tree ~ habitat_depth_level2, data = hg.df)
kruskal.test(PD_whole_tree ~ climate, data = hg.df)
kruskal.test(PD_whole_tree ~ collection_substrate, data = hg.df)
kruskal.test(PD_whole_tree~substrata_group, data = hg.df)
kruskal.test(PD_whole_tree~type_habitat_nonbay_bol, data = hg.df)
kruskal.test(PD_whole_tree~salinity_tolerance, data = hg.df)
kruskal.test(PD_whole_tree~swim_performance, data = hg.df)
kruskal.test(PD_whole_tree~swim_mode, data = hg.df)
kruskal.test(PD_whole_tree~p_class, data = hg.df)

#faiths
cor.test(hg.df$trophic_likely, hg.df$PD_whole_tree, method='spearman')
cor.test(hg.df$biomass_even, hg.df$PD_whole_tree, method='spearman')
cor.test(hg.df$biomass_quartile, hg.df$PD_whole_tree, method='spearman')
cor.test(hg.df$swim_acceleration, hg.df$PD_whole_tree, method='spearman')
cor.test(hg.df$swim_endurance, hg.df$PD_whole_tree, method='spearman')
cor.test(hg.df$ratio_gi_to_tl, hg.df$PD_whole_tree, method='spearman')
cor.test(hg.df$ratio_dorsal_to_tl, hg.df$PD_whole_tree, method='spearman')
cor.test(hg.df$ratio_gape_to_tl, hg.df$PD_whole_tree, method='spearman')
cor.test(hg.df$ph_cells_per_g_log, hg.df$PD_whole_tree, method='spearman')

##stats test of categorical cell density (log) md==ph_cells_per_g_log

# Gill
kruskal.test(ph_cells_per_g_log ~ habitat_depth_level1, data = gill.df)
kruskal.test(ph_cells_per_g_log ~ habitat_depth_level2, data = gill.df)
kruskal.test(ph_cells_per_g_log ~ climate, data = gill.df)
kruskal.test(ph_cells_per_g_log ~ collection_substrate, data = gill.df)
kruskal.test(ph_cells_per_g_log~substrata_group, data = gill.df)
kruskal.test(ph_cells_per_g_log~type_habitat_nonbay_bol, data = gill.df)
kruskal.test(ph_cells_per_g_log~salinity_tolerance, data = gill.df)
kruskal.test(ph_cells_per_g_log~swim_performance, data = gill.df)
kruskal.test(ph_cells_per_g_log~swim_mode, data = gill.df)
kruskal.test(ph_cells_per_g_log~p_class, data = gill.df)

swimmode_g<-ggplot(gill.df, aes(swim_mode, h_cells_per_g_log))+
  geom_boxplot()
swimmode_g

#Spearman correlation
cor.test(gill.df$trophic_likely, gill.df$ph_cells_per_g_log, method='spearman')
cor.test(gill.df$biomass_even, gill.df$ph_cells_per_g_log, method='spearman')
cor.test(gill.df$biomass_quartile, gill.df$ph_cells_per_g_log, method='spearman')
cor.test(gill.df$swim_acceleration, gill.df$ph_cells_per_g_log, method='spearman')
cor.test(gill.df$swim_endurance, gill.df$ph_cells_per_g_log, method='spearman')
cor.test(gill.df$ratio_gi_to_tl, gill.df$ph_cells_per_g_log, method='spearman')
cor.test(gill.df$ratio_dorsal_to_tl, gill.df$ph_cells_per_g_log, method='spearman')
cor.test(gill.df$ratio_gape_to_tl, gill.df$ph_cells_per_g_log, method='spearman')
cor.test(gill.df$ph_cells_per_g_log, gill.df$ph_cells_per_g_log, method='spearman')

# SKIN
kruskal.test(ph_cells_per_g_log ~ habitat_depth_level1, data = skin.df)
kruskal.test(ph_cells_per_g_log ~ habitat_depth_level2, data = skin.df)
kruskal.test(ph_cells_per_g_log ~ climate, data = skin.df)
kruskal.test(ph_cells_per_g_log ~ collection_substrate, data = skin.df)
kruskal.test(ph_cells_per_g_log~substrata_group, data = skin.df)
kruskal.test(ph_cells_per_g_log~type_habitat_nonbay_bol, data = skin.df)
kruskal.test(ph_cells_per_g_log~salinity_tolerance, data = skin.df)
kruskal.test(ph_cells_per_g_log~swim_performance, data = skin.df)
kruskal.test(ph_cells_per_g_log~swim_mode, data = skin.df)
kruskal.test(ph_cells_per_g_log~p_class, data = skin.df)

#cells
cor.test(skin.df$trophic_likely, skin.df$ph_cells_per_g_log, method='spearman')
cor.test(skin.df$biomass_even, skin.df$ph_cells_per_g_log, method='spearman')
cor.test(skin.df$biomass_quartile, skin.df$ph_cells_per_g_log, method='spearman')
cor.test(skin.df$swim_acceleration, skin.df$ph_cells_per_g_log, method='spearman')
cor.test(skin.df$swim_endurance, skin.df$ph_cells_per_g_log, method='spearman')
cor.test(skin.df$ratio_gi_to_tl, skin.df$ph_cells_per_g_log, method='spearman')
cor.test(skin.df$ratio_dorsal_to_tl, skin.df$ph_cells_per_g_log, method='spearman')
cor.test(skin.df$ratio_gape_to_tl, skin.df$ph_cells_per_g_log, method='spearman')
cor.test(skin.df$ph_cells_per_g_log, skin.df$ph_cells_per_g_log, method='spearman')


#MG
kruskal.test(ph_cells_per_g_log ~ habitat_depth_level1, data = mg.df)
kruskal.test(ph_cells_per_g_log ~ habitat_depth_level2, data = mg.df)
kruskal.test(ph_cells_per_g_log ~ climate, data = mg.df)
kruskal.test(ph_cells_per_g_log ~ collection_substrate, data = mg.df)
kruskal.test(ph_cells_per_g_log~substrata_group, data = mg.df)
kruskal.test(ph_cells_per_g_log~type_habitat_nonbay_bol, data = mg.df)
kruskal.test(ph_cells_per_g_log~salinity_tolerance, data = mg.df)
kruskal.test(ph_cells_per_g_log~swim_performance, data = mg.df)
kruskal.test(ph_cells_per_g_log~swim_mode, data = mg.df)
kruskal.test(ph_cells_per_g_log~p_class, data = mg.df)

#cells
cor.test(mg.df$trophic_likely, mg.df$ph_cells_per_g_log, method='spearman')
cor.test(mg.df$biomass_even, mg.df$ph_cells_per_g_log, method='spearman')
cor.test(mg.df$biomass_quartile, mg.df$ph_cells_per_g_log, method='spearman')
cor.test(mg.df$swim_acceleration, mg.df$ph_cells_per_g_log, method='spearman')
cor.test(mg.df$swim_endurance, mg.df$ph_cells_per_g_log, method='spearman')
cor.test(mg.df$ratio_gi_to_tl, mg.df$ph_cells_per_g_log, method='spearman')
cor.test(mg.df$ratio_dorsal_to_tl, mg.df$ph_cells_per_g_log, method='spearman')
cor.test(mg.df$ratio_gape_to_tl, mg.df$ph_cells_per_g_log, method='spearman')
cor.test(mg.df$ph_cells_per_g_log, mg.df$ph_cells_per_g_log, method='spearman')

#HG
kruskal.test(ph_cells_per_g_log ~ habitat_depth_level1, data = hg.df)
kruskal.test(ph_cells_per_g_log ~ habitat_depth_level2, data = hg.df)
kruskal.test(ph_cells_per_g_log ~ climate, data = hg.df)
kruskal.test(ph_cells_per_g_log ~ collection_substrate, data = hg.df)
kruskal.test(ph_cells_per_g_log~substrata_group, data = hg.df)
kruskal.test(ph_cells_per_g_log~type_habitat_nonbay_bol, data = hg.df)
kruskal.test(ph_cells_per_g_log~salinity_tolerance, data = hg.df)
kruskal.test(ph_cells_per_g_log~swim_performance, data = hg.df)
kruskal.test(ph_cells_per_g_log~swim_mode, data = hg.df)
kruskal.test(ph_cells_per_g_log~p_class, data = hg.df)

#cells
cor.test(hg.df$trophic_likely, hg.df$ph_cells_per_g_log, method='spearman')
cor.test(hg.df$biomass_even, hg.df$ph_cells_per_g_log, method='spearman')
cor.test(hg.df$biomass_quartile, hg.df$ph_cells_per_g_log, method='spearman')
cor.test(hg.df$swim_acceleration, hg.df$ph_cells_per_g_log, method='spearman')
cor.test(hg.df$swim_endurance, hg.df$ph_cells_per_g_log, method='spearman')
cor.test(hg.df$ratio_gi_to_tl, hg.df$ph_cells_per_g_log, method='spearman')
cor.test(hg.df$ratio_dorsal_to_tl, hg.df$ph_cells_per_g_log, method='spearman')
cor.test(hg.df$ratio_gape_to_tl, hg.df$ph_cells_per_g_log, method='spearman')
cor.test(hg.df$ph_cells_per_g_log, hg.df$ph_cells_per_g_log, method='spearman')


##Shannon analysis
# Gill
kruskal.test(shannon ~ habitat_depth_level1, data = gill.df)
kruskal.test(shannon ~ habitat_depth_level2, data = gill.df)
kruskal.test(shannon ~ climate, data = gill.df)
kruskal.test(shannon ~ collection_substrate, data = gill.df)
kruskal.test(shannon~substrata_group, data = gill.df)
kruskal.test(shannon~type_habitat_nonbay_bol, data = gill.df)
kruskal.test(shannon~salinity_tolerance, data = gill.df)
kruskal.test(shannon~swim_performance, data = gill.df)
kruskal.test(shannon~swim_mode, data = gill.df)
kruskal.test(shannon~p_class, data = gill.df)

#Spearman correlation
cor.test(gill.df$trophic_likely, gill.df$shannon, method='spearman')
cor.test(gill.df$biomass_even, gill.df$shannon, method='spearman')
cor.test(gill.df$biomass_quartile, gill.df$shannon, method='spearman')
cor.test(gill.df$swim_acceleration, gill.df$shannon, method='spearman')
cor.test(gill.df$swim_endurance, gill.df$shannon, method='spearman')
cor.test(gill.df$ratio_gi_to_tl, gill.df$shannon, method='spearman')
cor.test(gill.df$ratio_dorsal_to_tl, gill.df$shannon, method='spearman')
cor.test(gill.df$ratio_gape_to_tl, gill.df$shannon, method='spearman')
cor.test(gill.df$ph_cells_per_g_log, gill.df$shannon, method='spearman')

# SKIN
kruskal.test(shannon ~ habitat_depth_level1, data = skin.df)
kruskal.test(shannon ~ habitat_depth_level2, data = skin.df)
kruskal.test(shannon ~ climate, data = skin.df)
kruskal.test(shannon ~ collection_substrate, data = skin.df)
kruskal.test(shannon~substrata_group, data = skin.df)
kruskal.test(shannon~type_habitat_nonbay_bol, data = skin.df)
kruskal.test(shannon~salinity_tolerance, data = skin.df)
kruskal.test(shannon~swim_performance, data = skin.df)
kruskal.test(shannon~swim_mode, data = skin.df)
kruskal.test(shannon~p_class, data = skin.df)

#shann
cor.test(skin.df$trophic_likely, skin.df$shannon, method='spearman')
cor.test(skin.df$biomass_even, skin.df$shannon, method='spearman')
cor.test(skin.df$biomass_quartile, skin.df$shannon, method='spearman')
cor.test(skin.df$swim_acceleration, skin.df$shannon, method='spearman')
cor.test(skin.df$swim_endurance, skin.df$shannon, method='spearman')
cor.test(skin.df$ratio_gi_to_tl, skin.df$shannon, method='spearman')
cor.test(skin.df$ratio_dorsal_to_tl, skin.df$shannon, method='spearman')
cor.test(skin.df$ratio_gape_to_tl, skin.df$shannon, method='spearman')
cor.test(skin.df$ph_cells_per_g_log, skin.df$shannon, method='spearman')

#MG
kruskal.test(shannon ~ habitat_depth_level1, data = mg.df)
kruskal.test(shannon ~ habitat_depth_level2, data = mg.df)
kruskal.test(shannon ~ climate, data = mg.df)
kruskal.test(shannon ~ collection_substrate, data = mg.df)
kruskal.test(shannon~substrata_group, data = mg.df)
kruskal.test(shannon~type_habitat_nonbay_bol, data = mg.df)
kruskal.test(shannon~salinity_tolerance, data = mg.df)
kruskal.test(shannon~swim_performance, data = mg.df)
kruskal.test(shannon~swim_mode, data = mg.df)
kruskal.test(shannon~p_class, data = mg.df)

#shann
cor.test(mg.df$trophic_likely, mg.df$shannon, method='spearman')
cor.test(mg.df$biomass_even, mg.df$shannon, method='spearman')
cor.test(mg.df$biomass_quartile, mg.df$shannon, method='spearman')
cor.test(mg.df$swim_acceleration, mg.df$shannon, method='spearman')
cor.test(mg.df$swim_endurance, mg.df$shannon, method='spearman')
cor.test(mg.df$ratio_gi_to_tl, mg.df$shannon, method='spearman')
cor.test(mg.df$ratio_dorsal_to_tl, mg.df$shannon, method='spearman')
cor.test(mg.df$ratio_gape_to_tl, mg.df$shannon, method='spearman')
cor.test(mg.df$ph_cells_per_g_log, mg.df$shannon, method='spearman')

#HG
kruskal.test(shannon ~ habitat_depth_level1, data = hg.df)
kruskal.test(shannon ~ habitat_depth_level2, data = hg.df)
kruskal.test(shannon ~ climate, data = hg.df)
kruskal.test(shannon ~ collection_substrate, data = hg.df)
kruskal.test(shannon~substrata_group, data = hg.df)
kruskal.test(shannon~type_habitat_nonbay_bol, data = hg.df)
kruskal.test(shannon~salinity_tolerance, data = hg.df)
kruskal.test(shannon~swim_performance, data = hg.df)
kruskal.test(shannon~swim_mode, data = hg.df)
kruskal.test(shannon~p_class, data = hg.df)

#shann
cor.test(hg.df$trophic_likely, hg.df$shannon, method='spearman')
cor.test(hg.df$biomass_even, hg.df$shannon, method='spearman')
cor.test(hg.df$biomass_quartile, hg.df$shannon, method='spearman')
cor.test(hg.df$swim_acceleration, hg.df$shannon, method='spearman')
cor.test(hg.df$swim_endurance, hg.df$shannon, method='spearman')
cor.test(hg.df$ratio_gi_to_tl, hg.df$shannon, method='spearman')
cor.test(hg.df$ratio_dorsal_to_tl, hg.df$shannon, method='spearman')
cor.test(hg.df$ratio_gape_to_tl, hg.df$shannon, method='spearman')
cor.test(hg.df$ph_cells_per_g_log, hg.df$shannon, method='spearman')



##numeric variables
#linear model: 
# lmH<-lm(target_variable~predictor_variable, data=df.gill)
#summary(lmH)

#KW tetst
#kruskal.test(weight ~ group, data = my_data)
#post hoc testing...
#pairwise.wilcox.test(PlantGrowth$weight, PlantGrowth$group,p.adjust.method = "BH")

ggplot(gill.df, aes(x=biomass_quartile, y=chao1))+
  geom_point()+
  geom_smooth(method=lm)

lmTrophic_g<-lm(trophic_likely~chao1, data=gill.df)
summary(lmTrophic_g)

lmSA_g<-lm(swim_acceleration~chao1, data=gill.df)
summary(lmSA_g)

lmSE_g<-lm(swim_endurance~chao1, data=gill.df)
summary(lmSE_g)

lmFL_g<-lm(fl_cm~chao1,data=gill.df)
summary(lmFL_g)

lmGape_g<-lm(gape_cm~chao1, data=gill.df)
summary(lmGape_g)

lmGITL_g<-lm(ratio_gi_to_tl~chao1, data=gill.df)
summary(lmGITL_g)

lmDTL_g<-lm(ratio_dorsal_to_tl~chao1, data=gill.df)
summary(lmDTL_g)

lmGaTL_g<-lm(ratio_gape_to_tl~chao1, data=gill.df)
summary(lmGaTL_g)

lmK_g<-lm(condition_factor~chao1, data=gill.df)
summary(lmK_g)

lmLog_g<-lm(ph_cells_per_g_log~chao1, data=gill.df)
summary(lmLog_g)

#not using
lmBQ_g<-lm(biomass_quartile~chao1, data=gill.df)
summary(lmBQ_g)

lmBE_g<-lm(biomass_even~chao1, data=gill.df)
summary(lmBE_g)






#changes columns from string into numeric
#gill.df2$chlorophyll<-as.numeric(gill.df2$chlorophyll)


##multifigure plot 
# reference: https://michaelgastner.com/R_for_QR/multi-panel-plots.html
#fig1<-
#fig2<-
#par(mfrow=c(5,5))