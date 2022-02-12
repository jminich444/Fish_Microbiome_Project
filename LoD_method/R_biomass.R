###Notebook for calculating microbial biomass from reads when using standards
## Step I: Katharoseq -- do not currently have a notebook for this other than jupyter notebook created by Lisa/Pedro
  # Step I. determine the lower read count to exclude samples "original katharoseq"
  #Step 1) Go to Qiita analysis ID 43679 ()

  #Step 2) filter samples which were all run using the same methods 
    #a) samples from your actual study:  title='FMP101'

  #Step 3) In Qiita, annotate features using the  (pre-fitted sklearn-based taxonomy classifier) use silva for marine samples [classify_sklearn] and run command <visualize taxonomy with an interactive barplot> [barplot]

  #Step 4) download the level-7 csv file

  #Step 5) sum the total number of reads per sample and save into a new column named <total_reads> along with creating a new column called <log_total_reads> which is the log of the <total_reads>

  #Step 6) filter your control samples which include known cell counts  

    #b) samples which are positive controls:  control_2='3'
      # data dict: control_2=0 'actual samples'; 1= 'dna extraction negative controls where water is added'; 2='dna extraction negative controls where are specific for sampling e.g a swab'; 3='dna extraction positive controls where cells of isolates or mock communities are used'; 4='library prep negative controls'; 5='library prep positive controls DNA titrations'
      #c) exclude zymo mock community standards as they do not have known cell counts 
      #control_type='KL mock'

  #step 7) manually identify the sOTUs/ASVs which are part of the controls 'targets' (known strains or mock community isolates)
      # in this case this is a Paracoccus and Bacillus mock

  #step 8) sum the total reads per the control targets and save as new column name <sum_controls>

  #step 9) find the relative abundance of the controls by taking sum_controls/total_reads and save as new column <control_rel_abund>
  
  #step 10) plot total_reads (x axis) by the control_rel_abund and fit with allosteric sigmoidal (using Prism or notebook)

  #step 11) use wolfram or an equation solver to solve for x when y=0.9 # this is a personal preference and in original katharoseq we used 0.5 but since I've refined to use 0.9 - this would indicate: for a control sample, the number of reads where at least 90% of the relative abundance is made up of the expected control targets (bacillus and paracoccus)
      # When y=0.9, x=1150 reads in this example

  #step 12) for all samples used in the given analysis and processed the same as controls, exclude any sample with less than 1150 reads

  #step 13) exclude any positive control with less than 1150 reads

########################################################################
## Step II. Estimate the microbial biomass of samples by using the read counts from controls with known cellular inputs

##key metadata columns
#total_reads
#control_cell_into_extraction [total amount of cells estimated to be put into DNA extraction tube]
#control_2 [3=DNA positive controls]
#control_type [KL_mock = standards with known cell counts vs zymo = mock community only useful for composition]

library(tidyverse)
library(ggplot2)
getwd()

#navigate to folder containing file with sample names, metadata, and a column 'total_reads' which is the total ASVs per sample
setwd("Desktop/LoD_method/")

#read in the file containing metadata with total reads per sample
biomass<-read.csv("FMP_reads.csv")
head(biomass)

    glimpse(biomass) #gives overview of variable types in columns
    dim(biomass) #gives the shape of the object 'matrix or df'

#remove samples with reads less than the desired cutoff as calculated previously (eg. 1150)
biomass %>% filter(total_reads > 1150)

dim(biomass)

#log10 transform the total_reads and store in new column <log_total_reads>
biomass$log_total_reads<- log10(biomass$total_reads)

#subset by positive controls only...create new df 
control.df<-filter(biomass,control_2=="3")

    #verify subset worked
    testf<-ggplot(control.df, aes(control_type,log_total_reads))+
      geom_boxplot()
    testf

#subset by positive controls with known cell counts (KL_mock)
control.df<-filter(control.df,control_type=="KL_mock")

#log transform the total cells into extraction and save as new column <log_control_cell_into_extraction>
control.df$log_control_cell_into_extraction<- log10(control.df$control_cell_into_extraction)

    #verify subset worked
    testf<-ggplot(control.df, aes(log_total_reads,log_control_cell_into_extraction))+
      geom_point() + geom_smooth(method="lm")
    testf

#Run a linear model on the standards log(total_reads) by log (cells) to obtain equation 
lm_standards<-lm(control.df$log_total_reads~control.df$log_control_cell_into_extraction)
summary(lm_standards)

#hack to export slope and y-int 
  #https://davetang.org/muse/2012/02/10/manual-linear-regression-analysis-using-r/
#create a slope function
slope<-function(x,y){
  mean_x<-mean(x)
  mean_y<-mean(y)
  nom<-sum((x-mean_x)*(y-mean_y))
  denom<-sum((x-mean_x)^2)
  m<-nom/denom
  return(m)
}
#slope formula is: covariance (x,y)/variance(x)
slope2<-function(x,y){
  return(cov(x,y)/var(x))
}
intercept<-function(x,y,m){
  b<-mean(y)-(m*mean(x))
  return(b)
}
my_slope<-slope(control.df$log_total_reads,control.df$log_control_cell_into_extraction)
my_intercept<-intercept(control.df$log_total_reads,control.df$log_control_cell_into_extraction,my_slope)

#result of slope of standards
my_slope
#3.500501

#result of y-intercept of standards
my_intercept
#-11.01769

#go back to original file and calculate the log cells per extraction save as new column <estimated_biomass_per_pcrrxn>
#enter the PCR reaction volume
pcr_template_vol<-5
#enter the DNA extraction volume
dna_extract_vol<-60

#estimates number of cells in the total PCR reaction volume (e.g in this case: cells per 5ul)
#save into new column called <estimated_biomass_per_pcrrxn>
biomass$estimated_biomass_per_pcrrxn<-10^((biomass$log_total_reads*my_slope)+my_intercept)

#normalize this volume up to the DNA extraction volume (DNA extraction volume / PCR volume) (cells per 60 ul)
#save into new column called <estimated_biomass_per_dnarxn>
biomass$estimated_biomass_per_dnarxn<-biomass$estimated_biomass_per_pcrrxn*(dna_extract_vol/pcr_template_vol)

#normalize this volume now based on the amount of actual physical tissue or primary material placed into the DNA extraction tube: <extraction_mass_g>
#save as <estimated_cells_per_g>
biomass$estimated_cells_per_g<-biomass$estimated_biomass_per_dnarxn/biomass$extraction_mass_g

#log transform biomass estimate and save as new column <log_estimated_cells_per_g>
biomass$log_estimated_cells_per_g<-log10(biomass$estimated_cells_per_g)

#export file 
write.csv(biomass,"biomass.csv", row.names=FALSE)
