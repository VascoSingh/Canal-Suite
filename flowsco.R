# Install packages if required
if(!require('flowCore')) {install.packages('flowCore')}
if(!require('Biobase')) {install.packages('Biobase')}

# Load packages
library('flowCore')
library('Biobase')

# Set working directory
getwd()
setwd("Home/Documents/GitHib/BaxterAlgorithms/VS_Functions") # set to your desired working folder (directory)
PrimaryDirectory <- getwd()

# Find file names of .csv files in the current working directory
FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")
FileNames

# Chose which .csv file to read into 'data' -- rename 'sample_data.csv' to whatever file you want to read
data <- read.csv("BigData.csv", row.names = 1) # if the first column contains names for each row, then change "row.names = NULL" to "row.names = 1" 
data

# Give a name for your output .fcs file (don't foret to add '.fcs' at the end of the name)
fcsfilename <- "testdata.fcs"


##### END USER INPUT #####


# Convert data to matrix
data <- as.matrix(data)

# Check data and data column names -- for this script to work, the first row must be the column names
head(data)
dimnames(data)[[2]]

# Create FCS file metadata - column names with descriptions
metadata <- data.frame(name=dimnames(data)[[2]],
                       desc=paste('this is column',dimnames(data)[[2]],'from your CSV') #Change?
)

metadata

# Create FCS file metadata - ranges, min, and max settings
metadata$range <- apply(apply(data,2,range),2,diff)
metadata$minRange <- apply(data,2,min)
metadata$maxRange <- apply(data,2,max)

metadata$range
metadata$minRange
metadata$maxRange

# Create flowframe with tSNE data
data.ff <- new("flowFrame",
               exprs=data,
               parameters=AnnotatedDataFrame(metadata)
)

data.ff

# Create an 'output' folder
setwd(PrimaryDirectory)
dir.create("Output", showWarnings = FALSE)
setwd("Output")

# Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
write.FCS(data.ff, fcsfilename)

setwd(PrimaryDirectory)



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#Analyzing flow data

#Load the packages
library(flowCore)
library(flowAI)
library(ggcyto)
library(flowWorkspace)
library(openCyto)
library(gridExtra)

#Change to path of your output file
myfile <- "/Users/vasco/Documents/GitHub/BaxterAlgorithms/VS_Functions/Output/testdata.fcs"
library(flowCore)
fcsfile <- read.FCS(myfile)
fcsfile
names(fcsfile) #View the parameters of your data
keyword(fcsfile) #View the metadata of flowFrame

#Cleaning
fcsfile_clean <- flow_auto_qc(fcsfile)
fcsfile_clean

#Transformation
trans<-estimateLogicle(fcsfile_clean, colnames(fcsfile_clean[,1:12]))
fcsfile_clean_trans <- transform(fcsfile_clean, trans)

#Visualize the Results
autoplot(fcsfile_clean_trans)
autoplot(fcsfile_clean_trans, x="Area", y="MinIntensity", bin = 64)

# [1] "<Area> Area"                   "<Centroid_1> Centroid_1"       "<Centroid_2> Centroid_2"      
# [4] "<BoundingBox_1> BoundingBox_1" "<BoundingBox_2> BoundingBox_2" "<BoundingBox_3> BoundingBox_3"
# [7] "<BoundingBox_4> BoundingBox_4" "<EquivDiameter> EquivDiameter" "<Extent> Extent"              
# [10] "<MeanIntensity> MeanIntensity" "<MinIntensity> MinIntensity"   "<MaxIntensity> MaxIntensity"  
# [13] "<WellNum> WellNum"             "<TimeNum> TimeNum"             "<Cell> Cell"                  
# [16] "<AnaPass> AnaPass"             "<ImgPlane> ImgPlane"          


#Basic Gating
ggcyto(fcsfile_clean_trans, aes(x="Area", y="MinIntensity"))+geom_hex(bins=64);

install.packages('flowWorkspace')

rg1<-rectangleGate("Area"=c(0.6, Inf), filterId = "NoneDebris")
ggcyto(fcsfile_clean_trans, aes(x="Area", y="MinIntensity"))+geom_hex(bins=64)+geom_gate(rg1);
rg1


#GGCYTO -> Flow Cytometry Visualization
p<-ggcyto(fcsfile_clean_trans, aes(x="Area", y="MinIntensity"), subset="singlets")
p<-p + geom_hex(bins=64)


