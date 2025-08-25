require("devtools")
devtools::install_github("Resende-Lab/SimpleMating")
require("SimpleMating")
library(SimpleMating)

setwd("C:/Users/ehtis/Downloads/simplemate")

############################### Prepare input files ####################################

## GSNP markers - Genotypic file. I am converting HapMap to numeric, compatible for SimpleMating 

# Read HapMap file
geno_raw <- read.delim("NMRIL_1973_geno_diploid.hmp.txt", stringsAsFactors = FALSE)

# HapMap format reminder:
# Col1 = rs# (marker ID), Col2 = alleles, Col3-11 = metadata, Col12+ = taxa genotypes

# Extract marker names and allele info
marker_names <- geno_raw[, 1]   # marker IDs
allele_info  <- geno_raw[, 2]   # e.g. "A/C"

# Extract genotype matrix (taxa in columns)
geno_matrix_raw <- geno_raw[, -(1:11)]

# Get taxa names from column headers
taxa_names <- colnames(geno_matrix_raw)

# Function to recode genotypes into PopVar coding
recode_geno <- function(geno, major, minor) {
  if (geno == paste0(major, major)) return(0)  # homozygous major
  if (geno == paste0(minor, minor)) return(2)   # homozygous minor
  if (geno %in% c(paste0(major, minor), paste0(minor, major))) return(1) # heterozygous
  return(NA)  # missing or unrecognized
}

# Recode all markers
geno_coded <- matrix(NA, nrow = length(taxa_names), ncol = length(marker_names))
for (j in seq_along(marker_names)) {
  alleles <- strsplit(allele_info[j], "/")[[1]]
  major <- alleles[1]
  minor <- alleles[2]
  geno_coded[, j] <- sapply(geno_matrix_raw[j, ], recode_geno, major = major, minor = minor)
}

# Convert to data frame and add taxa as first column
geno_coded_df <- data.frame(Taxa = taxa_names, geno_coded, stringsAsFactors = FALSE)

# Add first row as marker names
colnames(geno_coded_df) <- c("Taxa", marker_names)

# Save SimpleMating-ready file as CSV
write.csv(geno_coded_df, "geno.csv", row.names = FALSE, quote = FALSE)

# now import the file  
geno <- read.csv("geno.csv", row.names = 1)
# Convert to numeric matrix without losing row or column names
geno <- data.matrix(geno)
str(geno)
geno[1:5, 1:5]


# Phenotypic file. BLUPs

# pheno <- read.delim("pheno.txt", stringsAsFactors = FALSE)
# pheno <- read.delim("pheno.txt", row.names = 1)

# twice Taxa columns are due to SimpleMating requirements
pheno <- read.delim("pheno.txt", header = TRUE); rownames(pheno) <- pheno[[1]]
str(pheno)
head(pheno, 5)


# Genetic Map file. 

# Read your geno file 
geno <- read.csv("geno.csv", stringsAsFactors = FALSE)

# helper: part before underscore (e.g., "SCM002812.1" from "SCM002812.1_177620")
pre_underscore <- sub("^(.*)_.*$", "\\1", marker_names)

# extract the digits just before the dot (e.g., "002812" from "SCM002812.1")
digits_before_dot <- sub("^.*?(\\d+)\\..*$", "\\1", pre_underscore)

# get last two digits of that string (so "002812" -> "12")
last_two <- ifelse(nchar(digits_before_dot) >= 2,
                   substr(digits_before_dot, nchar(digits_before_dot) - 1, nchar(digits_before_dot)),
                   digits_before_dot)

# convert to numeric
last_two_num <- as.integer(last_two)

# automatic offset:
# we assume the smallest code maps to chromosome 1. 
# e.g., if codes are 12,13,... then min=12 -> offset = 11 -> 12-11 = 1
offset <- min(last_two_num, na.rm = TRUE) - 1
chromosome <- last_two_num - offset

# extract position (numeric after underscore). If no underscore -> NA
position <- ifelse(grepl("_", marker_names),
                   as.numeric(sub(".*_", "", marker_names)),
                   NA_real_)

# build dataframe (ensures Marker is a column, not rownames)
map_df <- data.frame(
  Chromosome = chromosome,
  Position = position,
  Marker = marker_names,
  stringsAsFactors = FALSE
)

# quick sanity checks
bad_chr <- which(is.na(map_df$Chromosome) | map_df$Chromosome <= 0)
bad_pos <- which(is.na(map_df$Position))

if (length(bad_chr) > 0) {
  warning("Some markers failed to parse a chromosome. Inspect them:")
  print(head(map_df[bad_chr, ], 10))
}
if (length(bad_pos) > 0) {
  warning("Some markers failed to parse a position. Inspect them:")
  print(head(map_df[bad_pos, ], 10))
}

# write CSV (exactly 3 columns, no rownames)
write.csv(map_df, "map_data.csv", row.names = FALSE, quote = FALSE)

# read your genetic map 
map_data <- read.csv("map_data.csv")
str(map_data)
head(map_data, 5)


# import markers effects for same SNP markers (additive and dominance) 
# They were calculated using GBLUP package using another code. Link: https://github.com/EhtishamSK/Add-Dom 

marker_effects <- read.csv("marker_effects.csv", row.names = 1)
# remove the markers ID columns here to match the requirements of Simple Mating package
marker_effects <- marker_effects[ , -1]
head(marker_effects, 10)


###################### Run Simple Mateing Models ###############################

# Planning crosses. Argument MateDesgin offer 4 options 

# 1. "full" = parents should be crossed with all parents, considering reciprocal crosses
# 2. "full_p" = all parents should be crossed with all parents, considering self of parents, and reciprocal crosses
# 3. "half" = neither reciprocal nor parents are included
# 4. "half_p" = all parents could be crossed with all parents, considering self of parents but not reciprocal crosses

crosPlan <- planCross(TargetPop = pheno$Taxa, MateDesign = "half")

# total crosses = 2145 see first 10 crosses 
head(crosPlan, 10)

# Thin population based on relatedness, if your markers are related. 
?relateThinning

# Creating criterion data frame
Crit <- data.frame(Id = pheno$Taxa, Criterion = pheno$AUDPC)
head(Crit, 10)
# Creating relationship matrix
ScaleMarkers <- scale(geno[, -1], scale = FALSE) # Exclude Taxa column
ScaleMarkers[1:5, 1:5]

relMat <- (ScaleMarkers %*% t(ScaleMarkers)) / ncol(ScaleMarkers)
relMat[1:5, 1:5]
# Thinning
parents2keep <- relateThinning(K = relMat, 
                               Criterion = Crit, 
                               threshold = 0.08, 
                               max.per.cluster = 5)


# None of the genotypes are thinned. Double check if values are small for relMat
summary(relMat)
hist(relMat)

# Estimate performance of progenies based on the parents' performance for single trait 
?getMPV

#### some checks to get before using MPV function ######
# Combine all parents in the MatePlan
all_parents <- unique(c(crosPlan$Parent1, crosPlan$Parent2))
# See if any parents are missing in Criterion
missing_ids <- setdiff(all_parents, Crit$Id)
missing_ids
# Check if rownames exist
rownames(relMat)
colnames(relMat)
# If NULL, assign them from Crit
rownames(relMat) <- Crit$Id
colnames(relMat) <- Crit$Id
#Reorder Crit to match relMat
Crit <- Crit[match(rownames(relMat), Crit$Id), ]

############## Run model for Single trait (ST) mid-parental value (mpv) #################
ST_mpv <- getMPV(MatePlan = crosPlan,
                 Criterion = Crit,
                 K = relMat)
head(ST_mpv, 20)
tail(ST_mpv, 20)

# Save output from getMPV
write.csv(ST_mpv, file = "SingleTrait_MPV.csv", row.names = FALSE)


############# calculate total genetic variance ###########################
?getTGV
# Though it was done in the start of the code. Do it again to make sure, there are no errors 
# This step is necessary for getTGV to run without errors, as marker_effects must have marker IDs as row names
# matching the marker columns in geno, with additive and dominance effects in columns 1 and 2, respectively.
marker_effects <- read.csv("marker_effects.csv", row.names = 1)
head(marker_effects, 10)

# Single trait Total genetic variance 
ST_tgv <- getTGV(MatePlan = crosPlan,
                 Markers = geno,
                 addEff = marker_effects[, 1],
                 domEff = marker_effects[, 2],
                 K = relMat, 
                 ploidy = 2)

head(ST_tgv, 10)

# Save output from getTGV
write.csv(ST_tgv, file = "SingleTrait_TGV.csv", row.names = FALSE)


############## Predicts usefulness component for a set of crosses #########################
?getUsefAD


# Check dimensions
dim(geno)              # should be nInd x nMarkers
dim(marker_effects)    # should be nMarkers x 2
dim(relMat)            # should be nInd x nInd
dim(map_data)          # should be nMarkers x 3 (or more)

# Check matching names
all(colnames(geno) == rownames(marker_effects))   # or check order
all(colnames(geno) %in% map_data$Marker)        # TRUE?
all(rownames(geno) %in% rownames(relMat))         # TRUE?

all(colnames(geno) == map_data$Marker)            # must be TRUE, not just %in%
all(rownames(geno) == rownames(relMat))


# now run the model 
usef_AD <- getUsefAD(MatePlan = crosPlan,
                     Markers = geno,
                     addEff = marker_effects$Add_Effect,
                     domEff = marker_effects$Dom_Effect,
                     K = relMat,
                     Map.In = map_data, # Chromosome and marker name
                     linkDes = NULL, # Linkage disequilibrium matrix
                     propSel = 0.05,
                     Method = "NonPhased")

# The model generates 2 outputs 

# 1. Cross.ID, parents (Parent 1 and Parent 2), progeny mean (Mean), progeny variance (Variance), 
# standard deviation (SD, progeny variance squared root), and usefulness for the trait/traits (Usefulness)
head(usef_AD[[1]], 10)

# 2. input for the optimization function, with the number of criteria used being the usefulness of the cross (Y).
head(usef_AD[[2]], 10)


# Save first output (cross summaries)
write.csv(usef_AD[[1]], file = "Cross_Summaries.csv", row.names = FALSE)

# Save second output (optimization input)
write.csv(usef_AD[[2]], file = "Optimization_Input.csv", row.names = FALSE)

##################### optimization function ####################

maxGain = selectCrosses(data = usef_AD[[2]],
                        n.cross = 200,
                        max.cross = 20,
                        min.cross = 1,
                        culling.pairwise.k = 0.5)

maxGain[[1]]
head(maxGain[[2]], 10)

# Save maxGain[[1]] to CSV (selected crosses)
write.csv(maxGain[[1]], file = "selected_crosses.csv", row.names = FALSE)

#Save maxGain[[2]] to CSV (data or metrics)
write.csv(maxGain[[2]], file = "crosses_data.csv", row.names = FALSE)

# plot 
maxGain[[3]]
