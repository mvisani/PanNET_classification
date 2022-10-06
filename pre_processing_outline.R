# interesting project:
# 030
# 043

source("ChAMP_functions_Lionel_Philipp_220727.R")


# build raw data 
METH.rgSet = read.metharray.exp(targets = sampleSheet, # column Basename containing file locations
                                extended = TRUE, 
                                force = force)

# Build sample sheet / meta data 

# load raw data 
METH.rgSet = readRDS(raw_data)
# Check order of samples (identified by Sample_UID)

# get detection p values
detectionP(METH.rgSet)
getNBeads(Meth.rgSet)

# NOOB 
METH = champ.load_extended(rgSet_object = METH.rgSet,
                           sampleSheet = METH.sampleSheet,
                           method = "minfi", 
                           filterDetP = F, # include low p value probes
                           filterBeads = F, # include probes detected in few beads
                           detPcut = 1e-2,
                           arraytype = "EPIC", # set accorgingly
                           preproc = "Noob", 
                           dyeMethod = "single", 
                           dataToInclude = c("loadQC", "mset"),
                           force = F)

# convert to beta values (meth / unmeth + meth)
METH.beta = getBeta(METH$mset, "Illumina")[METH.shared_probes, ]

# Combine 450K and EPIC beta matrices (intersection)

# put probes in ascending order
METH.norm = champ.norm(METH.beta, 
                       arraytype = "450K", 
                       cores = 1,
                       resultsDir = "results")



