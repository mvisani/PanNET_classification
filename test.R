

# ChAMP Load With Preprocess Noob -----------------------------------------

# modified champ.load function

# ---- initial modification: 
# preprocessNoob and pattern matching for the sample sheet, which allows having multiple sheets in the same folder!
# --------------------------
# ---- 210421 addition: ---- 
# adding the percentage of failed probes as an attribute value to the data object. 
# This allows for easy extraction of this value for QC
# --------------------------
# ---- 220122 addition: ----
# In the original import function the sample sheet needs to be located in the data folder
# This is not very practical and makes it harder to separate raw and processed data 
# The function was altered to accept the sample sheet as a table and the raw data location as a variable
# This overrides the initial modifications made to the script
# --------------------------
# ---- 220615 addition: ----
# the preprocessing is extended to allow functional and ENmix normalization 
# A Function for plotting the EPIC bisulfite control probes is added
# --------------------------
# ---- 220727 addition: ----
# The loading function now accepts an rgSet to skip loading the idat files 
# This allows loading the idat files only one and then performing pre-processing and filtering on load
# The loading method needs to be set to minfi and a .Rds file of an RgSet needs to be provided
# The data needs to be rgSet_extended to allow filtering by p values and number of beads
# The plotting function for the bisulfite control probes is modified to make it easier to understand
# --------------------------

champ.load_extended <-function (directory = NULL, sampleSheet = NULL, method = "ChAMP", rgSet_object = NULL,
                                methValue = "B", autoimpute = TRUE, 
                                filterDetP = TRUE, ProbeCutoff = 0, SampleCutoff = 0.1, 
                                detPcut = 0.01, filterBeads = TRUE, beadCutoff = 0.05, filterNoCG = TRUE, 
                                filterSNPs = TRUE, population = NULL, filterMultiHit = TRUE, 
                                filterXY = TRUE, force = FALSE, arraytype = "450K",
                                preproc = "Raw", dyeMethod = "single", 
                                dataToInclude = c("loadQC", "mset", "mset_unfilt", "rgSet", "pd", "intensity", "beta", "detP")) 
{
  message("[===================================]")
  message("[<<<< ChAMP.LOAD_EXTENDED START >>>>]")
  message("-------------------------------------")
  message("function modified for output of QC parameters")
  mybeadcount <- function(x) {
    nb <- getNBeads(x)
    typeIadd <- getProbeInfo(x, type = "I")
    typeImatchA <- match(typeIadd$AddressA, rownames(nb))
    typeImatchB <- match(typeIadd$AddressB, rownames(nb))
    typeIIadd <- getProbeInfo(x, type = "II")
    typeIImatch <- match(typeIIadd$Address, rownames(nb))
    nbcg <- nb
    locusNames <- getManifestInfo(x, "locusNames")
    bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames), 
                      dimnames = list(locusNames, sampleNames(x)))
    TypeII.Name <- getProbeInfo(x, type = "II")$Name
    bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$AddressA, 
    ]
    TypeI <- getProbeInfo(x, type = "I")
    bcB <- bc_temp
    bcA <- bc_temp
    bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB, ]
    bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA, ]
    bcB3 <- which(bcB < 3)
    bcA3 <- which(bcA < 3)
    bcA2 <- bcA
    bcB2 <- bcB
    bcA2[bcA3] <- NA
    bcA2[bcB3] <- NA
    bc <- data.frame(bcA2)
    bc
  }
  myLoad <- list(mset = NA, mset_unfilt = NA, rgSet = NA, pd = NA, intensity = NA, 
                 beta = NA, detP = NA, loadQC = NA)
  if (method == "minfi") {
    if (is.null(rgSet_object)){
      directory = sampleSheet$Basename[1]
      message("\n[ Loading Data with Minfi Method ]")
      message("----------------------------------")
      message("Loading data from ", directory)
      # suppressWarnings(targets <- read.metharray.sheet(myDir)) # changed
      myLoad$rgSet <- read.metharray.exp(targets = sampleSheet, 
                                  extended = TRUE, 
                                  force = force)
    } else{
      message("using preloaded rgSet")
      myLoad$rgSet = rgSet_object
      # Only the extended rgSet contains all required information
      stopifnot(class(myLoad$rgSet) == "RGChannelSetExtended")
    }
    
    if (arraytype == "EPIC") 
      myLoad$rgSet@annotation <- c(array = "IlluminaHumanMethylationEPIC", 
                            annotation = "ilm10b4.hg19")
    sampleNames(myLoad$rgSet) = myLoad$rgSet[[1]]
    myLoad$pd <- pData(myLoad$rgSet)
    if (preproc == "Raw"){
      message("No preprocessing was done! This is the default setting for ChAMP")
      myLoad$mset <- preprocessRaw(myLoad$rgSet) 
    } else if (preproc == "Noob"){ # Noob background and dye bias correction
      message("Normexp Convolution model is used for background correction! \n")
      myLoad$mset <- preprocessNoob(myLoad$rgSet, dyeMethod = dyeMethod, verbose = TRUE)
    } else if (preproc == "Funnorm"){ # functional normalization (includes noob)
      message("Functional normalization is used for background correction \n")
      myLoad$mset = preprocessFunnorm(myLoad$rgSet, keepCN = F, ratioConvert = F) # the return type is a genomic methyl set 
    } else if (preproc == "ENmix"){
      message("EMmix background, dye bias and probe type correction")
      myLoad$mset = preprocessENmix(myLoad$rgSet) # according to the documentation the return type is again an RGset
      myLoad$mset = preprocessRaw(myLoad$mset) 
    }
    if(is.element("mset_unfilt", dataToInclude))
      myLoad$mset_unfilt = myLoad$mset
    else
      myLoad$mset_unfilt = "empty"
    myLoad$detP <- detectionP(myLoad$rgSet)
    message("<< Read DataSet Success. >>\n")
    if (methValue == "B") 
      tmp = getBeta(myLoad$mset, "Illumina")
    else tmp = getM(myLoad$mset)
    tmp[myLoad$detP >= detPcut] <- NA
    message("The fraction of failed positions per sample\n \n            (You may need to delete samples with high proportion of failed probes\n): ")
    numfail <- matrix(colMeans(is.na(tmp)))
    rownames(numfail) <- colnames(myLoad$detP)
    colnames(numfail) <- "Failed CpG Fraction."
    print(numfail)
    RemainSample <- which(numfail < SampleCutoff)
    if (any(numfail >= SampleCutoff)) 
      message("The detSamplecut parameter is : ", SampleCutoff, 
              "\nSamples : ", paste(rownames(numfail)[which(numfail >= 
                                                              SampleCutoff)], collapse = ","), " will be deleted.\n", 
              "There are ", length(RemainSample), " samples left for analysis.\n")
    myLoad$rgSet <- myLoad$rgSet[, RemainSample]
    myLoad$detP <- myLoad$detP[, RemainSample]
    myLoad$mset <- myLoad$mset[, RemainSample]
    myLoad$pd <- myLoad$pd[RemainSample, ]
    tmp <- tmp[, RemainSample]
    if (filterDetP) {
      mset.f = myLoad$mset[rowSums(is.na(tmp)) <= ProbeCutoff * 
                      ncol(myLoad$detP), ]
      
      # As an addition to the original function the filtered probes are saved in a variable to be able to save the information with the data 
      high.detPvalue <- dim(myLoad$mset)[1] - dim(mset.f)[1]
      
      if (ProbeCutoff == 0) {
        message("Filtering probes with a detection p-value above ", 
                detPcut, " in one or more samples has removed ", 
                high.detPvalue, " probes from the analysis. If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples.")
      }
      else {
        message("Filtering probes with a detection p-value above ", 
                detPcut, " in at least ", ProbeCutoff * 100, 
                "% of samples has removed ", high.detPvalue, " probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample file to identify potentially bad samples.")
      }
      myLoad$mset = mset.f
      tmp <- tmp[rowSums(is.na(tmp)) <= ProbeCutoff * 
                   ncol(myLoad$detP), ]
      message("<< Filter DetP Done. >>\n")
    }
    if (sum(is.na(tmp)) == 0) {
      message("\nThere is no NA values in your matrix, there is no need to imputation.\n")
    }
    else {
      message("\nThere are ", sum(is.na(tmp)), " NA remain in filtered Data Set. Impute can be done for remain NAs, but not suitable for small number samples. For small Data Set (like only 20 samples), we suggest you set parameter ProbeCutoff as 0 in champ.load() here, which would remove all NA involved probe no matter how many samples of those probes are NA.\n")
    }
    if (autoimpute & sum(is.na(tmp)) > 0) {
      message("Impute will be conducted here for remain ", 
              sum(is.na(tmp)), "  NAs. Note that if you don't do this, NA values would be kept in your data set. You may use champ.impute() function to do more complex imputation as well.")
      message("\nImpute function is working now, it may need couple minutes...")
      zz <- file("ImputeMessage.Rout", open = "wt")
      sink(zz)
      sink(zz, type = "message")
      tmp <- impute.knn(tmp, k = 5)$data
      sink(type = "message")
      sink()
      message("<< Imputation Done. >>\n")
    }
    if (filterBeads) {
      bc = mybeadcount(myLoad$rgSet)
      bc2 = bc[rowSums(is.na(bc)) < beadCutoff * (ncol(bc)), 
      ]
      mset.f2 = myLoad$mset[featureNames(myLoad$mset) %in% row.names(bc2), 
      ]
      
      # As an addition to the original function the filtered probes are saved in a variable to be able to save the information with the data 
      low.beadcount <- dim(myLoad$mset)[1] - dim(myLoad$mset.f2)[1]
      
      tmp <- tmp[rownames(tmp) %in% row.names(bc2), ]
      message("Filtering probes with a beadcount <3 in at least ", 
              beadCutoff * 100, "% of samples, has removed ", 
              low.beadcount, " from the analysis.")
      myLoad$mset = mset.f2
      message("<< Filter Beads Done. >>\n")
    }
    if (filterNoCG) {
      mset.f2 = dropMethylationLoci(myLoad$mset, dropCH = T)
      
      # As an addition to the original function the filtered probes are saved in a variable to be able to save the information with the data 
      no.CpG <- dim(myLoad$mset)[1] - dim(mset.f2)[1]
      
      tmp <- tmp[rownames(tmp) %in% featureNames(mset.f2), 
      ]
      message("Filtering non-cg probes, has removed ", 
              no.CpG, " from the analysis.")
      myLoad$mset <- mset.f2
      message("<< Filter NoCG Done. >>\n")
    }
    if (filterSNPs) {
      if (arraytype == "450K") {
        if (is.null(population)) {
          message("Using general 450K SNP list for filtering.")
          data(hm450.manifest.hg19)
          maskname <- rownames(hm450.manifest.hg19)[which(hm450.manifest.hg19$MASK_general == 
                                                            TRUE)]
        }
        else if (!population %in% c("AFR", "EAS", "EUR", 
                                    "SAS", "AMR", "GWD", "YRI", "TSI", "IBS", 
                                    "CHS", "PUR", "JPT", "GIH", "CHB", "STU", 
                                    "ITU", "LWK", "KHV", "FIN", "ESN", "CEU", 
                                    "PJL", "ACB", "CLM", "CDX", "GBR", "BEB", 
                                    "PEL", "MSL", "MXL", "ASW")) {
          message("Using general 450K SNP list for filtering.")
          data(hm450.manifest.hg19)
          maskname <- rownames(hm450.manifest.hg19)[which(hm450.manifest.hg19$MASK_general == 
                                                            TRUE)]
        }
        else {
          message("Using ", population, " specific 450K SNP list for filtering.")
          data(hm450.manifest.pop.hg19)
          maskname <- rownames(hm450.manifest.pop.hg19)[which(hm450.manifest.pop.hg19[, 
                                                                                      paste("MASK_general_", population, sep = "")] == 
                                                                TRUE)]
        }
      }
      else {
        if (is.null(population)) {
          message("Using general EPIC SNP list for filtering.")
          data(EPIC.manifest.hg19)
          maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general == 
                                                           TRUE)]
        }
        else if (!population %in% c("AFR", "EAS", "EUR", 
                                    "SAS", "AMR", "GWD", "YRI", "TSI", "IBS", 
                                    "CHS", "PUR", "JPT", "GIH", "CHB", "STU", 
                                    "ITU", "LWK", "KHV", "FIN", "ESN", "CEU", 
                                    "PJL", "ACB", "CLM", "CDX", "GBR", "BEB", 
                                    "PEL", "MSL", "MXL", "ASW")) {
          message("Using general EPIC SNP list for filtering.")
          data(EPIC.manifest.hg19)
          maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general == 
                                                           TRUE)]
        }
        else {
          message("Using ", population, " specific EPIC SNP list for filtering.")
          data(EPIC.manifest.pop.hg19)
          maskname <- rownames(EPIC.manifest.pop.hg19)[which(EPIC.manifest.pop.hg19[, 
                                                                                    paste("MASK_general_", population, sep = "")] == 
                                                               TRUE)]
        }
      }
      mset.f2 = myLoad$mset[!featureNames(myLoad$mset) %in% maskname, 
      ]
      
      # As an addition to the original function the filtered probes are saved in a variable to be able to save the information with the data 
      snp.probes <- dim(myLoad$mset)[1] - dim(mset.f2)[1]
      
      tmp <- tmp[!rownames(tmp) %in% maskname, ]
      message("Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper, 2016, has removed ", 
              snp.probes, " from the analysis.")
      myLoad$mset = mset.f2
      message("<< Filter SNP Done. >>\n")
    }
    if (filterMultiHit) {
      data(multi.hit)
      mset.f2 = myLoad$mset[!featureNames(myLoad$mset) %in% multi.hit$TargetID, 
      ]
      
      # As an addition to the original function the filtered probes are saved in a variable to be able to save the information with the data 
      multi.probes <- dim(myLoad$mset)[1] - dim(myLoad$mset.f2)[1]
      
      tmp <- tmp[!rownames(tmp) %in% multi.hit$TargetID, 
      ]
      message("Filtering probes that align to multiple locations as identified in Nordlund et al, has removed ", 
              multi.probes, " from the analysis.")
      myLoad$mset = mset.f2
      message("<< Filter MultiHit Done. >>\n")
    }
    if (filterXY) {
      if (arraytype == "EPIC") 
        data(probe.features.epic)
      else data(probe.features)
      autosomes = probe.features[!probe.features$CHR %in% 
                                   c("X", "Y"), ]
      mset.f2 = myLoad$mset[featureNames(myLoad$mset) %in% row.names(autosomes), 
      ]
      
      # As an addition to the original function the filtered probes are saved in a variable to be able to save the information with the data 
      XY.probes <- dim(myLoad$mset)[1] - dim(mset.f2)[1]
      
      tmp <- tmp[rownames(tmp) %in% row.names(autosomes), 
      ]
      message("Filtering probes on the X or Y chromosome has removed ", 
              XY.probes, " from the analysis.")
      myLoad$mset = mset.f2
      message("<< Filter XY chromosome Done. >>\n")
    }
    message(paste(if (methValue == "B") 
      "[Beta"
      else "[M", "value is selected as output.]\n"))
    myLoad$beta <- tmp
    myLoad$intensity <- minfi::getMeth(myLoad$mset) + minfi::getUnmeth(myLoad$mset)
    myLoad$detP <- myLoad$detP[which(row.names(myLoad$detP) %in% row.names(myLoad$beta)), 
    ]
    if (min(myLoad$beta, na.rm = TRUE) <= 0) 
      myLoad$beta[myLoad$beta <= 0] <- min(myLoad$beta[myLoad$beta > 
                                                0])
    message("Zeros in your dataset have been replaced with smallest positive value.\n")
    if (max(myLoad$beta, na.rm = TRUE) >= 0) 
      myLoad$beta[myLoad$beta >= 1] <- max(myLoad$beta[myLoad$beta < 
                                                1])
    message("One in your dataset have been replaced with largest value below 1.\n")
    message("The analysis will be proceed with ", dim(myLoad$beta)[1], 
            " probes and ", dim(myLoad$beta)[2], " samples.\n")
    message("Current Data Set contains ", sum(is.na(myLoad$beta)), 
            " NA in ", if (methValue == "B") 
              "[Beta]"
            else "[M]", " Matrix.\n")
    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
    
    message("adding QC information")
    # This is a new addition to the ChAMP load function
    # Originally I added is information via attibutes but it is much easier to simply append it to the output list
    myLoad$loadQC = list(fraction_failed_CpG = as.data.frame(numfail),
                  high_detection_pvalue = high.detPvalue,
                  low_bead_count = low.beadcount,
                  non_CpG_probes = no.CpG,
                  probes_containing_SNPs = snp.probes,
                  probes_mapping_multiple_locations = multi.probes,
                  probes_mapping_XY_chromosomes = XY.probes)
    
    
    return(myLoad[dataToInclude])
  }
  else {
    message("\n[ Loading Data with ChAMP Method ]")
    message("----------------------------------")
    message("Note that ChAMP method will NOT return rgSet or mset, they object defined by minfi. Which means, if you use ChAMP method to load data, you can not use SWAN or FunctionNormliazation method in champ.norm() (you can use BMIQ or PBC still). But All other function should not be influenced.\n")
    
    warning("Champ load is temporarily disabled because it does not readily work with the current data structure")
    
    #myImport <- champ.import(directory, arraytype = arraytype)
    #if (methValue == "B") 
    #  myLoad <- champ.filter(beta = myImport$beta, M = NULL, 
    #                         pd = myImport$pd, intensity = myImport$intensity, 
    #                         Meth = NULL, UnMeth = NULL, detP = myImport$detP, 
    #                         beadcount = myImport$beadcount, autoimpute = autoimpute, 
    #                         filterDetP = filterDetP, ProbeCutoff = ProbeCutoff, 
    #                         SampleCutoff = SampleCutoff, detPcut = detPcut, 
    #                         filterBeads = filterBeads, beadCutoff = beadCutoff, 
    #                         filterNoCG = filterNoCG, filterSNPs = filterSNPs, 
    #                         population = population, filterMultiHit = filterMultiHit, 
    #                         filterXY = filterXY, arraytype = arraytype)
    #else myLoad <- champ.filter(beta = NULL, M = myImport$M, 
    #                            pd = myImport$pd, intensity = myImport$intensity, 
    #                            Meth = NULL, UnMeth = NULL, detP = myImport$detP, 
    #                            beadcount = myImport$beadcount, autoimpute = autoimpute, 
    #                            filterDetP = filterDetP, ProbeCutoff = ProbeCutoff, 
    #                            SampleCutoff = SampleCutoff, detPcut = detPcut, 
    #                            filterBeads = filterBeads, beadCutoff = beadCutoff, 
    #                            filterNoCG = filterNoCG, filterSNPs = filterSNPs, 
    #                            population = population, filterMultiHit = filterMultiHit, 
    #                            filterXY = filterXY, arraytype = arraytype)
    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
    
    return(myLoad)
  }
}

# truncate decimals without rounding --------------------------------------

trunc_number_n_decimals <- function(numberToTrunc, nDecimals){
  numberToTrunc <- numberToTrunc + (10^-(nDecimals+5))
  splitNumber <- strsplit(x=format(numberToTrunc, digits=20, format=f), split="\\.")[[1]]
  decimalPartTrunc <- substr(x=splitNumber[2], start=1, stop=nDecimals)
  truncatedNumber <- as.numeric(paste0(splitNumber[1], ".", decimalPartTrunc))
  return(truncatedNumber)
}



# find top X most variables metabolites, genes, or similar ----------------

# requires: matrixstat

topXvariable <- function(df, topX = 1000){
  
  require(matrixStats)
  
  # calculate row var and sort
  Var <- rowVars(df)
  names(Var) <- rownames(df)
  Var <- sort(Var, decreasing = T)
  
  # find top X
  VarX <- MetaboliteVar[c(1:topX)]
  plot(VarX, main = paste("Top ", topX))
  df_out <- df[rownames(df) %in% names(VarX),]
}

# Find names of probes that exceed a chosen delta-beta threshold ----------

# this function has to be used to filter M-value based on the delta-beta cutoff,
# since this is not possible on a log-scale

NamesDeltaBetaCut <- function(pd, Beta, group1_Name, group2_Name, deltaBetaThreshold = 0.2){
  
  # Mark delta-beta > |0.2| probes for filtration after conversion to M-values
  
  Group1 <- paste("^",pd[which(pd$Sample_Group == group1_Name), "Sample_Name"],"$",
                  collapse = "|", sep = "")
  Group2 <- paste("^",pd[which(pd$Sample_Group == group2_Name), "Sample_Name"],"$",
                  collapse = "|", sep = "")
  cat(group1_Name, ":", Group1, "\n")
  cat(group2_Name, ":", Group2, "\n")
  
  Mean_Group1 <- rowMeans(Beta[,grep(Group1,colnames(Beta))])
  Mean_Group2 <- rowMeans(Beta[,grep(Group2,colnames(Beta))])
  
  
  # find name of the probes that are above delta-Beta threshold (absolute Value) 
  cut <- names(which((abs(Mean_Group1 - Mean_Group2) > deltaBetaThreshold) == T))
  
  # db <- abs(Mean_Group1 - Mean_Group2)
  # cut <- (db > deltaBetaThreshold)
  # cut2 <- which(cut == T)
  # cut3 <- names(cut2)
  # 
  
  return(cut) 
}


# Plot Beta- or/and M-values ---------------------------------------------

plot_BetaOrM <- function(data1, data2, pd, title = ""){
  
  require(RColorBrewer)
  par(mfrow=c(1,2))
  
  
  densityPlot(data1, sampGroups=pd$Sample_Group, main=deparse(substitute(data1)), 
              legend=FALSE)
  legend("top", legend = levels(factor(pd$Sample_Group)), 
         text.col=brewer.pal(8,"Pastel1"))
  densityPlot(data2, sampGroups=pd$Sample_Group, main=deparse(substitute(data2)), 
              legend=FALSE)
  legend("topleft", legend = levels(factor(pd$Sample_Group)), 
         text.col=brewer.pal(8,"Pastel1"))
  
  mtext(title, side = 3, line = -1, outer = TRUE)
  
}


# PLotting detection p value ----------------------------------------------

# (adpated from https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html)

plot_detP <- function(RGset, pd, detP, threshold){
  
  pd$ID <- paste(pd$Sample_Group,pd$Sample_Name,sep=".")
  
  sampleNames(RGset) <- pd$ID
  
  nGroup <- length(unique(pd$Sample_Group))
  
  if (nGroup < 3){
    nGroup <- 3
  }
  
  # examine mean detection p-values across all samples to identify any failed samples
  pal <- brewer.pal(nGroup,"Pastel1")
  
  par(mfrow=c(1,2))
  text("Mean Detection P-Value")
  
  barplot(colMeans(detP), col=pal[factor(pd$Sample_Group)], las=2, 
          cex.names=0.8, ylab="Mean detection p-values")
  abline(h=threshold,col="red", lty = 3)
  
  legend("topleft", legend=levels(factor(pd$Sample_Group)), fill=pal,
         bg="white")
  
  barplot(colMeans(detP), col=pal[factor(pd$Sample_Group)], las=2, 
          cex.names=0.8, ylim=c(0,(threshold+threshold)), ylab="Mean detection p-values")
  abline(h=threshold,col="red", lty = 3)
  
  legend("topleft", legend=levels(factor(pd$Sample_Group)), fill=pal, 
         bg="white")
  
}




# Convert Beta- to M-values -----------------------------------------------

Beta_To_M <- function(beta_matrix){
  return(log2(beta_matrix/(1-beta_matrix)))
}

# This function is closely modeled along the minfi function controlStripPlot
# For type I and II probes different plots are generated to better capture the different layouts
# --- NOTE ---
# For now this function is only for EPIC arrays 450K arrays have other probe numberings
# ------------
plotBSCProbes = function(rgSet, 
                         control_type = c("BISULFITE CONVERSION I", 
                                          "BISULFITE CONVERSION II")){
  control_type = match.arg(control_type)
  control_probes = getProbeInfo(rgSet, type = "Control") %>% 
    as.data.frame() %>% 
    filter(Type == control_type)
  # Fixing color names
  control_probes$Color[control_probes$Color == "Lime"] = "limegreen"
  
  control_red = getRed(rgSet[control_probes$Address, , drop = F]) %>% 
    as.data.frame() %>% 
    rownames_to_column("Address") %>% 
    mutate(detection_channel = "red")
  control_green = getGreen(rgSet[control_probes$Address, , drop = F]) %>% 
    as.data.frame() %>% 
    rownames_to_column("Address") %>% 
    mutate(detection_channel = "green")
  
  plotData = bind_rows(control_red,
                       control_green) %>% 
    pivot_longer(names_to = "sample", 
                 values_to = "intensity", 
                 - c(Address, detection_channel)) %>% 
    mutate(sample = factor(sample, 
                           levels = colnames(rgSet))) %>% 
    left_join(.,
              control_probes, 
              by = "Address")
  
  # The probe layout is different for type I and II probes 
  # Therefore the data is processed in two different ways
  if (control_type == "BISULFITE CONVERSION I"){
    plotData = plotData %>% 
      mutate(probe_channel = ifelse(grepl("[CU][12]$", ExtendedType),
                                    "green", 
                                    "red"))
    
    plot_green = plotData %>% 
      filter(detection_channel == "green" & probe_channel == "green") %>% 
      mutate(log_intensity = log2(intensity)) %>% 
      ggplot(aes(x = log_intensity, 
                 y = sample,
                 color = ExtendedType)) + 
      geom_point() +
      scale_color_manual(values = setNames(control_probes$Color,
                                           nm = control_probes$ExtendedType)) +
      coord_cartesian(xlim = c(6, 16)) +
      theme_classic() + 
      theme(panel.grid.major.y = element_line(color = "grey75"),
            axis.title.y = element_blank()) + 
      labs(x = "log2 ( intensity )",
           title = paste0("Bisulfite conversion type I ",
                          "green channel\n",
                          "high: C1, C2; low: U1, U2")) 
    
    plot_red = plotData %>% 
      filter(detection_channel == "red" & probe_channel == "red") %>% 
      mutate(log_intensity = log2(intensity)) %>%
      ggplot(aes(x = log_intensity, 
                 y = sample,
                 color = ExtendedType)) + 
      geom_point() +
      scale_color_manual(values = setNames(control_probes$Color,
                                           nm = control_probes$ExtendedType)) +
      coord_cartesian(xlim = c(6, 16)) +
      theme_classic() + 
      theme(panel.grid.major.y = element_line(color = "grey75"),
            axis.title.y = element_blank()) + 
      labs(x = "log2 ( intensity )",
           title = paste0("Bisulfite conversion type I ",
                          "red channel\n",
                          "high: C3, C4, C5; low: U3, U4, U5")) + 
      guides(size = "none")
    
    BSC_ratios = plotData %>% 
      filter(detection_channel == probe_channel) %>% 
      mutate(ExtendedType = str_extract(ExtendedType, 
                                        pattern = "[CU][0-9]")) %>% 
      separate(ExtendedType, 
               c(NA, "probe_type", "probe_number"), 
               sep = "") %>% 
      dplyr::select(sample, 
                    probe_channel,
                    probe_type, 
                    probe_number, 
                    intensity) %>% 
      pivot_wider(names_from = probe_type, 
                  values_from = intensity) %>% 
      rowwise() %>% 
      mutate(fraction_converted = C / sum(C, U))
    
    return(list(plot_green = plot_green,
                plot_red = plot_red,
                raw_intensities = plotData,
                BSC_ratios = BSC_ratios))
  } else{
    
    plot_green_red = 
      plotData %>% 
      mutate(log_intensity = log2(intensity)) %>% 
      ggplot(aes(x = log_intensity, 
                 y = sample, 
                 color = detection_channel)) + 
      geom_point(show.legend = F) + 
      scale_color_manual(values = setNames(
        RColorBrewer::brewer.pal(n = 3, name = "Set1")[c(1,3)],
        nm = c("red", "green"))) +
      coord_cartesian(xlim = c(6, 16)) + 
      theme_classic() + 
      theme(panel.grid.major.y = element_line(color = "grey75"),
            axis.title.y = element_blank()) + 
      labs(x = "log2 ( intensity )") 
    
    BSC_ratios = plotData %>% 
      separate(ExtendedType, 
               c(NA, "probe_number"), 
               sep = "-") %>% 
      dplyr::select(sample, 
                    detection_channel,
                    intensity,
                    probe_number) %>% 
      pivot_wider(names_from = detection_channel, 
                  values_from = intensity) %>% 
      rowwise() %>% 
      mutate(fraction_converted = red / sum(red, green)) 
    
    return(# In addition to the plot the plotData is saved 
      list(plot_green_red = plot_green_red,
           raw_intensities = plotData,
           BSC_ratios = BSC_ratios))
  }
}
