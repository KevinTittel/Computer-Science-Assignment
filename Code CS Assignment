# Computer Science for Business Analytics
# Individual Assignment
# Kevin Tittel (481044)

# Import packages and data -----------------------------------------------------
  install.packages("rjson") # To import the JSON data set
  library(rjson)
  install.packages("stringr") # To replace strings in the data set
  library(stringr)
  install.packages("stringi") # 
  library(stringi)
  install.packages("rlist")
  library(rlist)
  install.packages("purrr") # For flattening lists
  library(purrr)
  install.packages("base")
  install.packages("cluster") # To perform hierarchical clustering
  library(cluster)
  
  data <- unlist(fromJSON(file = "~/Downloads/TVs-all-merged.json"), recursive = FALSE)
  for (i in 1:length(data)) {
    data[[i]][["title"]] <- tolower(data[[i]][["title"]]) # Change to small letters
    data[[i]][["features"]] <- tolower(paste(stri_join_list(data[[i]][["featuresMap"]], 
                                      sep = "", collapse = NULL), collapse=" "))
    data[[i]][["title"]] <- str_replace_all(data[[i]][["title"]], 
                c("\"" = "inch", "inches" = "inch", "-inches" = "inch", 
                "-inch" = "inch", "-in" = "inch", "\\sin" = "inch", 
                " in$" = "inch", " inches" = "inch")) 
    data[[i]][["features"]] <- str_replace_all(data[[i]][["features"]], 
                c("\"" = "inch", "inches" = "inch", "-inches" = "inch", 
                "-inch" = "inch", "-in" = "inch", "\\sin" = "inch", 
                " in$" = "inch", " inches" = "inch")) 
    data[[i]][["title"]] <- str_replace_all(data[[i]][["title"]], 
                c("hertz" = "hz", "-hertz" = "hz", "-hz" = "hz", 
                  " hz$" = "hz", " hertz" = "hz", " hz" = "hz")) 
    data[[i]][["features"]] <- str_replace_all(data[[i]][["features"]], 
                c("hertz" = "hz", "-hertz" = "hz", "-hz" = "hz", 
                " hz$" = "hz", " hertz" = "hz", " hz" = "hz")) # Set equal words with same meaning
  }

# Set parameters ---------------------------------------------------------------
  numTVs <- length(data) # Number of TVs in the data set
  set.seed(123)

# Part 1: Create characteristic matrix -----------------------------------------
  # Create the unique and concise set of model words (i.e. characteristics) of the TVs 
    numTVs <- length(data) # Number of TVs in the data set
    modelWords <- vector(mode = "list", length = numTVs)
    for (i in 1:numTVs){
      modelWords[i] <- strsplit(data[[i]][["title"]], " ")
    }
    modelWords <- unlist(modelWords)
    modelWords <- unique(modelWords)
    modelWords <- modelWords[(regexpr("[a-zA-Z0-9]*(([0-9]+[^0-9, ]+)|([^0-9, ]+[0-9]+))[a-zA-Z0-9]*", modelWords) > 0)]
    modelWords <- c(modelWords, c("philips", "sharp", "aquos", "samsung", "nec", "tcl", "toshiba", "hisense",
                                  "sony", "bravia", "lg", "sanyo", "coby", "panasonic", "rca"))
    modelWords <- c(modelWords, c("led", "hdtv", "supersonic", "smart", "plasma", ))
    numModelWords <- length(modelWords)
    
  # Create a binary matrix, which appoints characteristics to the products 
    binaryMatrix <- matrix(data = 0, nrow = numModelWords, ncol = numTVs)
    for (i in 1:numModelWords) {
      for (j in 1:numTVs) {
        if (is.element(modelWords[i], unlist(strsplit(data[[j]][["title"]], " ")))) {
        # if (is.element(modelWords[i], unlist(strsplit(paste(data[[j]][["title"]], data[[j]][["features"]], collapse = ""), " ")))) {
          binaryMatrix[i,j] <- 1
        }
      }
    }

# Part 2: Signature Matrix using Min-Hashing -----------------------------------
  # Set parameters -------------------------------------------------------------
    numHashFunc <- 660 # Number of hash functions used for min-hashing
    
  # Compute the Min-Hash matrix
    rows <- 1:length(modelWords)
    minHashMatrix <- matrix(data = 0, nrow = numModelWords, ncol = numHashFunc)
    for (i in 1:numHashFunc) {
      a <- sample(1:numModelWords, 1)
      b <- sample(1:numModelWords, 1)
      minHashMatrix[,i] <- (a*rows + b)%%(numModelWords)
    }

  # Compute Signature Matrix
    signatureMatrix <- matrix(data = 1e+5, nrow = numHashFunc, ncol = numTVs)
    for (i in 1:numModelWords) {
      for (j in 1:numTVs) {
        if (binaryMatrix[i,j] == 1) {
          signatureMatrix[,j] <- pmin(signatureMatrix[,j], minHashMatrix[i,])
        }
      }
    }

# Part 3: LSH ------------------------------------------------------------------
  # Set necessary parameters and initialize matrices
    resultsMatrixLSH <- matrix(data = 0, nrow = 8, ncol = 11)
    rownames(resultsMatrixLSH) <- c("r", "b", "t", "ratio", "PC", "PQ", "F")
    numberRows <- c(2, 4, 5, 6, 10, 12, 15, 20, 30, 33, 44)
    
    # The higher r, the higher t, where r dominates b 
    # The higher t, the less candidate pairs, less FP (false positives) but more FN (false negatives)
    # Try: (b/r) in {(60/10); (55/12); (44/15); (33/20); (22/30); (20/33); (15/44)}
    
    # Find the total number of real duplicates first
    numRealDuplicates <- 0
    for (i in 1:(numTVs-1)) {
      for (j in (i+1):numTVs) {
        if (data[[i]][["modelID"]] == data[[j]][["modelID"]]) {
          numRealDuplicates <- numRealDuplicates + 1
        }
      }
    }
    
    for (l in 1:length(numberRows)) {
      numRows <- numberRows[l]
      numBands <- numHashFunc/numRows
      buckets <- matrix(data = 0, nrow = numBands, ncol = numTVs)
      candidatePairMatrix <- matrix(data = 0, nrow = numTVs, ncol = numTVs) 
      
      # Initialize candidate pair matrix
        for (i in 1:numBands) {
          buckets[i,] <- colSums(signatureMatrix[(((i-1)*numRows)+1):(i*numRows),])
        }
        for (i in 1:(numTVs-1)) {
          for (j in (i+1):numTVs) {
            if (sum((buckets[,i] == buckets[,j])) > 0) {
              candidatePairMatrix[i,j] <- 1
              candidatePairMatrix[j,i] <- 1
            } else {
            }
          }
        }
      
# Part 4: Evaluation LSH -------------------------------------------------------
  # Initialize counters for evaluation measures  
    numCorrectPred <- 0
    numFalsePred <- 0
    
  # Compute correct predictions and False Positives (FP)
    indices <- as.matrix(which(candidatePairMatrix == 1, arr.ind = TRUE, useNames = FALSE))
    indices <- indices[(indices[,1] < indices[,2]),]
    for (i in 1:nrow(indices)) {
      if (data[[indices[i,1]]][["modelID"]] == data[[indices[i,2]]][["modelID"]]) {
        numCorrectPred <- numCorrectPred + 1
      } else {
        numFalsePred <- numFalsePred + 1
      }
    }
  
  # Evaluation measures
    ratio <- (nrow(indices))/(((numTVs)*(numTVs-1))/2)
    PQ <- numCorrectPred / nrow(indices)
    PC <- numCorrectPred / numRealDuplicates
    Fstar <- (PQ * PC * 2) / (PQ + PC)
    
    resultsMatrixLSH[1,l] <- numRows
    resultsMatrixLSH[2,l] <- numBands
    resultsMatrixLSH[3,l] <- (1/numBands)^(1/numRows)
    resultsMatrixLSH[4,l] <- ratio
    resultsMatrixLSH[5,l] <- PC
    resultsMatrixLSH[6,l] <- PQ
    resultsMatrixLSH[7,l] <- Fstar

    
# Part 5: Similarity -----------------------------------------------------------
  # Initialize and fill distance matrix using Jaccard Similarity
    distanceMatrix <- matrix(data = NA, nrow = numTVs, ncol = numTVs)
    for (i in 1:(numTVs-1)) {
      for (j in (i+1):numTVs) {
        if (candidatePairMatrix[i,j] == 0) {
          distanceMatrix[i,j] <- 1e+5
          distanceMatrix[j,i] <- 1e+5
        } else {
          intersection <- sum(((binaryMatrix[,i] == 1) & (binaryMatrix[,j] == 1)))
          union <- sum(binaryMatrix[,i]) + sum(binaryMatrix[,j]) - intersection
          distanceMatrix[i,j] <- 1 - (intersection/union) # Jaccard Similarity
          distanceMatrix[j,i] <- 1 - (intersection/union)
        }
      }
    }
    distanceMatrix <- as.dist(distanceMatrix)

# Part 6: Clustering -----------------------------------------------------------
  # Set parameters
    clusterThreshold <- 0.40 # Parameter used in cutree() function
    
  # Hierarchical clustering
    clusters <- hclust(distanceMatrix, method = "complete")
    clusters <- cutree(clusters, h = clusterThreshold)
    uniqueClusters <- unique(clusters)
    newClusters <- vector(mode = "list", length = length(uniqueClusters))
    for (i in 1:length(uniqueClusters)) {
      newClusters[[i]] <- which(clusters == uniqueClusters[i])
    }

# Part 7: Final Evaluation -----------------------------------------------------
  # Initialize counters for TP and FP
    TP <- 0
    FP <- 0
    
  # Compute false positives   
    for (i in 1:length(newClusters)) {
      if (length(newClusters[[i]]) > 1) {
        pairs <- combn(newClusters[[i]], m = 2, simplify = FALSE)
        for (j in 1:length(pairs)) {
          if (data[[pairs[[j]][1]]][["modelID"]] == data[[pairs[[j]][2]]][["modelID"]]) {
            TP <- TP + 1
          } else {
            FP <- FP + 1
          }
        }
      }
    }
  
  # Compute measures
    precision <- TP/(TP+FP)
    recall <- TP/numRealDuplicates
    Fmeasure <- (2 * precision * recall)/(precision + recall)
    
    resultsMatrixLSH[8,l] <- Fmeasure
}

