library(MixSim)           # Generate overlap parameters
library(mclust)           # k-means
library(cluster)          # pam
library(fclust)           # fuzzy k-means
library(FPDclustering)    # PD Clustering
library(ContaminatedMixt) # CND clustering
library(MixGHD)           # ARI
library(ggplot2)          # Plot
library(dplyr)            # Plot
library(ggpubr)           # Plot

set.seed(1234)            # Random seed

##########################################
########## 1. Simulate Datasets ##########
##########################################

# 1.1 Data 1: use MixSim & simdataset to simulate data
# Define variables
Pover <- c(0.1, 0.4) # proportion of "inlier" data [overlap]
Pout <- c(0.05, 0.15)  # proportion of "outlier" data
n_datasets <- 20     # number of datasets to generate for each group

# Initialize list for final data
FinalData <- list()

# Loop through proportion of "inlier" data and "outlier" data
for (i in 1:2) {
  for (j in 1:2) {
    
    # Loop through number of datasets
    for (m in 1:n_datasets) {
      
      # Generate mixture model and simulate dataset
      repeat{
        Q <- MixSim(MaxOmega = Pover[i], K = 2, p = 2)  
        if (Q$fail == 0) break
      }
      
      # In case, randomness error happens, don't stop
      # Go back to the beginning of the loop
      # Until simulation succeeds
      retry_simulation = TRUE
      while (retry_simulation) {
        tryCatch({A <- simdataset(n = 300 * (1 - Pout[j]), Pi = Q$Pi, Mu = Q$Mu, S = Q$S,
                                  n.out = 300 *Pout[j], 
                                  int = c(-2, 3), alpha = 0.01); retry_simulation = FALSE},
                 warning = function(w) {}, error = function(e) {}, finally = {})
      }
      
      colors <- c("red", "green", "magenta")
      plot(A$X, xlab = "x1", ylab = "x2", type = "n")        # A = dataset
      for (k in 0:2){
        points(A$X[A$id == k, ], col = colors[k+1], pch = k, cex = 0.8)
        points(x = Q$Mu[k,1], y = Q$Mu[k,2], col = colors[k+1], bg = colors[k+1], 
               pch = 20+k, cex = 1.2)
      }
      
      # Assign groups for outliers
      OutlierPerGroup = round(300*Pout[j]/2, 0)
      FirstOutlierIndex = 300 * (1 - Pout[j]) + 1
      A$id[FirstOutlierIndex : (FirstOutlierIndex+OutlierPerGroup-1)] <- 1
      A$id[(FirstOutlierIndex + OutlierPerGroup) : 300] <- 2
      
      # Add mixture model and simulated dataset to final data list
      dataset_name <- paste0("Pover", i, "_Pout", j, "_dataset", m)
      print(dataset_name)
      FinalData[[dataset_name]] <- list("Q" = Q, "data" = A)
    }
  }
}

# 1.2 Simulate Data 2: Use rCN() to simulate data
N1 = 150
N2 = 150

P.mu = list("Lower Overlap" = matrix(c(0, 0, -3, 3), nrow=2, ncol=2),
            "Higher Overlap" = matrix(c(0, 0, -1, 3), nrow=2, ncol=2))
P.sigma = list(matrix(c(1, -0.5, -0.5, 1), nrow=2, ncol=2),
               matrix(c(1, 0.5, 0.5, 1), nrow=2, ncol=2))
P.eta = c(20, 30)
P.alpha = c(0.95, 0.85)  # 0.8, 0.9 
n_datasets <- 20

# Initialize list for final data
FinalData_rCN <- list()
for (i in 1:2) {
  Mu = P.mu[[i]]               # 1=Higher Overlap, 2=Lower Overlap
  for (j in 1:2) {
    Alpha = P.alpha[[j]]       # alpha_1=0.9, alpha_2=0.8
    Eta = P.eta[[j]]           # eta_1=20, eta_2=30
    for (m in 1:n_datasets) {  # Repeat 20 times for different datasets
      # Generating data described in (b)
      
      # Component 1, using mu[1], sigma[1]
      X1 = rCN(N1, mu = Mu[1, ], Sigma = P.sigma[[1]],
               alpha = Alpha, eta = Eta)
      # Component 2, using mu[2], sigma[2]
      X2 = rCN(N2, mu = Mu[2, ], Sigma = P.sigma[[2]],
               alpha = Alpha, eta = Eta)
      # Concatenate 2 components to make a full dataset
      A = list("X" = rbind(X1, X2), 
               "id" = c(rep(1, N1), rep(2, N2)))
      
      # Add mixture model and simulated dataset to final data list
      dataset_name <- paste0(names(P.mu)[i],
                             "_Alpha_", Alpha,
                             "_Eta_", Eta,
                             "_dataset_", m)
      print(dataset_name)
      FinalData_rCN[[dataset_name]] <- list("data" = A)
      
      # Plot data
      colors <- c("red", "green", "magenta")
      plot(A$X, xlab = "x1", ylab = "x2", type = "n")        # A = dataset
      for (k in 0:2){
        points(A$X[A$id == k, ], col = colors[k+1], pch = k, cex = 0.7)
        points(x = Mu[k,1], y = Mu[k,2], col = colors[k+1], bg = colors[k+1], 
               pch = 20+k, cex = 1.4)
      }
    }
  }
}

# 1.3 Concatenate 2 simulated datasets
combine_data <- c(FinalData, FinalData_rCN)

#######################################
########## 2. Initialization ##########
#######################################


# Task: Use k-means, k-medoids (PAM), Fuzzy k-means, PD Clustering, at random to get initialization labels

# Initialization using the combine_data
kmeans_cluster_labels <- rep(list(list()), 160)  # For each G/K, size: n*1=80*1
pam_cluster_labels <- rep(list(list()), 160)
fkm_cluster_labels <- rep(list(list()), 160)
pdc_cluster_labels <- rep(list(list()), 160)
random_labels <- rep(list(list()), 160)

for (i in 1:160) {
  for (j in 2:5) {
    print(paste0("Processing ", names(combine_data)[i], ", using G=", j))
    kmeans_result <- kmeans(combine_data[[i]]$data$X, centers = j, nstart = 10) # Perform k-means clustering
    kmeans_cluster_labels[[i]][[j]] <- kmeans_result$cluster                 # Extract cluster labels from k-means
    
    pam_result <- pam(x=combine_data[[i]]$data$X,k=j)                           # Perform PAM/k-medoids clustering
    pam_cluster_labels[[i]][[j]] <- pam_result$clustering                    # Extract cluster labels from PAM/k-medoids
    
    fkm_result <- FKM(X=combine_data[[i]]$data$X, k=j)
    fkm_cluster_labels[[i]][[j]] <- unname(fkm_result$clus[,1])
    
    # PDC contains random initialization, which may lead to divided-by-zero error
    # Use a loop to catch the error and retry until success.
    retry_pdc = TRUE
    pdc_result = NULL
    while (retry_pdc) {
      tryCatch({pdc_result = PDC(combine_data[[i]]$data$X, k=j); retry_pdc = FALSE},
               warning = function(w) {}, error = function(e) {}, finally = {})
    }
    pdc_cluster_labels[[i]][[j]] <- pdc_result$label
    
    random_labels[[i]][[j]]<-sample(1:j,300,replace=TRUE)
  }
}

###############################################
########## 3. Model-based Clustering ##########
###############################################

# Task: Perform CND clustering with 5 different initialization

# 3.1 Transform Initialization Matrix to meet CNmixt's parameter, 'start.z', requirement
# CNmixt's start.z needs a n*G matrix. 
# Thus, we need to do some transformation of initialization label matrix, 
# which is generated in the Initialization part.
# For example, 1 will be changed to be [1,0] or [0,1]
to_z_matrix <- function(labels, G) {           # 
  labels = as.matrix(labels)
  z <- matrix(0, nrow(labels), G)              # Create a 0 matrix to store z matrix
  for (i in 1:nrow(labels)) {                
    z[i,labels[i]]=1                           # labels[i] == ith observation's column index in matrix z
  }
  return(z)
}

# 3.2 Run CNmixt model & Store Output into five lists corresponding to five initialization algorism
# G=2:5 for each dataset
CND_kmeans_result <- rep(list(list()), 160)
CND_pam_result <- rep(list(list()), 160)
CND_fkm_result <- rep(list(list()), 160)
CND_pdc_result <- rep(list(list()), 160)
CND_random_result <- rep(list(list()), 160)

for (i in 1:160) {
  for (j in 2:5) {
    print(paste0("Processing ", names(FinalData)[i], ", using G=", j))
    dataset=combine_data[[i]]$data$X
    retry = TRUE
    while (retry) {
      tryCatch({   # random number related errors happened when looping CNmixt. The function of this code is: if code returns error, do not stop.
        CND_kmeans_result[[i]][[j]] <- CNmixt(dataset, model=c("VVV"), G=j, initialization="manual", start.z = to_z_matrix(kmeans_cluster_labels[[i]][[j]], j), verbose=FALSE)
        CND_pam_result[[i]][[j]] <- CNmixt(dataset, model=c("VVV"), G=j, initialization="manual", start.z = to_z_matrix(pam_cluster_labels[[i]][[j]], j), verbose=FALSE)
        CND_fkm_result[[i]][[j]] <- CNmixt(dataset, model=c("VVV"), G=j, initialization="manual", start.z = to_z_matrix(fkm_cluster_labels[[i]][[j]], j), verbose=FALSE)
        CND_pdc_result[[i]][[j]] <- CNmixt(dataset, model=c("VVV"), G=j, initialization="manual", start.z = to_z_matrix(pdc_cluster_labels[[i]][[j]], j), verbose=FALSE)
        CND_random_result[[i]][[j]] <- CNmixt(dataset, model=c("VVV"), G=j, initialization="manual", start.z = to_z_matrix(random_labels[[i]][[j]], j), verbose=FALSE)
        retry = FALSE},
        warning = function(w) {},
        error = function(e) {
          print(e)
          print("Retry CNmixt to bypass random error")
        },
        finally = {})
    }
  }
}

###################################################
########## 4. Data Extraction & Analysis ##########
###################################################

# Task: Get the ARI for the best model for each dataset & Store the result in five lists that correspond to five initialization algorisms
# 4.1 New Function to find the best model: return G/k
findBestModel <- function(CND_results) {
  best_model_all_groups = NA
  for (j in 2:5) {
    best_model = getBestModel(CND_results[[j]], criterion = "BIC")
    if (any(is.na(best_model_all_groups))) {  # For the first iteration, keep G/K=2 as the best model
      best_model_all_groups = best_model
    } else if (best_model$models[[1]]$IC$BIC > best_model_all_groups$models[[1]]$IC$BIC){
      best_model_all_groups = best_model
    }
  }
  return(best_model_all_groups)
}

# 4.2 Find ARI for each best model & each initialization algorism
best_ari_kmeans_result <- list()
best_ari_pam_result <- list()
best_ari_fkm_result <- list()
best_ari_pdc_result <- list()
best_ari_random_result <- list()

for (i in 1:160) {
  best_kmeans_model = findBestModel(CND_kmeans_result[[i]])
  best_ari_kmeans_result[[i]] = list(
    "ARI" = ARI(getCluster(best_kmeans_model),combine_data[[i]]$data$id),
    "G" = best_kmeans_model$models[[1]]$G,
    "Detection" = sum(best_kmeans_model$models[[1]]$detection[,2]=="bad")/300
  )
  
  best_pam_model = findBestModel(CND_pam_result[[i]])
  best_ari_pam_result[[i]] = list(
    "ARI" = ARI(getCluster(best_pam_model),combine_data[[i]]$data$id),
    "G" = best_pam_model$models[[1]]$G,
    "Detection" = sum(best_pam_model$models[[1]]$detection[,2]=="bad")/300
  )
  
  best_fkm_model = findBestModel(CND_fkm_result[[i]])
  best_ari_fkm_result[[i]] = list(
    "ARI" = ARI(getCluster(best_fkm_model),combine_data[[i]]$data$id),
    "G" = best_fkm_model$models[[1]]$G,
    "Detection" = sum(best_fkm_model$models[[1]]$detection[,2]=="bad")/300
  )
  
  best_pdc_model = findBestModel(CND_pdc_result[[i]])
  best_ari_pdc_result[[i]] = list(
    "ARI" = ARI(getCluster(best_pdc_model),combine_data[[i]]$data$id),
    "G" = best_pdc_model$models[[1]]$G,
    "Detection" = sum(best_pdc_model$models[[1]]$detection[,2]=="bad")/300
  )
  
  best_random_model = findBestModel(CND_random_result[[i]])
  best_ari_random_result[[i]] = list(
    "ARI" = ARI(getCluster(best_random_model),combine_data[[i]]$data$id),
    "G" = best_random_model$models[[1]]$G,
    "Detection" = sum(best_random_model$models[[1]]$detection[,2]=="bad")/300
  )
}

#############################
########## 5. Plot ##########
#############################

# Task: Plot best model's ARI, G & rCN data sets' detection

# 5.1 Extract ARI, G & Detection (rCN only) for the best model for each dataset & each criteria
# 5.1.1 ARI
C1_ARI <- matrix(NA, 20, 5)
C2_ARI <- matrix(NA, 20, 5)
C3_ARI <- matrix(NA, 20, 5)
C4_ARI <- matrix(NA, 20, 5)
C1_ARI_rCN <- matrix(NA, 20, 5)
C2_ARI_rCN <- matrix(NA, 20, 5)
C3_ARI_rCN <- matrix(NA, 20, 5)
C4_ARI_rCN <- matrix(NA, 20, 5)

for (i in 1:20) {
  C1_ARI[i,] <- c(best_ari_kmeans_result[[i]]$ARI, best_ari_pam_result[[i]]$ARI, best_ari_fkm_result[[i]]$ARI,
               best_ari_pdc_result[[i]]$ARI, best_ari_random_result[[i]]$ARI)
  C2_ARI[i,] <- c(best_ari_kmeans_result[[i+20]]$ARI, best_ari_pam_result[[i+20]]$ARI, best_ari_fkm_result[[i+20]]$ARI,
              best_ari_pdc_result[[i+20]]$ARI, best_ari_random_result[[i+20]]$ARI)
  C3_ARI[i,] <- c(best_ari_kmeans_result[[i+40]]$ARI, best_ari_pam_result[[i+40]]$ARI, best_ari_fkm_result[[i+40]]$ARI,
              best_ari_pdc_result[[i+40]]$ARI, best_ari_random_result[[i+40]]$ARI)
  C4_ARI[i,] <- c(best_ari_kmeans_result[[i+60]]$ARI, best_ari_pam_result[[i+60]]$ARI, best_ari_fkm_result[[i+60]]$ARI,
              best_ari_pdc_result[[i+60]]$ARI, best_ari_random_result[[i+60]]$ARI)
  C1_ARI_rCN[i,] <- c(best_ari_kmeans_result[[i+80]]$ARI, best_ari_pam_result[[i+80]]$ARI, best_ari_fkm_result[[i+80]]$ARI,
                  best_ari_pdc_result[[i+80]]$ARI, best_ari_random_result[[i+80]]$ARI)
  C2_ARI_rCN[i,] <- c(best_ari_kmeans_result[[i+100]]$ARI, best_ari_pam_result[[i+100]]$ARI, best_ari_fkm_result[[i+100]]$ARI,
                      best_ari_pdc_result[[i+100]]$ARI, best_ari_random_result[[i+100]]$ARI)
  C3_ARI_rCN[i,] <- c(best_ari_kmeans_result[[i+120]]$ARI, best_ari_pam_result[[i+120]]$ARI, best_ari_fkm_result[[i+120]]$ARI,
                      best_ari_pdc_result[[i+120]]$ARI, best_ari_random_result[[i+120]]$ARI)
  C4_ARI_rCN[i,] <- c(best_ari_kmeans_result[[i+140]]$ARI, best_ari_pam_result[[i+140]]$ARI, best_ari_fkm_result[[i+140]]$ARI,
                      best_ari_pdc_result[[i+140]]$ARI, best_ari_random_result[[i+140]]$ARI)
}

# 5.1.2 G
C1_G <- matrix(NA, 20, 5)
C2_G <- matrix(NA, 20, 5)
C3_G <- matrix(NA, 20, 5)
C4_G <- matrix(NA, 20, 5)
C1_G_rCN <- matrix(NA, 20, 5)
C2_G_rCN <- matrix(NA, 20, 5)
C3_G_rCN <- matrix(NA, 20, 5)
C4_G_rCN <- matrix(NA, 20, 5)

for (i in 1:20) {
  C1_G[i,] <- c(best_ari_kmeans_result[[i]]$G, best_ari_pam_result[[i]]$G, best_ari_fkm_result[[i]]$G,
                  best_ari_pdc_result[[i]]$G, best_ari_random_result[[i]]$G)
  C2_G[i,] <- c(best_ari_kmeans_result[[i+20]]$G, best_ari_pam_result[[i+20]]$G, best_ari_fkm_result[[i+20]]$G,
                  best_ari_pdc_result[[i+20]]$G, best_ari_random_result[[i+20]]$G)
  C3_G[i,] <- c(best_ari_kmeans_result[[i+40]]$G, best_ari_pam_result[[i+40]]$G, best_ari_fkm_result[[i+40]]$G,
                  best_ari_pdc_result[[i+40]]$G, best_ari_random_result[[i+40]]$G)
  C4_G[i,] <- c(best_ari_kmeans_result[[i+60]]$G, best_ari_pam_result[[i+60]]$G, best_ari_fkm_result[[i+60]]$G,
                  best_ari_pdc_result[[i+60]]$G, best_ari_random_result[[i+60]]$G)
  C1_G_rCN[i,] <- c(best_ari_kmeans_result[[i+80]]$G, best_ari_pam_result[[i+80]]$G, best_ari_fkm_result[[i+80]]$G,
                best_ari_pdc_result[[i+80]]$G, best_ari_random_result[[i+80]]$G)
  C2_G_rCN[i,] <- c(best_ari_kmeans_result[[i+100]]$G, best_ari_pam_result[[i+100]]$G, best_ari_fkm_result[[i+100]]$G,
                    best_ari_pdc_result[[i+100]]$G, best_ari_random_result[[i+100]]$G)
  C3_G_rCN[i,] <- c(best_ari_kmeans_result[[i+120]]$G, best_ari_pam_result[[i+120]]$G, best_ari_fkm_result[[i+120]]$G,
                    best_ari_pdc_result[[i+120]]$G, best_ari_random_result[[i+120]]$G)
  C4_G_rCN[i,] <- c(best_ari_kmeans_result[[i+140]]$G, best_ari_pam_result[[i+140]]$G, best_ari_fkm_result[[i+140]]$G,
                    best_ari_pdc_result[[i+140]]$G, best_ari_random_result[[i+140]]$G)
}

# 5.1.3 Detection for rCN data
C1_Detection <- matrix(NA, 20, 5)
C2_Detection <- matrix(NA, 20, 5)
C3_Detection <- matrix(NA, 20, 5)
C4_Detection <- matrix(NA, 20, 5)
C1_Detection_rCN <- matrix(NA, 20, 5)
C2_Detection_rCN <- matrix(NA, 20, 5)
C3_Detection_rCN <- matrix(NA, 20, 5)
C4_Detection_rCN <- matrix(NA, 20, 5)

for (i in 1:20) {
  C1_Detection[i,] <- c(best_ari_kmeans_result[[i]]$Detection, best_ari_pam_result[[i]]$Detection, best_ari_fkm_result[[i]]$Detection,
                            best_ari_pdc_result[[i]]$Detection, best_ari_random_result[[i]]$Detection)
  C2_Detection[i,] <- c(best_ari_kmeans_result[[i+20]]$Detection, best_ari_pam_result[[i+20]]$Detection, best_ari_fkm_result[[i+20]]$Detection,
                            best_ari_pdc_result[[i+20]]$Detection, best_ari_random_result[[i+20]]$Detection)
  C3_Detection[i,] <- c(best_ari_kmeans_result[[i+40]]$Detection, best_ari_pam_result[[i+40]]$Detection, best_ari_fkm_result[[i+40]]$Detection,
                            best_ari_pdc_result[[i+40]]$Detection, best_ari_random_result[[i+40]]$Detection)
  C4_Detection[i,] <- c(best_ari_kmeans_result[[i+60]]$Detection, best_ari_pam_result[[i+60]]$Detection, best_ari_fkm_result[[i+60]]$Detection,
                            best_ari_pdc_result[[i+60]]$Detection, best_ari_random_result[[i+60]]$Detection)
  C1_Detection_rCN[i,] <- c(best_ari_kmeans_result[[i+80]]$Detection, best_ari_pam_result[[i+80]]$Detection, best_ari_fkm_result[[i+80]]$Detection,
                      best_ari_pdc_result[[i+80]]$Detection, best_ari_random_result[[i+80]]$Detection)
  C2_Detection_rCN[i,] <- c(best_ari_kmeans_result[[i+100]]$Detection, best_ari_pam_result[[i+100]]$Detection, best_ari_fkm_result[[i+100]]$Detection,
                      best_ari_pdc_result[[i+100]]$Detection, best_ari_random_result[[i+100]]$Detection)
  C3_Detection_rCN[i,] <- c(best_ari_kmeans_result[[i+120]]$Detection, best_ari_pam_result[[i+120]]$Detection, best_ari_fkm_result[[i+120]]$Detection,
                      best_ari_pdc_result[[i+120]]$Detection, best_ari_random_result[[i+120]]$Detection)
  C4_Detection_rCN[i,] <- c(best_ari_kmeans_result[[i+140]]$Detection, best_ari_pam_result[[i+140]]$Detection, best_ari_fkm_result[[i+140]]$Detection,
                      best_ari_pdc_result[[i+140]]$Detection, best_ari_random_result[[i+140]]$Detection)
}

# 5.2 Plot ARI
# 5.2.1 FinalData ARI
C1_dataframe = data.frame(
  initial_method=c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(C1_ARI),
  G = c(C1_G),
  Detection = c(C1_Detection),
  Detected = c(C1_Detection)>0)

C2_dataframe = data.frame(
  initial_method=c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(C2_ARI),
  G = c(C2_G),
  Detection = c(C2_Detection),
  Detected = c(C2_Detection)>0)

C3_dataframe = data.frame(
  initial_method=c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(C3_ARI),
  G = c(C3_G),
  Detection = c(C3_Detection),
  Detected = c(C3_Detection)>0)

C4_dataframe = data.frame(
  initial_method=c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(C4_ARI),
  G = c(C4_G),
  Detection = c(C4_Detection),
  Detected = c(C4_Detection)>0)

plot_C1 <-ggplot(C1_dataframe, aes(x=initial_method, y=ARI, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +
  xlab("Initialization Algorisms") +
  ylab("ARI") +
  ggtitle("Scenario 1: Overlap 10% & Outlier 5% | n=300") +
  coord_cartesian(ylim = c(0, 1))


plot_C2 <-ggplot(C2_dataframe, aes(x=initial_method, y=ARI, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +
  xlab("Initialization Algorisms") +
  ylab("ARI") +
  ggtitle("Scenario 2: Overlap 10% & Outlier 15% | n=300") +
  coord_cartesian(ylim = c(0, 1)) 

plot_C3 <-ggplot(C3_dataframe, aes(x=initial_method, y=ARI, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +
  xlab("Initialization Algorisms") +
  ylab("ARI") +
  ggtitle("Scenario 3: Overlap 40% & Outlier 5% | n=300") +
  coord_cartesian(ylim = c(0, 1)) 

plot_C4 <-ggplot(C4_dataframe, aes(x=initial_method, y=ARI, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +
  xlab("Initialization Algorisms") +
  ylab("ARI") +
  ggtitle("Scenario 4: Overlap 40% & Outlier 15% | n=300") +
  coord_cartesian(ylim = c(0, 1))

ggarrange(plot_C1, plot_C2, plot_C3, plot_C4, ncol=2, nrow = 2)

mean(C4_dataframe[1:20,2])   # k-means
mean(C4_dataframe[21:40,2])  # PAM
mean(C4_dataframe[41:60,2])  # FKM
mean(C4_dataframe[61:80,2])  # PDC
mean(C4_dataframe[81:100,2]) # Random

# 5.2.2 FinalData_cRN ARI
C1_dataframe_rCN = data.frame(
  initial_method=c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(C1_ARI_rCN),
  G = c(C1_G_rCN),
  Detection = c(C1_Detection_rCN),
  Detected = c(C1_Detection_rCN)>0)

C2_dataframe_rCN = data.frame(
  initial_method=c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(C2_ARI_rCN),
  G = c(C2_G_rCN),
  Detection = c(C2_Detection_rCN),
  Detected = c(C2_Detection_rCN)>0)

C3_dataframe_rCN = data.frame(
  initial_method=c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(C3_ARI_rCN),
  G = c(C3_G_rCN),
  Detection = c(C3_Detection_rCN),
  Detected = c(C3_Detection_rCN)>0)

C4_dataframe_rCN = data.frame(
  initial_method=c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(C4_ARI_rCN),
  G = c(C4_G_rCN),
  Detection = c(C4_Detection_rCN),
  Detected = c(C4_Detection_rCN)>0)

plot_C1_rCN <-ggplot(C1_dataframe_rCN, aes(x=initial_method, y=ARI, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("ARI") +
  ggtitle("Scenario 1: Low Overlap & alpha=0.95, eta=20 | n=300") +
  coord_cartesian(ylim = c(0, 1))

plot_C2_rCN <-ggplot(C2_dataframe_rCN, aes(x=initial_method, y=ARI, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("ARI") +
  ggtitle("Scenario 2: Low Overlap & alpha=0.85, eta=30 | n=300") +
  coord_cartesian(ylim = c(0, 1))

plot_C3_rCN <-ggplot(C3_dataframe_rCN, aes(x=initial_method, y=ARI, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("ARI") +
  ggtitle("Scenario 3: High Overlap & alpha=0.95, eta=20 | n=300") +
  coord_cartesian(ylim = c(0, 1))

plot_C4_rCN <-ggplot(C4_dataframe_rCN, aes(x=initial_method, y=ARI, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("ARI") +
  ggtitle("Scenario 4: High Overlap & alpha=0.85, eta=30 | n=300") +
  coord_cartesian(ylim = c(0, 1)) 

ggarrange(plot_C1_rCN, plot_C2_rCN, plot_C3_rCN, plot_C4_rCN, ncol=2, nrow = 2)

# 5.3 Plot G
# 5.3.1 FinalData
C1_group_dataframe = as.data.frame(table(C1_dataframe$initial_method, C1_dataframe$G))
C2_group_dataframe = as.data.frame(table(C2_dataframe$initial_method, C2_dataframe$G))
C3_group_dataframe = as.data.frame(table(C3_dataframe$initial_method, C3_dataframe$G))
C4_group_dataframe = as.data.frame(table(C4_dataframe$initial_method, C4_dataframe$G))

names(C1_group_dataframe)[1] <- "Initialization"
names(C1_group_dataframe)[2] <- "G"
names(C1_group_dataframe)[3] <- "Frequency"

names(C2_group_dataframe)[1] <- "Initialization"
names(C2_group_dataframe)[2] <- "G"
names(C2_group_dataframe)[3] <- "Frequency"

names(C3_group_dataframe)[1] <- "Initialization"
names(C3_group_dataframe)[2] <- "G"
names(C3_group_dataframe)[3] <- "Frequency"

names(C4_group_dataframe)[1] <- "Initialization"
names(C4_group_dataframe)[2] <- "G"
names(C4_group_dataframe)[3] <- "Frequency"

plot_C1_G <- ggplot(C1_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 1: Overlap 10% & Outlier 5% | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C2_G <- ggplot(C2_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 2: Overlap 10% & Outlier 15% | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C3_G <- ggplot(C3_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 3: Overlap 40% & Outlier 5% | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C4_G <- ggplot(C4_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 4: Overlap 40% & Outlier 15% | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

ggarrange(plot_C1_G, plot_C2_G, plot_C3_G, plot_C4_G, ncol=2, nrow = 2)

# 5.3.2 FinalData_rCN
C1_rCN_group_dataframe = as.data.frame(table(C1_dataframe_rCN$initial_method, C1_dataframe_rCN$G))
C2_rCN_group_dataframe = as.data.frame(table(C2_dataframe_rCN$initial_method, C2_dataframe_rCN$G))
C3_rCN_group_dataframe = as.data.frame(table(C3_dataframe_rCN$initial_method, C3_dataframe_rCN$G))
C4_rCN_group_dataframe = as.data.frame(table(C4_dataframe_rCN$initial_method, C4_dataframe_rCN$G))

names(C1_rCN_group_dataframe)[1] <- "Initialization"
names(C1_rCN_group_dataframe)[2] <- "G"
names(C1_rCN_group_dataframe)[3] <- "Frequency"

names(C2_rCN_group_dataframe)[1] <- "Initialization"
names(C2_rCN_group_dataframe)[2] <- "G"
names(C2_rCN_group_dataframe)[3] <- "Frequency"

names(C3_rCN_group_dataframe)[1] <- "Initialization"
names(C3_rCN_group_dataframe)[2] <- "G"
names(C3_rCN_group_dataframe)[3] <- "Frequency"

names(C4_rCN_group_dataframe)[1] <- "Initialization"
names(C4_rCN_group_dataframe)[2] <- "G"
names(C4_rCN_group_dataframe)[3] <- "Frequency"

plot_C1_G_rCN <- ggplot(C1_rCN_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 1: Low Overlap & alpha=0.95, eta=20 | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C2_G_rCN <- ggplot(C2_rCN_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 2: Low Overlap & alpha=0.85, eta=30 | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C3_G_rCN <- ggplot(C3_rCN_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 3: High Overlap & alpha=0.95, eta=20 | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C4_G_rCN <- ggplot(C4_rCN_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 4: High Overlap & alpha=0.85, eta=30 | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

ggarrange(plot_C1_G_rCN, plot_C2_G_rCN, plot_C3_G_rCN, plot_C4_G_rCN, ncol=2, nrow = 2)

# 5.3.3 Plot Ooutlier Detection
# 5.3.3.1 Barplot
# MixSim Data
plot_C1_Detection <-ggplot(C1_dataframe, aes(x=initial_method, y=Detection, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("Outlier Detection") +
  ggtitle("Scenario 1: Low Overlap & Low Outlier| n=300") +
  coord_cartesian(ylim = c(0, 0.2)) 

plot_C2_Detection <-ggplot(C2_dataframe, aes(x=initial_method, y=Detection, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) + xlab("Initialization Algorisms") +
  ylab("Outlier Detection") +
  ggtitle("Scenario 2: Low Overlap & High Outlier | n=300") +
  coord_cartesian(ylim = c(0, 0.2)) 

plot_C3_Detection <-ggplot(C3_dataframe, aes(x=initial_method, y=Detection, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("Outlier Detection") +
  ggtitle("Scenario 3: High Overlap & Low Outlier | n=300") +
  coord_cartesian(ylim = c(0, 0.2))

plot_C4_Detection <-ggplot(C4_dataframe, aes(x=initial_method, y=Detection, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("Outlier Detection") +
  ggtitle("Scenario 4: High Overlap & High Outlier | n=300") +
  coord_cartesian(ylim = c(0, 0.2)) 

ggarrange(plot_C1_Detection, plot_C2_Detection, plot_C3_Detection, plot_C4_Detection, ncol=2, nrow = 2)

# rCN data
plot_C1_rCN_Detection <-ggplot(C1_dataframe_rCN, aes(x=initial_method, y=Detection, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("Outlier Detection") +
  ggtitle("Scenario 1: Low Overlap & alpha=0.95, eta=20 | n=300") +
  coord_cartesian(ylim = c(0, 0.2)) 

plot_C2_rCN_Detection <-ggplot(C2_dataframe_rCN, aes(x=initial_method, y=Detection, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("Outlier Detection") +
  ggtitle("Scenario 2: Low Overlap & alpha=0.85, eta=30 | n=300") +
  coord_cartesian(ylim = c(0, 0.2))

plot_C3_rCN_Detection <-ggplot(C3_dataframe_rCN, aes(x=initial_method, y=Detection, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("Outlier Detection") +
  ggtitle("Scenario 3: High Overlap & alpha=0.95, eta=20 | n=300") +
  coord_cartesian(ylim = c(0, 0.2))

plot_C4_rCN_Detection <-ggplot(C4_dataframe_rCN, aes(x=initial_method, y=Detection, fill=initial_method)) + 
  stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(alpha=1) + theme_bw() +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5,color = "black", size = 18, face = "bold"),
        axis.text = element_text(color = "black",size = 12),
        axis.title=element_text(color = "black",size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=rainbow(5)) +  xlab("Initialization Algorisms") +
  ylab("Outlier Detection") +
  ggtitle("Scenario 4: High Overlap & alpha=0.85, eta=30 | n=300") +
  coord_cartesian(ylim = c(0, 0.2)) 

ggarrange(plot_C1_rCN_Detection, plot_C2_rCN_Detection, plot_C3_rCN_Detection, plot_C4_rCN_Detection, ncol=2, nrow = 2)

### 5.3.3.2 mean of percentage of detected outliers
### MixSim data sets
mean(C1_dataframe[1:20,4])   # k-means
mean(C1_dataframe[21:40,4])  # PAM
mean(C1_dataframe[41:60,4])  # FKM
mean(C1_dataframe[61:80,4])  # PDC
mean(C1_dataframe[81:100,4]) # Random

mean(C2_dataframe[1:20,4])   # k-means
mean(C2_dataframe[21:40,4])  # PAM
mean(C2_dataframe[41:60,4])  # FKM
mean(C2_dataframe[61:80,4])  # PDC
mean(C2_dataframe[81:100,4]) # Random

mean(C3_dataframe[1:20,4])   # k-means
mean(C3_dataframe[21:40,4])  # PAM
mean(C3_dataframe[41:60,4])  # FKM
mean(C3_dataframe[61:80,4])  # PDC
mean(C3_dataframe[81:100,4]) # Random

mean(C4_dataframe[1:20,4])   # k-means
mean(C4_dataframe[21:40,4])  # PAM
mean(C4_dataframe[41:60,4])  # FKM
mean(C4_dataframe[61:80,4])  # PDC
mean(C4_dataframe[81:100,4]) # Random

### rCN data
mean(C1_dataframe_rCN[1:20,4])   # k-means
mean(C1_dataframe_rCN[21:40,4])  # PAM
mean(C1_dataframe_rCN[41:60,4])  # FKM
mean(C1_dataframe_rCN[61:80,4])  # PDC
mean(C1_dataframe_rCN[81:100,4]) # Random

mean(C2_dataframe_rCN[1:20,4])   # k-means
mean(C2_dataframe_rCN[21:40,4])  # PAM
mean(C2_dataframe_rCN[41:60,4])  # FKM
mean(C2_dataframe_rCN[61:80,4])  # PDC
mean(C2_dataframe_rCN[81:100,4]) # Random

mean(C3_dataframe_rCN[1:20,4])   # k-means
mean(C3_dataframe_rCN[21:40,4])  # PAM
mean(C3_dataframe_rCN[41:60,4])  # FKM
mean(C3_dataframe_rCN[61:80,4])  # PDC
mean(C3_dataframe_rCN[81:100,4]) # Random


mean(C4_dataframe_rCN[1:20,4])   # k-means
mean(C4_dataframe_rCN[21:40,4])  # PAM
mean(C4_dataframe_rCN[41:60,4])  # FKM
mean(C4_dataframe_rCN[61:80,4])  # PDC
mean(C4_dataframe_rCN[81:100,4]) # Random

### 5.3.3.3 Frequency of outlier detection table 
# MixSim data
table(C1_dataframe$G, C1_dataframe$Detected)
table(C2_dataframe$G, C2_dataframe$Detected)
table(C3_dataframe$G, C3_dataframe$Detected)
table(C4_dataframe$G, C4_dataframe$Detected)

# rCN data
table(C1_dataframe_rCN$G, C1_dataframe_rCN$Detected)
table(C2_dataframe_rCN$G, C2_dataframe_rCN$Detected)
table(C3_dataframe_rCN$G, C3_dataframe_rCN$Detected)
table(C4_dataframe_rCN$G, C4_dataframe_rCN$Detected)