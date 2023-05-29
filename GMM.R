# Load libraries
library(MixSim)           # Generate overlap parameters
library(mclust)           # k-means
library(cluster)          # pam
library(fclust)           # fuzzy k-means
library(FPDclustering)    # PD Clustering
library(ContaminatedMixt) # CND clustering
library(MixGHD)           # ARI


set.seed(1234)            # Random seed

##########################################
########## 1. Simulate Datasets ##########
##########################################

# Define variables
Pover <- c(0.1, 0.4) # proportion of overlap data
Pout <- c(0.05, 0.15)  # proportion of "outlier" data
n_datasets <- 20     # number of datasets to generate for each group

# Initialize list for final data
FinalData <- list()

# Loop through proportion of overlap data and "outlier" data
for (i in 1:2) {
  
  
  # Q <- MixSim(MaxOmega = Pover[i], K = 2, p = 2, hom = TRUE, sph = TRUE)
  
  for (j in 1:2) {
    
    # print(Q$Pi)
    
    # Loop through number of datasets
    for (m in 1:n_datasets) {
      
      # Generate mixture model and simulate dataset
      repeat{
        # homogeneous cluster; spherical covariance matrix structure
        Q <- MixSim(MaxOmega = Pover[i], K = 2, p = 2, hom = TRUE, sph = TRUE) 
        if (Q$fail == 0) break
      }
      
      A <- simdataset(n = 300 * (1 - Pout[j]), Pi = Q$Pi, Mu = Q$Mu, S = Q$S,
                      n.out = 300 *Pout[j], int = c(-1, 2), alpha = 0.01)
      
      colors <- c("red", "green", "blue")
      
      plot(A$X, xlab = "x1", ylab = "x2", type = "n") # A = dataset
      
      for (k in 0:2){
        points(A$X[A$id == k, ], col = colors[k+1], pch = k, cex = 0.8)
        
        points(x = Q$Mu[k,1], y = Q$Mu[k,2], col = colors[k+1], bg = colors[k+1], 
               pch = 20+k, cex = 1.2)
      }
      
      # Assign groups for outliers
      OutlierPerGroup = round(300 * (Pout[j]*0.3), 0)
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


#######################################
########## 2. Initialization ##########
#######################################

# Task: Use k-means, k-medoids (PAM), Fuzzy k-means, PD Clustering, at random to get initialization parameters
kmeans_cluster_labels <- rep(list(list()), 80)
pam_cluster_labels <- rep(list(list()), 80)
fkm_cluster_labels <- rep(list(list()), 80)
pdc_cluster_labels <- rep(list(list()), 80)
random_labels <- rep(list(list()), 80)

for (i in 1:80) {
  for (j in 2:5) {
    
    print(paste0("Processing ", names(FinalData)[i], ", using G=", j))
    
    kmeans_result <- kmeans(FinalData[[i]]$data$X, centers = j, nstart = 10) # Perform k-means clustering
    kmeans_cluster_labels[[i]][[j]] <- kmeans_result$cluster                 # Extract cluster labels from k-means
    
    
    pam_result <- pam(x=FinalData[[i]]$data$X, k=j)                           # Perform PAM/k-medoids clustering
    pam_cluster_labels[[i]][[j]] <- pam_result$clustering                    # Extract cluster labels from PAM/k-medoids
   
    
    fkm_result <- FKM(X=FinalData[[i]]$data$X, k=j)
    fkm_cluster_labels[[i]][[j]] <- unname(fkm_result$clus[,1])
    
    
    random_labels[[i]][[j]] <- sample(1:j, 300, replace=TRUE)
  }
}

for (i in 1:80) {
  for (g in 2:5) {
    
    print(paste0("Processing ", names(FinalData)[i], ", using G=", j))
    
    # PDC contains random initialization, which may lead to divided-by-zero error
    # Use a loop to catch the error and retry until success.
    retry_pdc = TRUE
    pdc_result = NULL
    while (retry_pdc) {
      tryCatch({pdc_result = PDC(FinalData[[i]]$data$X, k = g); retry_pdc = FALSE},
               warning = function(w) {}, error = function(e) {}, finally = {})
    }
    
    pdc_cluster_labels[[i]][[g]] <- pdc_result$label
    
  }
}

###############################################
########## 3. Gaussian Mixture Model ##########
###############################################

library(teigen)

# Results from K-means Initialization

kmeans.init.results = list()

print("K-means Gaussian Results Processing....")
for (i in 1:80) {
  
  print(paste0("Iteration: ", i))
  kmeans.init.results[[i]] <- teigen(FinalData[[i]]$data$X, Gs = 2:5,
                                           models = "UUUU", gauss = TRUE,
                                           init = kmeans_cluster_labels[[i]] )
  
}

# Results from Pam Initialization
pam.init.results = list()

print("PAM Gaussian Results Processing....")
for (i in 1:80) {
 
  print(paste0("Iteration: ", i))
  pam.init.results[[i]] <- teigen(FinalData[[i]]$data$X, Gs = 2:5, 
                                         models = "UUUU", gauss = TRUE,
                                         init = pam_cluster_labels[[i]])
  
}

# Results from Fuzzy C-Means Initialization
fkm.init.results = list()

print("Fuzzy C-Means Gaussian Results Processing....")
for (i in 1:80) {
  
  print(paste0("Iteration: ", i))
  fkm.init.results[[i]] <- teigen(FinalData[[i]]$data$X, Gs = 2:5,
                                         models = "UUUU", gauss = TRUE,
                                         init = fkm_cluster_labels[[i]])
  
}

# Results from PD Clustering Gaussian Initialization
pdc.init.results = list()

print("PD Clustering Init Results Processing....")
for (i in 1:80) {
  
  print(paste0("Iteration: ", i))
  pdc.init.results[[i]] = teigen(FinalData[[i]]$data$X, Gs = 2:5,
                                        models = "UUUU", gauss = TRUE,
                                        init = pdc_cluster_labels[[i]])
}


# Results from Random Initialization 
rand.init.results = list()

print("Random Init Gaussian Results Processing....")
for (i in 1:80) {
  
  print(paste0("Iteration: ", i))
  rand.init.results[[i]] <- teigen(FinalData[[i]]$data$X, Gs = 2:5, 
                                        models = "UUUU", gauss = TRUE,
                                        init = random_labels[[i]])
  
}


#######################################
########## 4. Finding ARI #############
#######################################

best.ARI.kmeans = list()
best.ARI.pam = list()
best.ARI.fkm = list()
best.ARI.pdc = list()
best.ARI.rand = list()

for (i in 1:80) {
  
  best.ARI.kmeans[[i]] = list(
    "ARI" = ARI(kmeans.init.results[[i]]$classification, FinalData[[i]]$data$id),
    "G" = kmeans.init.results[[i]]$G
  )
  
  best.ARI.pam[[i]] = list(
    "ARI" = ARI(pam.init.results[[i]]$classification, FinalData[[i]]$data$id),
    "G" = pam.init.results[[i]]$G
  )
  
  best.ARI.fkm[[i]] = list(
    "ARI" = ARI(fkm.init.results[[i]]$classification, FinalData[[i]]$data$id),
    "G" = fkm.init.results[[i]]$G
  )
  
  best.ARI.pdc[[i]] = list(
    "ARI" = ARI(pdc.init.results[[i]]$classification, FinalData[[i]]$data$id),
    "G" = pdc.init.results[[i]]$G
  )
  
  best.ARI.rand[[i]] = list(
    "ARI" = ARI(rand.init.results[[i]]$classification, FinalData[[i]]$data$id),
    "G" = rand.init.results[[i]]$G
  )
}


#####################################################
########## 5. Save ARI values in matrix #############
#####################################################

# dataset with 10% overlap and 5% outliers
ARI.matrix1 = matrix(NA, 20, 5)
# dataset with 10% overlap and 15% outliers
ARI.matrix2 = matrix(NA, 20, 5)
# dataset with 40% overlap and 5% outliers
ARI.matrix3 = matrix(NA, 20, 5)
# dataset with 40% overlap and 15% outliers
ARI.matrix4 = matrix(NA, 20, 5)

for (i in 1:20) {
  # dataset with 10% overlap and 5% outliers
  ARI.matrix1[i, ] = c(best.ARI.kmeans[[i]]$ARI, best.ARI.pam[[i]]$ARI, best.ARI.fkm[[i]]$ARI,
                      best.ARI.pdc[[i]]$ARI, best.ARI.rand[[i]]$ARI) 
  
  # dataset with 10% overlap and 15% outliers
  ARI.matrix2[i, ] = c(best.ARI.kmeans[[i+20]]$ARI, best.ARI.pam[[i+20]]$ARI, best.ARI.fkm[[i+20]]$ARI,
                       best.ARI.pdc[[i+20]]$ARI, best.ARI.rand[[i+20]]$ARI) 
  
  # dataset with 40% overlap and 5% outliers
  ARI.matrix3[i, ] = c(best.ARI.kmeans[[i+40]]$ARI, best.ARI.pam[[i+40]]$ARI, best.ARI.fkm[[i+40]]$ARI,
                       best.ARI.pdc[[i+40]]$ARI, best.ARI.rand[[i+40]]$ARI) 
  
  # dataset with 40% overlap and 15% outliers
  ARI.matrix4[i,] = c(best.ARI.kmeans[[i+60]]$ARI, best.ARI.pam[[i+60]]$ARI, best.ARI.fkm[[i+60]]$ARI,
                      best.ARI.pdc[[i+60]]$ARI, best.ARI.rand[[i+60]]$ARI) 
}

# plot boxplot
df1 = as.data.frame(ARI.matrix1)
df2 = as.data.frame(ARI.matrix2)
df3 = as.data.frame(ARI.matrix3)
df4 = as.data.frame(ARI.matrix4)
colnames(df1) = c("K-Means", "PAM", "FKM", "PDC", "Random")
colnames(df2) = c("K-Means", "PAM", "FKM", "PDC", "Random")
colnames(df3) = c("K-Means", "PAM", "FKM", "PDC", "Random")
colnames(df4) = c("K-Means", "PAM", "FKM", "PDC", "Random")

# create a 4 x 4 frame with each cell containing a plot
par(mfrow=c(2,2))

boxplot(df1, main="Data: Overlap 10% and Outlier 5% (n=300)",
        ylab= "Adjusted Rand Index", xlab="Initialization Algorithms",
        col = rainbow(5))

boxplot(df2, main="Data: Overlap 10% and Outlier 15% (n=300)",
        ylab= "Adjusted Rand Index", xlab="Initialization Algorithms",
        col = rainbow(5))

boxplot(df3, main="Data: Overlap 40% and Outlier 5% (n=300)",
        ylab= "Adjusted Rand Index", xlab="Initialization Algorithms",
        col = rainbow(5))

boxplot(df4, main="Data: Overlap 40% and Outlier 15% (n=300)",
        ylab= "Adjusted Rand Index", xlab="Initialization Algorithms",
        col = rainbow(5))


# G
G.matrix1 <- matrix(NA, 20, 5)
G.matrix2 <- matrix(NA, 20, 5)
G.matrix3 <- matrix(NA, 20, 5)
G.matrix4 <- matrix(NA, 20, 5)

for (i in 1:20) {
  G.matrix1[i, ] = c(best.ARI.kmeans[[i]]$G, best.ARI.pam[[i]]$G, best.ARI.fkm[[i]]$G,
                       best.ARI.pdc[[i]]$G, best.ARI.rand[[i]]$G) 
  G.matrix2[i, ] = c(best.ARI.kmeans[[i+20]]$G, best.ARI.pam[[i+20]]$G, best.ARI.fkm[[i+20]]$G,
                       best.ARI.pdc[[i+20]]$G, best.ARI.rand[[i+20]]$G) 
  G.matrix3[i, ] = c(best.ARI.kmeans[[i+40]]$G, best.ARI.pam[[i+40]]$G, best.ARI.fkm[[i+40]]$G,
                       best.ARI.pdc[[i+40]]$G, best.ARI.rand[[i+40]]$G) 
  G.matrix4[i,] = c(best.ARI.kmeans[[i+60]]$G, best.ARI.pam[[i+60]]$G, best.ARI.fkm[[i+60]]$G,
                      best.ARI.pdc[[i+60]]$G, best.ARI.rand[[i+60]]$G) 
}

# create a matrix: recording method, ARI, and G
# Scenario 1
S1_dataframe = data.frame(
  initial_method = c(rep("K-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(ARI.matrix1),
  G = c(G.matrix1))

# Scenario 2
S2_dataframe = data.frame(
  initial_method = c(rep("K-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(ARI.matrix2),
  G = c(G.matrix2))

# Scenario 3
S3_dataframe = data.frame(
  initial_method = c(rep("k-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(ARI.matrix3),
  G = c(G.matrix3))

# Scenario 4\
S4_dataframe = data.frame(
  initial_method = c(rep("K-means", 20), rep("PAM", 20), rep("FKM", 20), rep("PDC", 20), rep("Random", 20)),
  ARI = c(ARI.matrix4),
  G = c(G.matrix4))


library(ggplot2)          # Plot
library(dplyr)            # Plot
library(ggpubr)           # Plot

C1_group_dataframe = as.data.frame(table(S1_dataframe$initial_method, S1_dataframe$G))
C2_group_dataframe = as.data.frame(table(S2_dataframe$initial_method, S2_dataframe$G))
C3_group_dataframe = as.data.frame(table(S3_dataframe$initial_method, S3_dataframe$G))
C4_group_dataframe = as.data.frame(table(S4_dataframe$initial_method, S4_dataframe$G))


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
  theme(plot.title=element_text(hjust=0.5,color = "red", size = 14, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 1: Overlap 10% & Outlier 5% | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C2_G <- ggplot(C2_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "red", size = 14, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 2: Overlap 10% & Outlier 15% | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C3_G <- ggplot(C3_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "red", size = 14, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 3: Overlap 40% & Outlier 5% | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

plot_C4_G <- ggplot(C4_group_dataframe, aes(fill=G, y=Frequency, x=Initialization)) + 
  geom_bar(position="stack", stat="identity") +
  theme(plot.title=element_text(hjust=0.5,color = "red", size = 14, face = "bold"),
        axis.text = element_text(face="bold", size = 12),
        axis.title=element_text(size=14,face="bold"))  +
  xlab("Initialization Method") +
  ylab("ARI") +
  ggtitle("Scenario 4: Overlap 40% & Outlier 15% | n=300") +
  scale_fill_brewer(palette="Paired") +
  geom_text(aes(label = Frequency/20), size = 5, hjust = 0.5, vjust = 1.5, position = "stack")

ggarrange(plot_C1_G, plot_C2_G, plot_C3_G, plot_C4_G, ncol=2, nrow = 2)




