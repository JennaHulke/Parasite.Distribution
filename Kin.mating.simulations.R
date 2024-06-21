###simulations for running function.kin.distribution
### this takes the observed data, creates random distribution of clones/siblings
### with given intensities per host
### simulations done without replacement

obs.data <- read.csv("test.csv")
source("function.kin.distribution.R") #required to run the kinRelations function
### get pc and ps for observed data
results_obs <- kinRelations(obs.data)

sample.host <- obs.data[,1]##IDs number of parasites found within host by host ID, repeated numbers are multiple parasites within the host

sim.data <- data.frame(matrix(data= NA, nrow= nrow(obs.data), ncol = 4))
colnames(sim.data) <-c("host.ID", "parasite.ID", "Clonal.Group", "Sibling.Group")
sim.data[,2] <- obs.data[,2]
sim.data[,3] <- obs.data[,3]
sim.data[,4] <- obs.data[,4]

##bootstraps from observed data
sim.summary <- matrix(data= NA, nrow= 10000, ncol = 4)
colnames(sim.summary) <- c("pc_weighted", "ps_weighted", "pec", "pes")
for (j in 1:10000){
  ##samples observed parasites and reassigns them to observed hosts
    host <- sample(sample.host, 
                       size = nrow(obs.data), 
                       replace = FALSE)
    sim.data[,1] <- host ##keeps parasite ID, clone and sib ID, reassigns these to a host
    sim.data <-sim.data[order(sim.data[,1],decreasing=FALSE),]##reorders
    sim.data <- as.data.frame(sim.data)
    results <- kinRelations(sim.data)#runs kinRelations function
    #write.csv(sim.data, file = paste("Sim_data/sim", j, ".csv", sep =""))
    sim.summary[j,1] <- results[[1]]
    sim.summary[j,2] <- results[[2]]
    sim.summary[j,3] <- results[[3]]
    sim.summary[j,4] <- results[[4]]
    
    }
mean(sim.summary[,1])  
mean(sim.summary[,2])
mean(sim.summary[,3])
mean(sim.summary[,4])

hist(sim.summary[,1])
hist(sim.summary[,2])
hist(sim.summary[,3])
