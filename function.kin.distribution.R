kinRelations <- function(file) {
  relation <- file
  ### shows number of parasites found in each host/gland
  a <- aggregate(relation$host.ID, list(relation$host.ID), FUN = length)
  host <- length(a$x) ###number of glands/hosts
  tots <- a$x ##intensity by host
  
  ### set up dataset by host for pc and ps
  kin <- matrix(data = NA, nrow = host, ncol = 15)
  colnames(kin) <- c("host", "number of worms", "number of clonemates", "number of sibs", 
                     "number unrelated", "total number pairwise", "Clonemate.dyads", "faux_sib.dyads", 
                     "Sib.dyads", "Pc", "Ps", "unrelated#","x", "unrelated.clonemates.sibs", "14equalstotal")
  kin[, 1] <- a$Group.1 ##host/gland number
  kin[, 2] <- a$x ###intensity per host/gland
  
  start_2 <- 1
  stop_2 <- as.numeric(kin[1, 2])
  ###Identifies identical Clonemates found within a host
  for (j in 1:(nrow(kin)-1)) { ### loops through host using the intensity of infection
    b <- relation$Clonal.Group[start_2:stop_2] #lists parasite ID found in same host
    b_counts <- as.numeric(table(b))## IDs repeated parasite IDs (clonemates) found in a host
    b_counts <- b_counts[b_counts > 1]### IDs when there are multiple parasites with same IDs
    kin[j, 3] <- sum(b_counts) ##fills kin matrix
    ###Identifies full-siblings found within the same host
    c <- relation$Sibling.Group[start_2:stop_2]
    c_counts <- as.numeric(table(c))## IDs repeated parasite IDs (siblings) found in a host
    c_counts <- c_counts[c_counts > 1]### IDs when there are multiple parasites with same IDs
    kin[j, 4] <- sum(c_counts)##fills kin matrix
    ###Identifies unrelated individuals found within a host
    d <- relation$Sibling.Group[start_2:stop_2]
    d_counts <- as.numeric(table(d))
    d_counts <- d_counts[d_counts == 1]
    kin[j, 5] <- sum(d_counts)
    
    start_2 <- stop_2 + 1 #moves to next host ID to "start"
    stop_2 <- start_2 + kin[j + 1, 2] - 1 #identifies host many parasites are in that host to "stop"
  }

  kin[, 6] <- (kin[, 2] * (kin[, 2] - 1)) / 2  ### counts total pairwise comparisons
  kin[, 7] <- (kin[, 3] * (kin[, 3] - 1)) / 2  ### counts clonal pairwise comparisons
  kin[, 8] <- (kin[, 4] * (kin[, 4] - 1)) / 2  ### counts faux sib dyads (includes clonemates) pairwise comparisons
  kin[, 9] <- kin[, 8] - kin[, 7] ###takes faux sib pairwise comparisons and minuses clonemates pairwise comparsions to get true sib dyads
  
  ##pc * intenstiy
  kin[, 10] <- (kin[, 7]/kin[, 6]) * kin[, 2] #### average within host pc times the intensity
  kin[, 11] <- (kin[, 9]/kin[, 6]) * kin[, 2]#### average within host ps
  
  kin[is.na(kin)] <- 0
  f <- c()
  f <- subset(kin[,2],kin[,2]>1)
  pc_weighted <- sum(kin[, 10]) / sum(f) ### weighted average pc by host with intensities greater 1
  ps_weighted <- sum(kin[, 11]) / sum(f) ### weighted average ps by host
  
  kin[,12] <- (kin[,5]*(kin[,2]-1))-((kin[,5]*(kin[,5]-1))/2)##calculates the number of unrelated
  kin[,13] <- (kin[,5]*(kin[,4]))+((kin[,5]*(kin[,5]-1))/2)##serves as a check for the number of unrelated
  
  kin[,14] <- kin[,13]+kin[,9]+kin[,7]#counts total number of pairwise comparisons from results and checks with initial count 
  kin[,15] <- kin[,14]==kin[,6] ## summarizes kin[,14] with 1 being true, 0 being false
  ###total pairwise comparisons of parasites of component pop
  pair.comp <- (sum(kin[,2]) *(sum(kin[,2])-1))/2 
  
  ###calculates pec
  e <-relation$Clonal.Group
  e_counts <- as.numeric(table(e))##counts the number of time a parasite ID is repeated (clonemates)
  tot.clone <- e_counts[e_counts > 1]##reduces down to number of repeated IDs
  tot.clone.pairwise <- sum(tot.clone*(tot.clone -1)/2)##counts number of pairwise comparsions from reduced
  pec <- tot.clone.pairwise/(sum(kin[,2])*(sum(kin[,2])-1)/2)##calculates pec
  ###calculates pes
  g <-relation$Sibling.Group
  g_counts <- as.numeric(table(g))##counts the number of time a parasite ID is repeated (siblings)
  tot.sib <- g_counts[g_counts > 1]##reduces down to number of repeated IDs
  tot.sib.pairwise <- sum(tot.sib*(tot.sib -1)/2)##counts number of pairwise comparsions from reduced
  pes <- tot.sib.pairwise/(sum(kin[,2])*(sum(kin[,2])-1)/2)##calculates pes
  
  return(list(pc_weighted = pc_weighted, ps_weighted = ps_weighted, pec = pec, pes = pes))
  
  
}
