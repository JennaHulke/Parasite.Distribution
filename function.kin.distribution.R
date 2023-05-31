kinRelations <- function(file) {
  relation <- file
  ### get number of hosts
  a <- aggregate(relation$host.ID, list(relation$host.ID), FUN = length)
  host <- length(a$x)
  tots <- a$x
  
  ### set up dataset by host for pc and ps
  kin <- matrix(data = NA, nrow = host, ncol = 13)
  colnames(kin) <- c("host", "number of worms", "number of clonemates", "number of sibs", "number unrelated", "total number pairwise", "mating_clones", "mating_siblings", "diff", "Pc", "Ps", "unrelated#","x")
  kin[, 1] <- a$Group.1
  kin[, 2] <- a$x
  
  start_2 <- 1
  stop_2 <- as.numeric(kin[1, 2])
  
  for (j in 1:(nrow(kin)-1)) { ### loops through host
    b <- relation$Clonal.Group[start_2:stop_2]
    b_counts <- as.numeric(table(b))
    b_counts <- b_counts[b_counts > 1]
    kin[j, 3] <- sum(b_counts)
    
    c <- relation$Sibling.Group[start_2:stop_2]
    c_counts <- as.numeric(table(c))
    c_counts <- c_counts[c_counts > 1]
    kin[j, 4] <- sum(c_counts)
    
    d <- relation$Sibling.Group[start_2:stop_2]
    d_counts <- as.numeric(table(d))
    d_counts <- d_counts[d_counts == 1]
    kin[j, 5] <- sum(d_counts)
    
    start_2 <- stop_2 + 1
    stop_2 <- start_2 + kin[j + 1, 2] - 1
  }
  
  # Calculate values for the last row of 'kin'
  b <- relation$Clonal.Group[start_2:stop_2]
  b_counts <- as.numeric(table(b))
  b_counts <- b_counts[b_counts > 1]
  kin[host, 3] <- sum(b_counts)
  
  c <- relation$Sibling.Group[start_2:stop_2]
  c_counts <- as.numeric(table(c))
  c_counts <- c_counts[c_counts > 1]
  kin[host, 4] <- sum(c_counts)
  
  d <- relation$Sibling.Group[start_2:stop_2]
  d_counts <- as.numeric(table(d))
  d_counts <- d_counts[d_counts == 1]
  kin[host, 5] <- sum(d_counts)
  
  kin[, 6] <- (kin[, 2] * (kin[, 2] - 1)) / 2  ### counts total pairwise comparisons
  kin[, 7] <- (kin[, 3] * (kin[, 3] - 1)) / 2  ### counts clonal pairwise comparisons
  kin[, 8] <- (kin[, 4] * (kin[, 4] - 1)) / 2  ### counts sibling pairwise comparisons
  kin[, 9] <- kin[, 8] - kin[, 7]
  kin[,9][kin[,9] < 0] <-0
  
  ##pc * intenstiy
  kin[, 10] <- (kin[, 7]/kin[, 6]) * kin[, 2] #### average within host pc times the intensity
  kin[, 11] <- (kin[, 9]/kin[, 6]) * kin[, 2]#### average within host ps
  
  kin[is.na(kin)] <- 0
  f <- c()
  f <- subset(kin[,2],kin[,2]>1)
  pc_weighted <- sum(kin[, 10]) / sum(f) ### weighted average pc by host with intensities greater 1
  ps_weighted <- sum(kin[, 11]) / sum(f) ### weighted average ps by host
  
  kin[,12] <- (kin[,5]*(kin[,2]-1))-((kin[,5]*(kin[,5]-1))/2)
  kin[,13] <- (kin[,5]*(kin[,4]))/((kin[,5]*(kin[,5]-1))/2)
  
  ###total pairwise comparisons of parasites of component pop
  pair.comp <- (sum(kin[,2]) *(sum(kin[,2])-1))/2 
  
  pec <- sum(kin[,7])/sum(f)
  return(list(pc_weighted = pc_weighted, ps_weighted = ps_weighted, pec = pec))
  
  
}
