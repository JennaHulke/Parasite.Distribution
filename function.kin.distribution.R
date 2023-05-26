kinRelations <- function(file) {
  relation <- file
  ### get number of hosts
  a <- aggregate(relation$host.ID, list(relation$host.ID), FUN = length)
  host <- length(a$x)
  tots <- a$x
  
  ### set up dataset by host for pc and ps
  kin <- matrix(data = NA, nrow = host, ncol = 10)
  colnames(kin) <- c("host", "number of worms", "number of clonemates", "number of sibs", "total number pairwise", "mating_clones", "mating_siblings", "diff", "Pc", "Ps")
  kin[, 1] <- a$Group.1
  kin[, 2] <- a$x
  
  start_2 <- 1
  stop_2 <- kin[1, 2]
  
  for (j in 1:(nrow(kin)-1)) { ### loops through host
    b <- relation$Clonal.Group[start_2:stop_2]
    b_counts <- table(b)
    b_counts <- b_counts[b_counts > 1]
    kin[j, 3] <- sum(b_counts)
    
    c <- relation$Sibling.Group[start_2:stop_2]
    c_counts <- table(c)
    c_counts <- c_counts[c_counts > 1]
    kin[j, 4] <- sum(c_counts)
    
    start_2 <- stop_2 + 1
    stop_2 <- start_2 + kin[j + 1, 2] - 1
  }
  
  # Calculate values for the last row of 'kin'
  b <- relation$Clonal.Group[start_2:stop_2]
  b_counts <- table(b)
  b_counts <- b_counts[b_counts > 1]
  kin[host, 3] <- sum(b_counts)
  
  c <- relation$Sibling.Group[start_2:stop_2]
  c_counts <- table(c)
  c_counts <- c_counts[c_counts > 1]
  kin[host, 4] <- sum(c_counts)
  
  kin[, 5] <- (kin[, 2] * (kin[, 2] - 1)) / 2  ### counts total pairwise comparisons
  kin[, 6] <- (kin[, 3] * (kin[, 3] - 1)) / 2  ### counts clonal pairwise comparisons
  kin[, 7] <- (kin[, 4] * (kin[, 4] - 1)) / 2  ### counts clonal pairwise comparisons
  kin[, 8] <- kin[, 7] - kin[, 6]
  
  kin[, 9] <- (kin[, 6] / kin[, 5]) * kin[, 2] #### within host pc
  kin[, 10] <- (kin[, 8] / kin[, 5]) * kin[, 2]#### within host ps
  
  kin[is.na(kin)] <- 0
  host_pc <- sum(kin[, 9]) / sum(kin[, 2]) ### weighted average pc by host
  host_ps <- sum(kin[, 10]) / sum(kin[, 2])### weighted average ps by host
  
  pair.comp <- (sum(kin[,2]) *(sum(kin[,2])-1))/2
  total.clon.pair.comp <- (sum(kin[,3]) *(sum(kin[,3])-1))/2
  total.sib.pair.comp <- (sum(kin[,4]) *(sum(kin[,4])-1))/2
  sib.no.clones <- total.sib.pair.comp - total.clon.pair.comp
  
  total_pc <- total.clon.pair.comp/sum(kin[, 2])
  total_ps <- sib.no.clones/sum(kin[,2])
  return(list(host_pc = host_pc, host_ps = host_ps, total_pc = total_pc, total_ps = total_ps))
  
  
}

# Call the function with the appropriate argument
file <- read.csv("test.csv") ### mydata
  result <- kinRelations(file)

