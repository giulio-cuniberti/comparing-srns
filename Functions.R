## libraries #######################################################################################

library(pracma)
library(lpSolve)



## auxiliary functions #############################################################################

populateAB <- function (conservation.laws, preorder.matrix) {
  m <- nrow(preorder.matrix)  # number of inequalities
  d <- ncol(preorder.matrix)  # number of species
  if (ncol(conservation.laws) != d) {
    stop('number of columns of matrices "conservation.laws" and "preorder.matrix" does not match')
  }
  A <- matrix(rep(0, (m+1)*d), nrow = m+1, ncol = d)
  B <- matrix(rep(0, (m+1)*d), nrow = m+1, ncol = d)
  N <- rbind(preorder.matrix, rep(0, d))  # preorder.matrix with extra 0 row
  # compute the matrices
  for (i in (m+1):1) {
    for (j in 1:d) {
      D <- cbind(t(conservation.laws), -t(conservation.laws), t(preorder.matrix), -N[i, ])
      e <- rep(0, d)
      e[j] <- 1
      if (A[m+1, j] == 1) {
        A[i, j] <- 1  # no need to solve lp problem in this case
      } else if (lp('min', rep(0, ncol(D)), D, rep("=", d), e)$status == 0) {
        A[i, j] <- 1
      }
      if (B[m+1, j] == 1) {
        B[i, j] <- 1  # no need to solve lp problem in this case
      } else if (lp('min', rep(0, ncol(D)), D, rep("=", d), -e)$status == 0) {
        B[i, j] <- 1
      }
    }
  }
  return(list(A = A, B = B))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

checkConditions <- function(source.complexes, reaction.vectors, conservation.laws, 
                            preorder.matrix) {
  n <- nrow(source.complexes)  # number of reactions
  d <- ncol(source.complexes)  # number of species
  m <- nrow(preorder.matrix)  # number of inequalities
  if (nrow(reaction.vectors) != n) {
    stop('number of rows of matrices "source.complexes" and "reaction.vectors" does not match')
  }
  if ((ncol(reaction.vectors) != d) || (ncol(conservation.laws) != d) || 
      (ncol(preorder.matrix) != d)) {
    stop('number of columns of matrices "source.complexes", "reaction.vectors", "conservation.laws" 
          and "preorder.matrix" does not match')
  }
  # compute matrices A and B
  AandB <- populateAB(conservation.laws, preorder.matrix)
  A <- AandB$A
  B <- AandB$B
  # create the output variables
  result <- 1
  ratesIneq <- matrix(rep(0, 2*n), nrow = 2)
  speciesIneq <- rbind(A[m+1, ], -B[m+1, ])
  firstFail <- ''
  # check conditions for each reaction
  for (r in 1:n) {
    vec <- preorder.matrix %*% reaction.vectors[r, ]
    pos <- vec*(vec>0)
    neg <- -vec*(vec<0)
    check <- 0
    # conditions (a)
    if (all(pos == 0)) {
      check <- 1
    } else if (all(A[m+1, ]*source.complexes[r, ] == source.complexes[r, ])) {
      check <- 1
      ratesIneq[1, r] <- 1
    } else {
      for (i in 1:m) {
        e <- rep(0, m)
        e[i] <- 1
        if (all(pos == e) && all(A[i, ]*source.complexes[r, ] == source.complexes[r, ])) {
          check <- 1
          ratesIneq[1, r] <- 1
        }
      }
    }
    if (check != 1) {
      result <- 0
      firstFail <- paste('reaction', r, '/ conditions (a)')
      break
    }
    # conditions (b)
    if (all(neg == 0)) {
      check <- 2
    } else if (all(B[m+1, ]*source.complexes[r, ] == source.complexes[r, ])) {
      check <- 2
      ratesIneq[2, r] <- -1
    } else {
      for (i in 1:m) {
        e <- rep(0, m)
        e[i] <- 1
        if (all(neg == e) && all(B[i, ]*source.complexes[r, ] == source.complexes[r, ])) {
          check <- 2
          ratesIneq[2, r] <- -1
        }
      }
    }
    if (check != 2) {
      result <- 0
      firstFail <- paste('reaction', r, '/ conditions (b)')
      break
    }
  }
  return(list(result = result, rate.inequalities = ratesIneq, 
              species.inequalities = speciesIneq, first.fail = firstFail))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

populateRC <- function(source.complexes, product.complexes) {
  n <- nrow(source.complexes)  # number of reactions
  d <- ncol(source.complexes)  # number of species
  if (nrow(product.complexes) != n) {
    stop('number of rows of matrices "source.complexes" and "product.complexes" does not match')
  }
  if ((ncol(product.complexes) != d)) {
    stop('number of columns of matrices "source.complexes" and "product.complexes" does not match')
  }
  reaction.vectors <- product.complexes-source.complexes         # matrix of reaction vectors
  Ct <- null(reaction.vectors)    # matrix of conservation laws by columns
  # matrix of conservation laws by rows (or 0 row if there are no conservation laws)
  if (is.null(Ct)) {
    conservation.laws <- matrix(rep(0, d), nrow = 1)
  } else {
    conservation.laws <- t(Ct)
  }
  # rounding should prevent NUMERICAL ISSUES (such as getting 1e-16 instead of 0) 
  conservation.laws <- round(conservation.laws, digits = 10)
  return(list(R = reaction.vectors, C = conservation.laws))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

toBinary <- function(decimal) {
  if (decimal == 0) {
    return('0')
  }
  binary <- ''
  while (decimal > 0) {
    binary <- paste0(decimal %% 2, binary)
    decimal <- decimal %/% 2
  }
  return(binary)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

toSymbols <- function(inequalities.matrix) {
  vec <- c()
  n <- ncol(inequalities.matrix)
  for (i in 1:n) {
    if (all(inequalities.matrix[, i] == c(1,-1))) {
      vec[i] <- '='
    } else if (all(inequalities.matrix[, i] == c(1,0))) {
      vec[i] <- '+'
    } else if (all(inequalities.matrix[, i] == c(0,-1))) {
      vec[i] <- '-'
    } else if (all(inequalities.matrix[, i] == c(0,0))) {
      vec[i] <- '?'
    } else {
      vec[i] <- 'ERROR'
    }
  }
  return(t(as.matrix(vec)))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

processResults <- function(orderings) {
  preorders <- list()
  equivalences <- list()
  orderings <- unique(orderings)
  l <- length(orderings)
  if (l == 0) {
    output <- list(preorders = preorders, equivalences = equivalences)
    return(output)
  } else {
    n <- ncol(orderings[[1]]$rate.inequalities)
    d <- ncol(orderings[[1]]$species.inequalities)
  }
  for (i in 1:(l-1)) {
    reverseIneq <- 0
    stop <- 0
    k <- 1
    while ((k <= n+d) && (stop == 0)) {
      if (all(cbind(orderings[[i]]$rate.inequalities, 
                    orderings[[i]]$species.inequalities)[, k] == c(0,-1))) {
        reverseIneq <- 1
        stop <- 1
      } else if (all(cbind(orderings[[i]]$rate.inequalities, 
                           orderings[[i]]$species.inequalities)[, k] == c(1,0))) 
      {
        stop <- 1
      }
      k <- k+1
    }
    if (reverseIneq == 1) {
      orderings[[i]]$rate.inequalities <- as.matrix(-orderings[[i]]$rate.inequalities[c(2,1), ])
      orderings[[i]]$species.inequalities <- 
        as.matrix(-orderings[[i]]$species.inequalities[c(2,1), ])
    }
    for (j in (i+1):l) {
      if (all(orderings[[j]]$rate.inequalities == -orderings[[i]]$rate.inequalities[c(2,1), ])) {
        orderings[[j]]$rate.inequalities <- -orderings[[j]]$rate.inequalities[c(2,1), ]
        orderings[[j]]$species.inequalities <- -orderings[[j]]$species.inequalities[c(2,1), ]
      }
    }
  }
  orderings <- unique(orderings)
  l <- length(orderings)
  ordersMatrix <- matrix(nrow = 2*(n+d), ncol = l)
  for (k in 1:l) {
    ordersMatrix[, k] <- c(abs(orderings[[k]]$rate.inequalities), 
                           abs(orderings[[k]]$species.inequalities))
  }
  sortedIndexes <- do.call(order, c(as.data.frame(t(ordersMatrix[1:(2*(n+d)), ])), 
                                    list(decreasing = TRUE)))
  orderings <- orderings[sortedIndexes]
  orderings[[l]] <- NULL
  orderings[[1]] <- NULL
  l <- length(orderings)
  if (l == 0) {
    output <- list(preorders = preorders, equivalences = equivalences)
    return(output)
  }
  for (k in 1:l) {
    orderings[[k]]$rate.inequalities <- toSymbols(orderings[[k]]$rate.inequalities)
    orderings[[k]]$species.inequalities <- toSymbols(orderings[[k]]$species.inequalities)
  }
  for (k in 1:l) {
    vec <- orderings[[k]]$species.inequalities
    if (any(vec == '+' | vec == '-')) {
      preorders <- append(preorders, list(orderings[[k]]))
    } else {
      equivalences <- append(equivalences, list(orderings[[k]]))
    }
  }
  output <- list(preorders = preorders, equivalences = equivalences)
  return(output)
}



## main functions ##################################################################################

checkOrdering <- function(source.complexes, product.complexes, preorder.matrix) {
  # compute matrices R and C
  RandC <- populateRC(source.complexes, product.complexes)
  reaction.vectors <- RandC$R
  conservation.laws <- RandC$C
  # everything can be now computed with the function checkConditions()
  output <- checkConditions(source.complexes, reaction.vectors, conservation.laws, preorder.matrix)
  output$rate.inequalities <- toSymbols(output$rate.inequalities)
  output$species.inequalities <- toSymbols(output$species.inequalities)
  if (output$result == 1) {
    cat('\n YES: the trajectories can be ordered almost surely as the rate constants vary \n\n')
  } else {
    cat('\n NO: the trajectories cannot be ordered almost surely as the rate constants vary \n\n')
  }
  return(output)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

findOrderings <- function(source.complexes, product.complexes) {
  n <- nrow(source.complexes)  # number of reactions
  d <- ncol(source.complexes)  # number of species
  orderings <- list()
  # compute matrices reaction.vectors and conservation.laws
  RandC <- populateRC(source.complexes, product.complexes)
  reaction.vectors <- RandC$R
  conservation.laws <- RandC$C
  # try all 'simple' and essentially different orders
  E <- matrix(nrow = 2*d, ncol = d)
  for (i in 1:d) {
    e <- rep(0, d)
    e[i] <- 1
    E[i, ] <- e
    E[d+i, ] <- -e
  }
  cat('\n\n')
  if (d <= 6) {steps <- 10}
  else if (d == 7) {steps <- 20}
  else if (d == 8) {steps <- 50}
  else if (d >= 9) {steps <- 100}
  t0 <- Sys.time()
  cat(paste0(c('> PROGRESS ', rep(' ', steps-5), '%      TIME NOW     END TIME\n'), collapse = ''))
  elapsed <- 0
  cat(paste0(c('[', rep('-', steps), ']    ', 0, '      ',
               format(Sys.time(), "%H:%M:%S"), '     ESTIMATE\n'), collapse = ''))
  lastPrinted <- 0
  for (i in (2^d-1):0) {
    for (j in i:0) {
      coef1 <- as.integer(strsplit(toBinary(i), NULL)[[1]])
      l1 <- length(coef1)
      coef2 <- as.integer(strsplit(toBinary(j), NULL)[[1]])
      l2 <- length(coef2)
      coeff <- rep(0, 2*d)
      coeff[(d-l1+1):d] <- coef1
      coeff[(2*d-l2+1):(2*d)] <- coef2
      preorder.matrix <- diag(coeff) %*% E
      preorder.matrix <- preorder.matrix[apply(preorder.matrix, 1, function(x) !all(x == 0)), ]
      if (!is.matrix(preorder.matrix)) {
        preorder.matrix <- matrix(preorder.matrix, ncol = d)
      }
      check <- checkConditions(source.complexes, reaction.vectors, conservation.laws, 
                               preorder.matrix)
      if (check$result == 1) {
        orderings[[length(orderings)+1]] <- list(rate.inequalities = check$rate.inequalities, 
                                                 species.inequalities = check$species.inequalities)
      }
      progress <- floor(steps*(sum((2^d):(i+1))-j)/sum((2^d):1))
      if (progress > lastPrinted) {
        elapsed <- Sys.time() - t0
        remaining <- elapsed/progress*(steps-progress)
        endtime <- Sys.time() + remaining
        if (100*progress/steps < 10) {
          cat(paste0(c('[', rep('|', progress), rep('-', steps-progress),
                       ']    ', 100*progress/steps, '      ', 
                       format(Sys.time(), "%H:%M:%S"), '     ', 
                       format(endtime, "%H:%M:%S"), '\n'), collapse = ''))
        } else if (100*progress/steps == 100) {
          cat(paste0(c('[', rep('|', progress), rep('-', steps-progress),
                       ']    ', 100*progress/steps, '    ', 
                       format(Sys.time(), "%H:%M:%S"), '     ', 
                       format(endtime, "%H:%M:%S"), '\n'), collapse = ''))
        } else {
          cat(paste0(c('[', rep('|', progress), rep('-', steps-progress),
                       ']    ', 100*progress/steps, '     ', 
                       format(Sys.time(), "%H:%M:%S"), '     ', 
                       format(endtime, "%H:%M:%S"), '\n'), collapse = ''))
        }
        lastPrinted <- progress
      }
    }
  }
  cat('\n')
  orderings <- processResults(orderings)
  cat(paste0(length(orderings$preorders), ' PREORDERS (', 2*length(orderings$preorders), 
             ' with their opposites) and ', length(orderings$equivalences), 
             ' EQUIVALENCES were found\n\n\n'))
  return(orderings)
}


