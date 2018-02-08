#' Function to use the solution of mixture modelling to phase the BAF
#'
#' This function is primarily used within the MiMMAl function
#' @param x This is a list element created in the MiMMAl function.
#' @param seed This is the seed used when mirroring individual values. Default = 1.
#' @export
#' @examples
#' getBAFphased()

##  This function takes the result of a mixture model and flips the BAF according to the probability that a given
##    BAF is part of a distribution with a mean < 0.5.
##
##  Coded by George Cresswell; 2016-11-08
getBAFphased = function(x, seed=1) {

  #x is the mix object
  mix.object <- x

  #This is only for two mixtures
  k <- 2

  #Sort the BAFs and make them x now
  x <- mix.object$x

  #Make an empty object to collect this
  densofdists = list()

  #Get the density
  for (i in 1:k) {

    #What is the density function for this distribution?
    densofdists[[i]] = dnorm(x,
                             mean = mix.object$mu[i],
                             sd = mix.object$sigma[i])

  }

  #Which means are less than 0.5?
  flip.dists = mix.object$mu < 0.5

  #Now what do we return?
  if(all(flip.dists)) {

    #Make a matrix tell us to flip 'em all!
    probs = cbind(x, rep(1, times = length(x)))

  } else if(all(!flip.dists)) {

    #Make a matrix tell us to do nuffin'
    probs = cbind(x, rep(0, times = length(x)))

  } else {

    #Now for the bread and bu'er
    flip.dists = order(flip.dists, decreasing = TRUE)

    #Now make the matrix we want to return
    prob.flip  = densofdists[[flip.dists[1]]] / (densofdists[[flip.dists[1]]]+densofdists[[flip.dists[2]]])

    #What are the probs?
    probs = cbind(x, prob.flip)

  }

  set.seed(seed)
  runifs = runif(length(x))

  probs = cbind(probs, runifs)

  flippedBAF = apply(probs, 1, function(x) {

    if(x[2] > x[3]) {

      return(1 - x[1])

    } else {return(x[1])}

  })

  return(flippedBAF)

}
