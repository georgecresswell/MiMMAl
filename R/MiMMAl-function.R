#' The main MiMMAl function.
#'
#' This function takes segmented BAF data as input and models the major and minor allele distributions.
#' @param samplename The output name for your sample.
#' @param inputfile The file to be read in as containing the columns; chr, pos, BAF, BAFseg.
#' @param min.snps Minimum number of SNPs to model. Default = 10.
#' @param sd.width This value is multiplied by the chosen sd and the tested sd range is determined as the sd±sd.width. Default = 1/3.
#' @param preset.sd Add a float value for the standard deviation if it is to be preset, without this it 'learns' the sd. Default = NULL.
#' @param seed Seed set for reproducibility. Default = 1.
#' @param baf.res 10^-baf.res is used as the size of the steps used in the grid search for testing BAF values. Default = 2.
#' @param use.ks.gate Require the data to pass a normality test before modelling. Default = TRUE.
#' @param plot.sd.den Plot the sd density plot? Default = TRUE.
#' @param plot.lit.plot Plot the global grid search? Default = TRUE.
#' @param plot.star.plot Plot the local grid search? Default = TRUE.
#' @param plot.transformed Plot the transformed data? Default = TRUE.
#' @import mixtools
#' @import ggplot2
#' @import cowplot
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#' @importFrom cowplot ggsave
#' @export
#' @examples
#' runMiMMAl()

#Define MiMMAl function
runMiMMAl = function(samplename,
                     inputfile,
                     min.snps    = 10,
                     sd.width    = 1/3,
                     preset.sd   = NULL,
                     seed        = 1,
                     baf.res     = 2,
                     use.ks.gate      = TRUE,
                     plot.sd.den      = TRUE,
                     plot.lit.plot    = TRUE,
                     plot.star.plot   = TRUE,
                     plot.transformed = TRUE) {

  #Read in the data
  BAFdata  = read.table(inputfile,
                        header = TRUE,
                        stringsAsFactors = FALSE)

  #Name the output file
  outputfile=paste0(samplename,".BAFphased.txt")

  #Add option to proceed with a preset standard deviation
  if(is.null(preset.sd)) {

    #Record the standard deviations calculated
    sds    = NULL

    #Run through the chromosomes to segment and model them to find the standard deviation of the array
    for (chr in unique(BAFdata[,1])) {

      BAFrawchr = BAFdata[BAFdata[,1]==chr,c(2,3)]
      BAFrawchr = BAFrawchr[!is.na(BAFrawchr[,2]),]

      BAF = BAFrawchr[,2]
      pos = BAFrawchr[,1]
      names(BAF) = rownames(BAFrawchr)
      names(pos) = rownames(BAFrawchr)

      print(paste("BAFlen=",length(BAF),sep=""))

      #
      BAFsegm = BAFdata[BAFdata[,1]==chr,c(4)]

      #Now take the positions of segment boundaries
      starts         = c(1, which(BAFsegm[-1] != BAFsegm[-length(BAFsegm)])+1)
      ends           = c(which(BAFsegm[-1] != BAFsegm[-length(BAFsegm)]), length(BAFsegm))
      segment.bounds = cbind(starts, ends)

      #Run through the segments and model them
      for(seg in 1:length(starts)) {

        #Take the BAF in this segment
        BAF.seg = BAF[starts[seg]:ends[seg]]

        #Only bother modelling if it greater than 10 SNPs, otherwise just take the mean
        if(length(BAF.seg) > min.snps) {

          #Set the smu
          smu = c(-mean(abs(BAF.seg-0.5)), mean(abs(BAF.seg-0.5)))

          #Model it as a mixture of two normal distributions with means that are reflected and the same sds
          set.seed(seed)
          mixmdl = normalmixEM(BAF.seg-0.5,
                               mu = smu,
                               maxrestarts = 10000,
                               verb = FALSE,
                               mean.constr = c("a", "-a"),
                               sd.constr = c("a", "a"))

          #Record the calculated standard deviation
          sds = c(sds, mixmdl$sigma[1])

          #Show me what you got!
          print(paste0("Standard deviation is: ",paste0(round(mixmdl$sigma[1], digits = 3), collapse = " ")))

        }

      }

    }

    #Calculate the densities of the standard deviations
    sds.density  = density(sds)

    #Now calculated the standard deviation of the signal seen in the array
    array.sd     = sds.density$x[which.max(sds.density$y)]

    #Make a dataframe for ggplot because it always fucking insists on it
    sds.den.df   = data.frame(sd=sds.density$x, density=sds.density$y)

    if(plot.sd.den) {
      #Output a png of the standard deviations
      png(paste0("Density_of_standard_deviation_",samplename,".png"), width = 600, height = 600)
      p = ggplot(sds.den.df, aes(x=sd, y=density)) +
        geom_line() +
        geom_vline(xintercept=array.sd, color="gray", linetype = "longdash") +
        geom_vline(xintercept=c(array.sd*(1-sd.width), array.sd*(1+sd.width)), color="red", linetype = "longdash") +
        ggtitle(paste0("Density of standard deviation ",samplename))
      print(p)
      dev.off()
    }

  } else {array.sd = preset.sd}

  #So what range of standard deviations will we use in the grid search?
  sigmai     = seq(round(array.sd*(1-sd.width), digits=3),
                   round(array.sd*(1+sd.width), digits=3),
                   length.out = 21)

  #What is the baf search range?
  bafsearch  = seq(0.5,
                   0,
                   by=-10^-baf.res)

  #Make a variable for recording BAFoutput
  BAFoutput = NULL

  #Run through the chromosomes to segment and model them
  for (chr in unique(BAFdata[,1])) {

    BAFrawchr = BAFdata[BAFdata[,1]==chr,c(2,3)]
    BAFrawchr = BAFrawchr[!is.na(BAFrawchr[,2]),]

    BAF = BAFrawchr[,2]
    pos = BAFrawchr[,1]
    names(BAF) = rownames(BAFrawchr)
    names(pos) = rownames(BAFrawchr)

    #
    BAFsegm = BAFdata[BAFdata[,1]==chr,c(4)]

    #Now take the positions of segment boundaries
    starts         = c(1, which(BAFsegm[-1] != BAFsegm[-length(BAFsegm)])+1)
    ends           = c(which(BAFsegm[-1] != BAFsegm[-length(BAFsegm)]), length(BAFsegm))
    segment.bounds = cbind(starts, ends)

    #Collect those results!
    esti.baf  = NULL
    BAFphased = NULL

    #Run through the segments and model them
    for(seg in 1:length(starts)) {

      #Take the BAF in this segment
      BAF.seg = BAF[starts[seg]:ends[seg]]

      #What is the probability that this is actually just normal?
      ks.normal.p.value = ks.test(BAF.seg, pnorm, mean = 0.5, sd = array.sd)$p.value

      #Only bother modelling if it greater than 10 SNPs, otherwise just take the mean
      if(ks.normal.p.value < 0.05 | !use.ks.gate) {

        #Do the grid search across the standard deviations and the full list of means
        grid.list = lapply(sigmai, function(sigmai) {

          #Calculate mu that produces best loglik
          logliks = lapply(bafsearch, function(mui) {

            n      = length(BAF.seg)
            k      = 2
            x      = BAF.seg
            mu     = c(0.5 - mui, 0.5 + mui)
            sigma  = c(sigmai, sigmai)
            lambda = c(0.5, 0.5)

            #Calculate the logliks for this mui
            res = .C("normpost", as.integer(n), as.integer(k),
                     as.double(x), as.double(mu), as.double(sigma),
                     as.double(lambda), res2 = double(n * k), double(3 * k),
                     post = double(n * k), loglik = double(1),
                     PACKAGE = "mixtools")

            return(res$loglik)

          })

          #Make the logliks a vector
          logliks        = unlist(logliks)

          return(logliks)

        })

        #Grid list
        grid.search = do.call(rbind, grid.list)

        #Name matrix
        rownames(grid.search) = sigmai
        colnames(grid.search) = bafsearch

        #What is the solution?
        maxin   = which(grid.search==grid.search[which.max(grid.search)],
                        arr.ind = TRUE)

        #Define this function outside
        cellfun = function(j, i, x, y, width, height, fill) {
          if(i == maxin[1] & j == maxin[2]) {
            grid.rect(x = x, y = y,
                      width = width,
                      height = height,
                      gp = gpar(col = "black", fill = NA))}}

        #Make the heatmap
        if(plot.lit.plot) {
          png(paste0(samplename,"_lit_plot_chromosome_",chr,"_segment_",seg,".png"), width = 9, height = 9, units = 'in', res = 500)
          print(Heatmap(grid.search,
                        name = "Loglikelihood",
                        column_title = paste0("Chromosome_",chr,"_segment_",seg),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        col = colorRamp2(c(0, max(grid.search)/2, max(grid.search)*0.999, max(grid.search)), c("blue", "red", "orange", "white")),
                        cell_fun = cellfun))
          dev.off()
        }

        #Define the limits of the local search for the standard deviation
        local.sigmai   = na.omit(sigmai[(maxin[1]-1):(maxin[1]+1)])

        #We will investigate 11 different values
        local.sigmai   = seq(local.sigmai[1],
                             local.sigmai[length(local.sigmai)],
                             length.out = 11)

        #Define the limits of the local search for the standard deviation
        local.bafsearch = na.omit(bafsearch[(maxin[2]-1):(maxin[2]+1)])

        #We will investigate 11 different values
        local.bafsearch = seq(local.bafsearch[1],
                              local.bafsearch[length(local.bafsearch)],
                              length.out = 21)

        #Do the grid search across the standard deviations and the full list of means
        local.grid.list = lapply(local.sigmai, function(sigmai) {

          #Calculate mu that produces best loglik
          logliks = lapply(local.bafsearch, function(mui) {

            n      = length(BAF.seg)
            k      = 2
            x      = BAF.seg
            mu     = c(0.5 - mui, 0.5 + mui)
            sigma  = c(sigmai, sigmai)
            lambda = c(0.5, 0.5)

            #Calculate the logliks for this mui
            res = .C("normpost", as.integer(n), as.integer(k),
                     as.double(x), as.double(mu), as.double(sigma),
                     as.double(lambda), res2 = double(n * k), double(3 * k),
                     post = double(n * k), loglik = double(1),
                     PACKAGE = "mixtools")

            return(res$loglik)

          })

          #Make the logliks a vector
          logliks        = unlist(logliks)

          return(logliks)

        })

        #Grid list
        local.grid.search = do.call(rbind, local.grid.list)

        #Name matrix
        rownames(local.grid.search) = local.sigmai
        colnames(local.grid.search) = round(local.bafsearch, digits = 4)

        #What is the solution?
        local.maxin   = which(local.grid.search==local.grid.search[which.max(local.grid.search)],
                              arr.ind = TRUE)

        #Define this function outside
        cellfun = function(j, i, x, y, width, height, fill) {
          if(i == local.maxin[1] & j == local.maxin[2]) {
            grid.rect(x = x, y = y,
                      width = width,
                      height = height,
                      gp = gpar(col = "black", fill = NA))}}

        #Make the heatmap
        if(plot.star.plot) {
          png(paste0(samplename,"_starburst_plot_chromosome_",chr,"_segment_",seg,".png"), width = 9, height = 9, units = 'in', res = 500)
          print(Heatmap(local.grid.search,
                        name = "Loglikelihood",
                        column_title = paste0("Chromosome_",chr,"_segment_",seg),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        col = colorRamp2(c(0, max(local.grid.search)/2, max(local.grid.search)*0.9999, max(local.grid.search)), c("blue", "red", "orange", "white")),
                        cell_fun = cellfun))
          dev.off()
        }

        #So what mean did we find? (the pengest mean)
        seg.mu         = c(0.5 - local.bafsearch[local.maxin[2]],
                           0.5 + local.bafsearch[local.maxin[2]])

        #Now let's create a list called mixmdl for this solution
        mixmdl         = list()

        #Record the mixmdl elements required (x, mu, sigma)
        mixmdl$x       = BAF.seg
        mixmdl$mu      = seg.mu
        mixmdl$sigma   = rep(local.sigmai[local.maxin[1]], times=2)

        #Variable for making a unique segment seed
        seed.variables = c(seed, seg, length(BAF.seg))

        #Save the segment seed
        seg.seed = NULL

        #Make a seed unique this sample segment that will make the segment unique but consistent
        for(svar in seed.variables) {

          #Set the seed and make a number
          set.seed(svar)
          seg.seed = c(seg.seed,
                       sample(0:9,
                              size = 2))

        }

        #Scramble it uniquely for the segment
        set.seed(seg)
        seg.seed = paste0(seg.seed[sample(1:length(seg.seed))],
                          collapse = "")

        #Transform the BAF using a probability function generated from the dists to get BAFphased
        BAFphasedseg = getBAFphased(mixmdl,
                                    seed=seg.seed)

        #Join the BAFs
        BAFphased    = c(BAFphased, BAFphasedseg)

        #Take the phased median as a the segment median, if on the off chance it is still less than 0.5, flip it
        esti.baf[seg] = ifelse(median(BAFphasedseg) < 0.5,
                               1 - median(BAFphasedseg),
                               median(BAFphasedseg))

      } else {

        #Make an announcment
        print(paste0("Called normality in segment ",seg," in chromosome ",chr," of ",samplename))

        #We take a 'pseudo-median', the medium in odd lengths, the near medium in even lengths but always a real value
        p.med.BAF.seg = sort(BAF.seg)[round((length(BAF.seg)+1)/2)]

        #Make the median of the segment just the median BAF
        esti.baf[seg] = max(c(p.med.BAF.seg, abs(p.med.BAF.seg-0.5)+0.5))

        #Use the BAFs for BAFphased
        BAFphased    = c(BAFphased, BAF.seg)

      }

    }

    #Now you have your segments!
    BAFphseg     = rep(esti.baf, times = ends - (starts-1))

    if(plot.transformed) {
      png(filename = paste(samplename,"_BAFphased_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
      plot(pos,
           BAFphased,
           pch=".",
           ylim = c(0,1),
           col="red",
           main=paste0(samplename," transformed BAF chr",chr))
      points(pos,
             BAFphseg,
             pch=".",
             col="blue")
      dev.off()
    }

    BAFoutputchr = cbind(rep(chr, length(BAFphseg)), pos, BAF, BAFphased, BAFphseg)
    BAFoutput = rbind(BAFoutput, BAFoutputchr)

  }

  colnames(BAFoutput) = c("Chromosome","Position","BAF","BAFphased","BAFseg")
  write.table(BAFoutput, outputfile, sep="\t", row.names=F, col.names=T, quote=F)

  #Reset seed
  set.seed(Sys.time())

}
