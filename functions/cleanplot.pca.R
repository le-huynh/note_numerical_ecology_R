# cleanplot.pca() function
#
# A function to draw a PCA biplot (scaling 1 or scaling 2) from an object
# of class "rda" containing PCA results, computed with vegan's rda() function.
# A circle of equilibrium contribution is drawn on scaling type 1 biplots.
#
# ARGUMENTS
#
# ##### General parameters
# res.pca          An rda{vegan} object. This object contains a PCA result if a
#                  single matrix was used in the rda(). The function will still
#                  operate with aRDA result but a Warning will be printed.
# ax1, ax2         Canonical axes to be drawn as abscissa and ordinate.
#                  Defaults: 1 and 2.
# scaling          Scaling type: only 1 or 2 are supported. Default: 1.
#
# ##### Items to be plotted
# plot.sites       If TRUE, the sites will be plotted as small circles.
# plot.spe         If TRUE, the species (or other response variables) will be
#                  plotted.
# label.sites      If TRUE, labels are added to the site symbols.
# label.spe        If TRUE, labels are added to the species arrows.
# cex.char1        Character size (for sites and response variables).
#
# ##### Label positions
# ## Positions: 1 = below the point, 2 = left, 3 = above, 4 = right. Default: 4.
# ## Note - Argument pos = NULL centres the label on the position of the object
# ##       (site point, species or environmental variable arrow, centroid) when
# ##       the object is not drawn.
# pos.sites        Position of site labels. 1 to 4, as above. Default: 2.
# pos.spe          Position of species labels. 1 to 4, as above. Default: 4.
#
# ##### Multipliers, selection of species to be plotted
# mult.spe         Multiplier for length of the species arrows. Default: 1.
# select.spe       Vector containing a selection of the species numbers to be
#                  drawn in the biplot, e.g. c(1, 2, 5, 8, 12). Draw all species
#                  if select.spe = NULL (default value). The species that are
#                  well represented in the RDA plot can be identified using
#                  goodness(pca.output.object, display = "species").
#
# ##### Position of the plot in frame, margins
# mar.percent      Factor to expand plot size to accomodate all objects and
#                  labels. Positive values increase the margins around the plot,
#                  negative values reduce them.
# optimum          If TRUE, the longest species arrow is stretched to
#                  a length equal to the distance to the origin of the site
#                  farthest from the origin of the plot of (ax1, ax2). This is
#                  an optimal combined representation of the sites and species.
#                  The lengths of the species arrows can be further modified
#                  using the arguments mult.spe.
# move.origin      Move plot origin right-left and up-down. 
#                  Default: move.origin = c(0,0).
#                  Ex. move.origin = c(-1, 0.5) moves origin by 1 unit left and
#                  0.5 unit up.
#
# ##### Varia
# silent           If FALSE, intermediate computation steps are printed.
#                  Default: TRUE. Interpretation: examine the code lines
#                  controlled by argument "silent".
#
# Reference
# Legendre, P. & L. Legendre. 2012. Numerical ecology, 3rd English edition.
#    Elsevier Science BV, Amsterdam.
#
# # Example 1 - Doubs river fish data provided with the NEwR book
# # Data also available in the R package ade4
#
# # Load the Doubs.RData file from the NEwR (2018) book material.
# load("Doubs.RData")
#
# # The fish community data at 30 sites are in object "spe" (30 sites x 27
# # species). Site 8 may be removed before PCA, where no fish had been caught.
# # We did not do it here.
# ### Not done ###   spe.29 <- spe[-8, ]
# # Examine the position of site 8 in the biplots produced by the example code
# # lines below.
#
# # Chord-transform the fish data
# library(vegan)
# spe.norm <- decostand(spe, "normalize")
# pca.out <- rda(spe.norm)
#
# # Scaling 1
# source("cleanplot.pca.R")
# dev.new(width = 12, height = 7, noRStudioGD = TRUE)
# par(mfrow = c(1, 2))
# cleanplot.pca(pca.out, scaling = 1, optimum = FALSE)
# cleanplot.pca(pca.out, scaling = 1, optimum = TRUE)
#
# # With optimum = TRUE (right-hand graph), the length of the longest species
# # arrow is equal to the distance between the origin and the site farthest
# # from the origin. Minor disadvantage: the radius of the equilibrium
# # contribution circle is not equal to sqrt(2/27) = 0.2721655 in that plot.
# # It has this value in the left-hand graph.
#
# # Argument silent = FALSE provides additional information in the R console, 
# # including the new circle radius.
# cleanplot.pca(pca.out,
#               scaling = 1,
#               optimum = FALSE,
#               silent = FALSE)
# cleanplot.pca(pca.out,
#               scaling = 1,
#               optimum = TRUE,
#               silent = FALSE)
#
# # Compare scaling 1 and scaling 2 biplots
# dev.new(width = 12, height = 7, noRStudioGD = TRUE)
# par(mfrow = c(1, 2))
# cleanplot.pca(pca.out, scaling = 1, optimum = TRUE)
# cleanplot.pca(pca.out, scaling = 2, optimum = TRUE)
#
#
# # Example 2, part 1
# # Cajo ter Braak's dune meadow vegetation data (20 sites x 30 species)
#
# library(vegan)
# data(dune)
# dune.hel <- decostand(dune, "hellinger")
# pca.dune.hel <- rda(dune.hel)
#
# # Compare scaling 1 and scaling 2 biplots
# dev.new(width = 12, height = 7, noRStudioGD = TRUE)
# par(mfrow = c(1, 2))
# cleanplot.pca(pca.dune.hel, scaling = 1, optimum = TRUE)
# cleanplot.pca(pca.dune.hel, scaling = 2, optimum = TRUE)
#
# # Example 2, part 2
#
# # The species in the previous plots are very crowded
# # A RDA of a community matrix by itself is a PCA of that matrix.
# # A demonstration is found in Legendre & Legendre (2012), p. 637.
# # The goodness function can only be called on RDA or CCA results
# rda.dune.hel <- rda(dune.hel ~ ., data = dune.hel)
# tmp <- goodness(rda.dune.hel)
# ( sp.sel <- which(tmp[, 2] >= 0.4) )   # 14 species selected for plotting
#
# # Scaling 2
# # Make the plots less crowded by only printing the best represented species
# # Right-hand graph: selected species only, and also example of moving the 
# # centroid of the points in the biplot.
# dev.new(width = 12, height = 7, noRStudioGD = TRUE)
# par(mfrow = c(1, 2))
# cleanplot.pca(pca.dune.hel, scaling = 2)  # All species
# cleanplot.pca(
#   pca.dune.hel,
#   scaling = 2,
#   select.spe = sp.sel,
#   move.origin = c(-0.3, 0)
# )
#
#
# License: GPL-2
# Authors: Francois Gillet, Daniel Borcard & Pierre Legendre
#          2016–2020

'cleanplot.pca' <-
  function(res.pca,
           ax1 = 1,
           ax2 = 2,
           scaling = 1,
           plot.sites = TRUE,
           plot.spe = TRUE,
           label.sites = TRUE,
           label.spe = TRUE,
           cex.char1 = 0.7,
           pos.sites = 2,
           pos.spe = 4,
           mult.spe = 1,
           select.spe = NULL,
           mar.percent = 0.1,
           optimum = TRUE,
           move.origin = c(0, 0),
           silent = TRUE) {
    
    ### Internal functions
    'stretch' <-
      function(sites, mat, ax1, ax2, n, silent = silent) {
        # Compute stretching factor for the species arrows
        # First, compute the longest distance to centroid for the sites
        tmp1 <- rbind(c(0, 0), sites[, c(ax1, ax2)])
        D <- dist(tmp1)
        target <- max(D[1:n])
        # Then, compute the longest distance to centroid for the species arrows
        if("matrix" %in% class(mat)) {
#        if (is.matrix(mat)) {
          p <- nrow(mat)   # Number of species to be drawn
          tmp2 <- rbind(c(0, 0), mat[, c(ax1, ax2)])
          D <- dist(tmp2)
          longest <- max(D[1:p])
        } else {
          tmp2 <- rbind(c(0, 0), mat[c(ax1, ax2)])
          longest <- dist(tmp2)
          # print(tmp2)
        }  # If a single row left in 'mat'
        #
        if (!silent)
          cat("target =",
              target,
              " longest =",
              longest,
              " fact =",
              target / longest,
              "\n")
        fact <- target / longest
      }
    
    'larger.plot' <-
      function(sit.sc,
               spe.sc,
               percent,
               move.origin,
               ax1,
               ax2) {
        # Internal function to expand plot limits (adapted from code by Pierre
        # Legendre)
        mat <- rbind(sit.sc, spe.sc)
        range.mat <- apply(mat, 2, range)
        rownames(range.mat) <- c("Min", "Max")
        z <- apply(range.mat, 2, function(x)
          x[2] - x[1])
        range.mat[1,] <- range.mat[1,] - z * percent
        range.mat[2,] <- range.mat[2,] + z * percent
        if (move.origin[1] != 0)
          range.mat[, ax1] <- range.mat[, ax1] - move.origin[1]
        if (move.origin[2] != 0)
          range.mat[, ax2] <- range.mat[, ax2] - move.origin[2]
        range.mat
      }
    
    "pcacircle" <-
      function (pca, mult.spe, fact.spe, silent = silent) {
        # This function draws a circle of equilibrium contribution on a PCA plot
        # generated from the result file of a vegan rda() analysis.
        eigenv <- pca$CA$eig
        p <- length(eigenv)
        n <- nrow(pca$CA$u)
        tot <- sum(eigenv)
        radius <- (2 / p) ^ 0.5 * mult.spe * fact.spe
        symbols(
          0,
          0,
          circles = radius,
          inches = FALSE,
          add = TRUE,
          fg = 2
        )
        if (!silent) {
          cat(
            "\nSpecies arrows and the radius of the equilibrium circle are stretched ",
            "by a factor of",
            mult.spe * fact.spe
          )
          cat(
            "\nThe radius of the equilibrium circle is thus",
            (2 / p) ^ 0.5,
            "*",
            mult.spe,
            "*",
            fact.spe,
            "=",
            radius,
            "\n"
          )
        }
      }
    ### End internal functions
    
    if (!class(res.pca)[1] == "rda")
      stop("The input file is not a vegan output object of class 'rda'",
           call. = FALSE)
    if (!(is.null(res.pca$CCA)))
      stop(
        "The input file contains an RDA, not a PCA result. ",
        "Use function triplot.rda from the NEwR (2018) book to produce an RDA triplot."
      )
    if (scaling != 1 &
        scaling != 2)
      stop("Function only available for scaling 1 or 2", call. = FALSE)
    
    k <- length(res.pca$CA$eig)         # n. of PCA eigenvalues
    n.sp <- length(res.pca$colsum)      # n. of species
    ahead <- 0.05   # Length of arrow heads
    aangle <- 30    # Angle of arrow heads
    # 'vec' will contain the selection of species to be drawn
    if (is.null(select.spe)) {
      vec <- 1:n.sp
    } else {
      vec <- select.spe
    }
    
    # Scaling 1: the species scores have norms of 1
    # Scaling 1: the site scores are scaled to variances = can.eigenvalues
    # Scaling 2: the species scores have norms of sqrt(can.eigenvalues)
    # Scaling 2: the site scores are scaled to variances of 1
    
    # This version reconstructs and uses the original RDA output of L&L 2012,
    # Section 11.1.3
    
    Tot.var = res.pca$tot.chi         # Total variance in response data Y
    eig.val = res.pca$CA$eig          # Eigenvalues of Y-hat
    Lambda = diag(eig.val)            # Diagonal matrix of eigenvalues
    eig.val.rel = eig.val / Tot.var   # Relative eigenvalues of Y-hat
    Diag = diag(sqrt(eig.val.rel))    # Diagonal matrix of sqrt(relative
                                      # eigenvalues)
    U.sc1 = res.pca$CA$v              # Species scores, scaling=1
    U.sc2 = U.sc1 %*% Lambda ^ (0.5)  # Species scores, scaling=2
    n = nrow(res.pca$CA$u)            # Number of observations
    Z.sc2 = res.pca$CA$u * sqrt(n - 1)# Site scores, scaling=2
    Z.sc1 = Z.sc2 %*% Lambda ^ (0.5)  # Site scores, scaling=1
    
    if (is.null(select.spe)) {
      vec <- 1:n.sp
    } else {
      vec <- select.spe
    }
    
    if (scaling == 1) {
      sit.sc <- Z.sc1
      spe.sc <- U.sc1[vec,]
    } else {
      # For scaling=2
      sit.sc <- Z.sc2
      spe.sc <- U.sc2[vec,]
    }
    if (is.null(rownames(sit.sc)))
      rownames(sit.sc) <- paste("Site", 1:n, sep = "")
    if (is.null(rownames(spe.sc)))
      rownames(spe.sc) <- paste("Sp", 1:n.sp, sep = "")
    
    fact.spe <- 1
    if (optimum) {
      fact.spe <-
        stretch(sit.sc[, 1:k], spe.sc[, 1:k], ax1, ax2, n, silent = silent)
    }
    if (!silent)
      cat("fact.spe =", fact.spe, "\n\n")
    spe.sc <- spe.sc * fact.spe * mult.spe
    
    lim <-
      larger.plot(
        sit.sc[, 1:k],
        spe.sc[, 1:k],
        percent = mar.percent,
        move.origin = move.origin,
        ax1 = ax1,
        ax2 = ax2
      )
    if (!silent)
      print(lim)
    
    # Draw the main plot
    mat <- rbind(sit.sc[, 1:k], spe.sc[, 1:k])
    plot(
      mat[, c(ax1, ax2)],
      type = "n",
      main = paste("PCA biplot - Scaling", scaling),
      xlim = c(lim[1, ax1], lim[2, ax1]),
      ylim = c(lim[1, ax2], lim[2, ax2]),
      xlab = paste("PCA ", ax1),
      ylab = paste("PCA ", ax2),
      asp = 1
    )
    abline(h = 0, v = 0, col = "grey60")
    
    # Draw the site scores
    if (plot.sites) {
      points(sit.sc[, ax1], sit.sc[, ax2], pch = 20)
      if (label.sites)
        text(
          sit.sc[, ax1],
          sit.sc[, ax2],
          labels = rownames(sit.sc),
          col = "black",
          pos = pos.sites,
          cex = cex.char1
        )
    } else {
      if (label.sites)
        text(
          sit.sc[, ax1],
          sit.sc[, ax2],
          labels = rownames(sit.sc),
          col = "black",
          pos = NULL,
          cex = cex.char1
        )
    }
    
    # Draw the species scores
    if (plot.spe) {
      arrows(
        0,
        0,
        spe.sc[, ax1],
        spe.sc[, ax2],
        length = ahead,
        angle = aangle,
        col = "red"
      )
      if (label.spe)
        text(
          spe.sc[, ax1],
          spe.sc[, ax2],
          labels = rownames(spe.sc),
          col = "red",
          pos = pos.spe,
          cex = cex.char1
        )
    } else {
      if (label.spe)
        text(
          spe.sc[, ax1],
          spe.sc[, ax2],
          labels = rownames(spe.sc),
          col = "red",
          pos = NULL,
          cex = cex.char1
        )
    }
    
    # If scaling = 1 draw circle of equilibrium contribution
    if (scaling == 1) {
      pcacircle(
        res.pca,
        mult.spe = mult.spe,
        fact.spe = fact.spe,
        silent = silent
      )
    }
  }
