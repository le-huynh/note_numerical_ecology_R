#' Non-parametric tests for multiple comparisons with ggplot
#'
#' Kruskal-Wallis or paired Wilcoxon test with post-hoc
#' multiple comparisons and boxplots with letters using ggplot2.
#' 
#' Function to perform (unpaired) Kruskal-Wallis (using agricolae::kruskal) 
#' or paired Wilcoxon test with post-hoc multiple comparisons and boxplots 
#' with letters.
#' ggplerk() is an adaptation of boxplerk() to ggplot2. This function requires 
#' ggplot2 and ggpubr libraries.
#' Mean values are added as big points.
#' 
#' @author François Gillet, 2021-08-22
#' @param Y the numeric response variable
#' @param X the explanatory factor (groups, treatments...) or a qualitative 
#' variable that can be converted to a factor
#' @param paired if TRUE and if 2 groups only, then Wilcoxon paired test is 
#' performed (objects are supposed to be ordered twice the same way in the 
#' data frame)
#' @param p.adj correction of p-values for multiple comparisons after 
#' Kruskal-Wallis test (default = "holm", else "none", "bonferroni",...)
#' @param varwidth if TRUE, adapt boxplot width to relative group size
#' @param notch if TRUE, add notches to help multiple comparisons
#' @param bcol fill colour(s) of boxplots; if bcol = "" then boxplots are 
#' filled according to a gradient of heat colors based on rank sums
#' @param xlab x-axis title (name of the explanatory factor X)
#' @param ylab y-axis title (name of the numeric response variable Y)
#' @return A summary table ($comparison) is printed and boxplots ($plot) are 
#' drawn with result of the test ($p.value) and, if it is significant, of  
#' post-hoc tests as letters (in decreasing order); the method applied to  
#' adjust p-values is also returned ($p.adjust)
#' @export
#' @examples
#' library(stats)
#' data(InsectSprays)
#'   
#' # Check ANOVA assumptions
#' shapiro.test(resid(aov(InsectSprays$count ~ InsectSprays$spray)))
#' bartlett.test(InsectSprays$count, InsectSprays$spray)
#' 
#' # Since assumptions are NOT met, perform the non-parametric tests:
#' ggplerk(
#'   InsectSprays$count,
#'   InsectSprays$spray,
#'   p.adj = "holm",
#'   ylab = "Count",
#'   xlab = "Spray")

ggplerk <- function(Y,
                    X,
                    paired = FALSE,
                    p.adj = "holm",
                    varwidth = TRUE,
                    notch = FALSE,
                    bcol = NULL,
                    xlab = "X factor",
                    ylab = "Y value") {
  
  if (!is.factor(X)) {
    X <- as.factor(X)
  }
  aa <- levels(X)
  data <- data.frame(X, Y)
  
  tt1 <- matrix(nrow = length(aa), ncol = 7)
  for (i in 1:length(aa)) {
    temp <- Y[X == aa[i]]
    tt1[i, 1] <- mean(temp, na.rm = TRUE)
    tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
    tt1[i, 3] <- sd(temp, na.rm = TRUE)
    tt1[i, 4] <- min(temp, na.rm = TRUE)
    tt1[i, 5] <- max(temp, na.rm = TRUE)
    tt1[i, 6] <- median(temp, na.rm = TRUE)
    tt1[i, 7] <- length(temp)
  }
  tt1 <- as.data.frame(tt1)
  row.names(tt1) <- aa
  colnames(tt1) <- c("mean", "se", "sd", "min", "max", "median", "n")
  
  # Kruskal-Wallis test
  comp <- agricolae::kruskal(Y, X, p.adj = p.adj)
  gror <- comp$groups[aa, ]
  tt1$rank <- gror$Y
  tt1$group <- gror$groups
  
  if (paired == TRUE & length(aa) == 2) {
    res <- wilcox.test(Y ~ X, paired = TRUE)
    pp <- res$p.value
    tt1$group <- rep("a", 2)
    if (pp <= 0.05) 
      tt1$group[which(tt1$rank == min(tt1$rank))] <- "b"
  }
  else {
    pp <- comp$statistics$p.chisq
  }
  
  if (is.null(bcol)) {
    bcol <- heat.colors(nrow(tt1), rev = TRUE)[rank(tt1$rank)]
  }
  
  p <- ggplot(data = data, aes(X, Y, aes(X, Y, fill = aa))) +
    geom_boxplot(fill = bcol,
                 varwidth = varwidth,
                 notch = notch,
                 outlier.colour = 1,
                 outlier.shape = 1,
                 outlier.alpha = 0.5,
                 na.rm = TRUE) +
    labs(x = xlab, y = ylab) +
    theme_pubr()
  
  if (pp <= 0.1) {
    p <- p + stat_summary(
      geom = "text",
      label = tt1$group,
      fun = "max",
      size = 5,
      colour = "steelblue4",
      vjust = -0.6)
  }
  
  if (pp < 0.001)
    st <- expression(italic(P) < 0.001)
  else {
    pp2 <- round(pp, 4)
    st <- bquote(italic(P) == .(pp2))
  }
  p <- p + labs(subtitle = st)
  
  p <- p + stat_summary(
    geom = "point",
    fun = "mean",
    pch = 21,
    size = 4,
    fill = "white",
    alpha = 0.75)
  
  list(comparison = tt1,
       p.value = pp,
       p.adjust = p.adj,
       plot = p)
}
