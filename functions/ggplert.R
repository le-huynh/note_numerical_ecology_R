#' Parametric tests for multiple comparisons with ggplot
#'
#' Unpaired ANOVA or paired t-test with post-hoc
#' multiple comparisons and boxplots with letters using ggplot2.
#' 
#' Function to perform (unpaired) ANOVA or paired t-test with post-hoc
#' multiple comparisons and boxplots with letters (using agricolae::LSD.test).
#' ggplert() is based on an adaptation to ggplot2 of boxplert() by Chinese 
#' students of Jiangshan Lai (http://blog.sciencenet.cn/blog-1835014-1151822.html).
#' This function requires ggplot2 and ggpubr libraries. 
#' Mean values are added as big points.
#' 
#' @author François Gillet, 2021-08-22
#' @param Y the numeric response variable
#' @param X the explanatory factor (groups, treatments...) or a qualitative 
#' variable that can be converted to a factor
#' @param paired if TRUE and if 2 groups only, then paired t-test is performed
#' (objects are supposed to be sorted twice the same way in the data frame)
#' @param p.adj correction of p-values for multiple comparisons after 
#' ANOVA test (default = "holm", else "none", "bonferroni",...)
#' @param varwidth if TRUE, adapt boxplot width to relative group size
#' @param notch if TRUE, add notches to help multiple comparisons
#' @param bcol fill colour(s) of boxplots; if bcol = "" then boxplots are 
#' filled according to a gradient of heat colors based on means
#' @param xlab x-axis title (name of the explanatory factor X)
#' @param ylab y-axis title (name of the numeric response variable Y)
#' @return A summary table ($comparison) is printed and boxplots ($plot) are 
#' drawn with result of the test ($p.value) and, if it is significant, of 
#' post-hoc tests as letters (in decreasing order); the method applied to  
#' adjustp-values is also returned ($p.adjust)
#' @export
#' @examples
#' data(sweetpotato)
#' 
#' # Check ANOVA assumptions
#' shapiro.test(resid(aov(sweetpotato$yield ~ sweetpotato$virus)))
#' bartlett.test(sweetpotato$yield, sweetpotato$virus)
#' 
#' # Since assumptions are met, perform the parametric tests:
#' ggplert(
#'   sweetpotato$yield,
#'   sweetpotato$virus,
#'   xlab = "Virus",
#'   ylab = "Yield",
#'   bcol = "",
#'   p.adj = "holm")

ggplert <- function(Y,
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

  # ANOVA test
  if (paired == TRUE & length(aa) == 2) {
    res <- t.test(Y ~ X, paired = TRUE)
    pp <- res$p.value
    tt1$group <- rep("a", 2)
    if (pp <= 0.05) 
      tt1$group[which(tt1$mean == min(tt1$mean))] <- "b"
  }
  else {
    model <- aov(Y ~ X)
    pp <- anova(model)$Pr[1]
    comp <- agricolae::LSD.test(model,
                                "X",
                                alpha = 0.05,
                                p.adj = p.adj,
                                group = TRUE)
    gror <- comp$groups[aa,]
    tt1$group <- gror$groups
  }
  
  if (is.null(bcol)) {
    bcol <- heat.colors(nrow(tt1), rev = TRUE)[rank(tt1$mean)]
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
