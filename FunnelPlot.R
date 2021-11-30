library(dplyr)
library(metafor)
library(ggplot2)
GSVA_results <- read.csv("/Users/sheerin.d/OneDrive\ -\ wehi.edu.au/2020/COVID-19/meta_analysis/og_GSVA_res.csv", row.names = 1)
# can ignore the getDOR function if we use the AUC vs. SE(AUC)
getDOR <- function(rocObject) {
  optimal <- suppressWarnings(pROC::coords(rocObject, "best", input = "threshold",
                                           ret = c("threshold", "specificity", "sensiotivity")))
  if (rocObject$direction == ">") {
    predicted <- ifelse(rocObject$original.predictor > optimal[1, 1],
                        rocObject$levels[1], rocObject$levels[2])
  } else {
    predicted <- ifelse(rocObject$original.predictor < optimal[1, 1],
                        rocObject$levels[1], rocObject$levels[2])
  }
  cfmat <- caret::confusionMatrix(as.factor(rocObject$original.response),
                                  as.factor(predicted))
  logDOR <- log(cfmat$table[1,1] * cfmat$table[2,2] / (cfmat$table[1,2] * cfmat$table[2,1]))
  SE_logDOR <- sqrt(1 / cfmat$table[1,1] + 1 / cfmat$table[1,2] +
                      1/ cfmat$table[2,1] + 1 / cfmat$table[2,2])
  return(data.frame(logDOR = logDOR, SE_logDOR = SE_logDOR))
}
# Based on biomarkers name
biomarkerName <- unique(GSVA_results$variable)
GSVA_results_biomarker <- lapply(1:length(biomarkerName), function(i) {
      x <- biomarkerName[i]
      GSVA_group_biomarker <- GSVA_results %>%
        dplyr::filter(variable == x)
      score <- GSVA_group_biomarker %>%
        dplyr::filter(TBStatus == "1. London Latent") %>%
        dplyr::select(value) %>%
        unlist(use.names = FALSE)
      outcome <- rep("1. London Latent", length(score))
      TBStatusList <- unique(GSVA_results$TBStatus)
      TBStatusListSub <- TBStatusList[TBStatusList != "1. London Latent"]
      result <- lapply(TBStatusListSub, function(y) {
          scoreAlt <- GSVA_group_biomarker %>%
            dplyr::filter(TBStatus == y) %>%
            dplyr::select(value) %>%
            unlist(use.names = FALSE)
          scoreFinal <- c(score, scoreAlt)
          outFinal <- c(outcome, rep(y, length(scoreAlt)))
          rocObject <- pROC::roc(outFinal, scoreFinal)
          dor <- getDOR(rocObject)
          data.frame(AUC = max(as.numeric(rocObject$auc), 1 - as.numeric(rocObject$auc)),
                     SE_AUC = sqrt(pROC::var(rocObject)),
                     logDOR = dor$logDOR, SE_logDOR = dor$SE_logDOR,
                     Signature = biomarkerName[i], Comparison = y)
      })
      do.call(rbind, result)
  })
GSVA_biomarker_all <- do.call(rbind, GSVA_results_biomarker)

# Based on signature
dfCI_biomarker <- lapply(1:length(GSVA_results_biomarker), function (i, effect, precision) {
  data <- GSVA_results_biomarker[[i]]
  effect <- data[, effect]
  estimate <- mean(effect[is.finite(effect)])

  precision <- data[, precision]
  se.seq <- seq(0, max(precision[is.finite(precision)]), 0.001)
  ll95 <- estimate - (1.96 * se.seq)
  ul95 <- estimate + (1.96 * se.seq)
  data.frame(ll95, ul95, se.seq, estimate, Signature = biomarkerName[i])
}, "AUC", "SE_AUC")

dfCI <- do.call(rbind, dfCI_biomarker)
scaleFUN <- function(x) sprintf("%.2f", x)
p <- ggplot(aes(x = AUC, y = SE_AUC), data = GSVA_biomarker_all) +
  geom_point(shape = 20) +
  xlab("AUC") + ylab("SE(AUC)") + facet_wrap(. ~ Signature) +
  geom_line(aes(x = ll95, y = se.seq), linetype = "dashed", data = dfCI) +
  geom_line(aes(x = ul95, y = se.seq), linetype = "dashed", data = dfCI) +
  geom_segment(aes(y = min(se.seq), x = estimate,
                   yend = max(se.seq), xend = estimate), data = dfCI) +
  scale_y_reverse() + scale_x_continuous(label = scaleFUN) + theme_bw()
p
pval_list <- lapply(1:length(GSVA_results_biomarker), function(i, estimate, precision) {
  estimate <- GSVA_results_biomarker[[i]][, estimate]
  estimate <- estimate[is.finite(estimate)]
  precision <- GSVA_results_biomarker[[i]][, precision]
  precision <- precision[is.finite(precision)]
  eggerTest <- metafor::regtest(estimate, sei = precision)

  data.frame(label = paste("p =", round(eggerTest$pval, 3)),
             Signature = biomarkerName[i])
}, "AUC", "SE_AUC")

pval_data <- do.call(rbind, pval_list)

p + geom_text(data = pval_data,
              mapping = aes(x = 0.5, y = 0, label = label), hjust = 0.1)

#### Based on comparison category (worse, strong correlation bet. AUC and SE(AUC))
# TBStatusName <- unique(GSVA_results$TBStatus)[-2]
# GSVA_results_TBStatus <- lapply(TBStatusName, function(x) {
#   indexTBStatus <- which(GSVA_biomarker_all$Comparison == x)
#   GSVA_biomarker_all[indexTBStatus, ]
# })
#
# dfCI_TBStatus <- lapply(1:length(GSVA_results_TBStatus), function (i) {
#   data <- GSVA_results_TBStatus[[i]]
#   estimate <- mean(data$AUC)
#   se.seq <- seq(0, max(data$SE_AUC), 0.001)
#   ll95 <- estimate - (1.96 * se.seq)
#   ul95 <- estimate + (1.96 * se.seq)
#   data.frame(ll95, ul95, se.seq, estimate, Comparison = TBStatusName[i])
# })
# dfCI2 <- do.call(rbind, dfCI_TBStatus)
#
# p2 <- ggplot(aes(x = AUC, y = SE_AUC), data = GSVA_biomarker_all) +
#   geom_point(shape = 20) +
#   xlab("AUC") + ylab("SE(AUC)") + facet_wrap(. ~ Comparison) +
#   geom_line(aes(x = ll95, y = se.seq), linetype = "dashed", data = dfCI2) +
#   geom_line(aes(x = ul95, y = se.seq), linetype = "dashed", data = dfCI2) +
#   geom_segment(aes(y = min(se.seq), x = estimate,
#                    yend = max(se.seq), xend = estimate), data = dfCI2) +
#   scale_y_reverse() + scale_x_continuous(label = scaleFUN) + theme_bw()
#
# pval_list2 <- lapply(1:length(GSVA_results_TBStatus), function(i) {
#   eggerTest <- metafor::regtest(GSVA_results_TBStatus[[i]]$AUC,
#                                 sei = GSVA_results_TBStatus[[i]]$SE_AUC)
#   data.frame(label = paste("p =", round(eggerTest$pval, 3)),
#              Comparison = TBStatusName[i])
# })
# pval_data2 <- do.call(rbind, pval_list2)
#
# p2 + geom_text(data = pval_data2,
#               mapping = aes(x = 0.5, y = 0, label = label), hjust = 0.1)
#
