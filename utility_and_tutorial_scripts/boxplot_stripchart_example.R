# Boxplot and stripchart for pre/post and change in value
# https://www.r-bloggers.com/visualizing-small-scale-paired-data-combining-boxplots-stripcharts-and-confidence-intervals-in-r/

# Set data:
pre <- imp_all_data_completed$Ln_IL6_0
post <- imp_all_data_completed$Ln_IL6_12

# Plot:
png('boxplot_scatterplot_pre_post.png', width = 7.3, height = 5, units = 'in', res = 600)
par(mfrow = c(1, 2))
# First graph
s <- seq(length(pre))
par(bty = "l")
boxplot(pre, post, main = "Raw data", 
        xlab = "Time", ylab = "Measure", 
        names=c("pre", "post"), col = c("lightblue", "lightgreen"))
stripchart(list(pre, post), vertical = T, 
           pch = 16, method = "jitter", 
           cex = 0.5, add = T)
segments(rep(0.95, length(pre))[s], pre[s], rep(2, length(pre))[s],
         post[s], col = 1, lwd = 0.5)
# Second graph
# Confidence intervals:
res <- t.test(post, pre, paired = T, conf.int = T)
# res <- wilcox.test(post, pre, paired = T, conf.int = T)
stripchart(post-pre, vertical = T, pch = 16, 
           method = "jitter", main = "Difference", 
           ylab = "Difference: Post – Pre", xlab = "Median +/- 95% CI")
points(1, res$estimate, col = "red", pch = 16, cex = 2)
arrows(1, res$conf.int[1], 1, res$conf.int[2], col = "red",
       code = 3, lwd = 3, angle = 90)
abline(h = 0, lty = 2) # Zero-effectline
dev.off()