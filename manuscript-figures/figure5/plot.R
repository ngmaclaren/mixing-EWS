library(sfsmisc)
library(latex2exp)
load("data.rda")
palette("Tableau 10")

pdf("demonstration.pdf", width = 8)
par(mar = c(4, 5, 0.5, 4)+0.5)
plot(NULL, xlim = range(us), ylim = range(c(avgvar_opt, avgvar_rand)),
     xlab = "u", ylab = TeX(r"(Average variance, $\hat{V}_S$)"),
     yaxt = "n", cex.axis = ticksize, cex.lab = labelsize)
eaxis(2, axTicks(2)[c(TRUE, FALSE)], cex.axis = ticksize, at.small = FALSE, sub10 = "10", las = 0)
lines(us[all_ls], avgvar_rand, lty = 1, lwd = 3, col = 2)
lines(us[all_ls], avgvar_opt, lty = 1, lwd = 3, col = 1)
par(new = TRUE)
matplot(us, plotX, type = "l", lty = 1, lwd = .5, col = "black",
        xlab = "", ylab = "", axes = FALSE)
axis(4, pretty(range(as.numeric(plotX))), cex.axis = ticksize)
mtext(TeX(r"($x_i^*$)"), 4, 3, cex = labelsize)
legend("topleft", bty = "n", lwd = 3, col = c(1, 2, "black"), cex = .75*labelsize,
       legend = c(TeX(r"($\hat{V}_S$, optimized node set)"),
                  TeX(r"($\hat{V}_S$, random node set)"), 
                  TeX(r"($x_i^*$)")))
dev.off()
