library(sfsmisc)
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
load("data.rda")
palette("Tableau 10")

pdf("./demonstrate-stoppingcriteria.pdf", width = 8)
par(mar = c(4, 4, 0.75, 4)+0.5)
plot(
    1:N, dvals, type = "l", col = 1, lwd = 3, lty = 1,
    xlab = "n", ylab = "d", cex.axis = ticksize, cex.lab = labelsize
)
abline(v = markns, col = "gray", lty = 3:2, lwd = .9*3)
par(new = TRUE)
plot(
    1:N, taus, type = "l", col = 2, lwd = 3, lty = 1,
    xlab = "", ylab = "", axes = FALSE
)
axis(4, pretty(range(taus)), cex.axis = ticksize)
mtext(expression(tau), 4, 3, cex = labelsize)
legend("bottomright", bty = "n", lwd = 3, col = 1:2, cex = labelsize,
       legend = c("d", expression(tau)))
dev.off()
