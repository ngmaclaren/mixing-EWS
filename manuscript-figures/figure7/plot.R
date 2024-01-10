library(latex2exp)
load("data.rda")

pdf("summaryfig-ba-doublewell.pdf", height = 4, width = 12)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
labelsize <- 2.5
ticksize <- 2
palette("Tableau 10")
for(i in seq_along(dl)) {

    props <- dl[[i]]$props
    dists <- dl[[i]]$dists
    
                                        # this is a panel
    plot(
        NULL, xlim = c(.5, 5.5), ylim = c(0, 2),
        xlab = TeX(r"($n$)", italic = TRUE), ylab = TeX(r"($p_1$, $p_2$)", italic = TRUE),
        cex.axis = ticksize, cex.lab = labelsize
    )
    abline(h = 1, lwd = 0.5, lty = 1, col = "black")

    matpoints(
        x = matrix(jitter(rep(1:5 - 0.15, each = nrow(props)), amount = .05), nrow = nrow(props)),
        y = 2*props,
        pch = 0, col = 1, cex = 1
    )

    points(
        x = 1:5 - 0.15, y = 2*colMeans(props),
        pch = 0, col = 1, cex = 5, lwd = 2.5
    )

    matpoints(
        x = matrix(jitter(rep(1:5 + 0.15, each = nrow(dists)), amount = .05), nrow = nrow(dists)),
        y = dists,
        pch = 1, col = 2, cex = 1
    )

    points(
        x = 1:5 + 0.15, y = colMeans(dists),
        pch = 1, col = 2, cex = 5, lwd = 2.5
    )
    
    mtext(paste0("(", letters[i], ")"), line = -2.5, adj = 0.02, font = 2, cex = 2)

    if(i == 1) {
        legend("topright", bty = "n", col = 1:2, pch = c(0, 1), pt.cex = 2, pt.lwd = 2,
               legend = c(TeX(r"($p_1$)", italic = TRUE), TeX(r"($p_2$)", italic = TRUE)), cex = 1.75)
    }
}
dev.off()
