                                        # This is τ against d
                                        # v2
plot_result <- function(result, simrun, filename = NULL) {
    library(latex2exp)
    palette("Tableau 10")
    S_d <- result$d
    S_ews <- result$tau
    S_stdev <- result$stdev
    largesd <- result$largesd
    
    neach <- lengths(S_d)
    ptcol <- 1:length(S_d)
    ptcols <- unlist(mapply(function(n1, n2) rep(n1, n2), ptcol, neach))

    ht <- 8
    wd <- 12
    labelsize <- 2.5
    ticksize <- 2
    if(save_plots) {
        pdf(filename, height = ht, width = wd, family = "sans")
    } else dev.new(height = ht, width = wd)
    par(mfrow = c(2, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
    xlim <- range(unlist(S_ews))
    ylim <- range(unlist(S_d))
    for(i in seq_along(S_d)) {
        x <- S_ews[[i]]
        y <- S_d[[i]]

        maxd <- which.max(y)
        max_y <- y[maxd]
        max_x <- x[maxd]
        
        plot(
            x, y,
            xlab = TeX(r"($\tau$)", italic = TRUE), ylab = TeX(r"($d$)", italic = TRUE),
            xlim = xlim, ylim = ylim,
            pch = 1, col = 1, cex = 1, lwd = 2,
            cex.axis = ticksize, cex.lab = labelsize#, font.lab = 3
        )
        points(
            x = largesd$tau[[i]],
            y = largesd$d[[i]],
            pch = 4, cex = 2, lwd = 2, col = 2
        )
                                        # Uncomment to compare node selected by greedy algorithm
        ## if(i < length(S_d)) {
        ##     add_greedy_point(i, simrun)
        ## } else {
        ##     add_greedy_point(N, simrun)
        ## }

        abline(v = max_x, col = "gray", lwd = .9, lty = 2)
        abline(h = max_y, col = "gray", lwd = .9, lty = 2)

        if(i == 1) {
            legend("topright", bty = "n", legend = c("Sampled node sets", "Largest SD node set"),
                   col = 1:2, pch = c(1, 4), pt.lwd = 2, cex = 1.75)
        }
        mtext(paste0("(", letters[i], ")"), line = -2.5, adj = 0.02, font = 2, cex = 2)
    }
    if(save_plots) dev.off()
}
