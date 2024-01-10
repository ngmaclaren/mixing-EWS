library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
library(sfsmisc)
save_plots <- TRUE
palette("Tableau 10")

load("stochastic-forcing.rda")
                                        # Panel A
pdf("./stochastic-forcing-demonstration.pdf", width = 8)
par(mar = c(4, 5, 1, 4)+0.8)
matplot(seq_len(nrow(plotX)), plotX, type = "l", lty = 1, lwd = .5, col = "black",
        xlab = "", ylab = "", axes = FALSE, ylim = c(0.5, 5.5), xaxs = "i", yaxs = "i")
axis(4, pretty(c(0.5, 5.5)),
     cex.axis = ticksize)
mtext(TeX(r"($x_i^*$)"), 4, 3, cex = labelsize)
par(new = TRUE)
plot(NULL, xlim = range(us), ylim = c(0, 0.008),
     xlab = "Sample index, k", ylab = TeX(r"(Average variance, $\hat{V}_S$)"),
     cex.axis = ticksize, cex.lab = labelsize, xaxs = "i", yaxs = "i")
lines(us, c(Xs_var[, tipping_node], rep(NA, 5)), lty = 1, lwd = 3, col = 2)
lines(us, c(avgvar_opt, rep(NA, 5)), lty = 1, lwd = 3, col = 1)
legend("topleft", bty = "n", lwd = 3, col = c(1, 2, "black"), cex = .75*labelsize,
       legend = c(TeX(r"($\hat{V}_S$, optimized node set)"),
                  TeX(r"($\hat{V}_S$, tipping node)"),
                  TeX(r"($x_i^*$)")))
dev.off()

                                        # Panel B
with(
    list(
        result = sr_result,
        simrun = list(Xs, Cs),
        N = ncol(Xs[[1]]),
        filename = "./stochastic-forcing-sixpanel.pdf"
    ), {
        S_d <- result$d
        S_ews <- result$tau
        S_stdev <- result$stdev
        largesd <- result$largesd
        neach <- lengths(S_d)
        ptcol <- 1:length(S_d)
        ptcols <- unlist(mapply(function(n1, n2) rep(n1, n2), ptcol, neach))
        plotNs <- c(1:5, N)
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
                cex.axis = ticksize, cex.lab = labelsize
            )
            abline(v = max_x, col = "gray", lwd = .9, lty = 2)
            abline(h = max_y, col = "gray", lwd = .9, lty = 2)
            mtext(paste("n =", plotNs[i]), line = -2.5, adj = 0.02, font = 1, cex = .75*labelsize)
        }
        if(save_plots) dev.off()
    }
)

load("sudden-stress.rda")
                                        # Panel C
pdf("./sudden-stress-demonstration.pdf", width = 8)
palette("Tableau 10")
par(mar = c(4, 5, 1, 4)+0.8)
matplot(us, plotX, type = "l", lty = 1, lwd = .5, col = "black",
        xlab = "", ylab = "", axes = FALSE, ylim = c(0.5, 6), xaxs = "i", yaxs = "i")
axis(4, pretty(c(0.5, 5.5)),
     cex.axis = ticksize)
mtext(TeX(r"($x_i^*$)"), 4, 3, cex = labelsize)
par(new = TRUE)
plot(NULL, xlim = range(us),
     ylim = c(0, 0.006),
     xlab = "Sample index, k", ylab = TeX(r"(Average variance, $\hat{V}_S$)"),
     cex.axis = ticksize, cex.lab = labelsize, xaxs = "i", yaxs = "i")
lines(us, c(Xs_var[, tipping_node], rep(NA, 10)), lty = 1, lwd = 3, col = 2)
lines(us, c(avgvar_opt, rep(NA, 10)), lty = 1, lwd = 3, col = 1)
legend("topleft", bty = "n", lwd = 3, col = c(1, 2, "black"), cex = .75*labelsize,
       legend = c(TeX(r"($\hat{V}_S$, optimized node set)"),
                  TeX(r"($\hat{V}_S$, tipping node)"),
                  TeX(r"($x_i^*$)")))
dev.off()

                                        # Panel D
with(
    list(
        result = sr_result,
        simrun = list(Xs, Cs),
        N = ncol(Xs[[1]]),
        filename = "./sudden-stress-sixpanel.pdf"
    ), {
        S_d <- result$d
        S_ews <- result$tau
        S_stdev <- result$stdev
        largesd <- result$largesd
        neach <- lengths(S_d)
        ptcol <- 1:length(S_d)
        ptcols <- unlist(mapply(function(n1, n2) rep(n1, n2), ptcol, neach))
        plotNs <- c(1:5, N)
        ht <- 8
        wd <- 12
        labelsize <- 2.5
        ticksize <- 2
        if(save_plots) {
            pdf(filename, height = ht, width = wd, family = "sans")
        } else dev.new(height = ht, width = wd)
        par(mfrow = c(2, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
        xlim <- range(unlist(S_ews))
        ylim <- c(0, 3.4)#range(unlist(S_d))
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
                yaxt = "n",
                cex.axis = ticksize, cex.lab = labelsize#, font.lab = 3
            )
            axis(2, axTicks(2)[c(TRUE, FALSE)], cex.axis = ticksize)
            abline(v = max_x, col = "gray", lwd = .9, lty = 2)
            abline(h = max_y, col = "gray", lwd = .9, lty = 2)
            mtext(paste("n =", plotNs[i]), line = -2.5, adj = 0.02, font = 1, cex = .75*labelsize)
        }
        if(save_plots) dev.off()
    }
)
