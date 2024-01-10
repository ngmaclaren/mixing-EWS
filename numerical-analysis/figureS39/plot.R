library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
palette("Tableau 10")

load("./pdata-alt-u.rda")
pdata_u <- pdata
load("./pdata-alt-D.rda")
pdata_D <- pdata

plotit <- function(plotdata, show_ylabel = TRUE) {
    ns <- 1:5
    p1 <- do.call(rbind, lapply(plotdata, function(x) colMeans(x$props)))
    p2 <- do.call(rbind, lapply(plotdata, function(x) colMeans(x$dists)))
    if(length(plotdata) == 3) ltys <- 1:3 else ltys = c(1, 3)
    matplot(
        ns, 2*t(p1), type = "o", pch = 0, lty = ltys, lwd = 2, cex = 2, col = 1, ylim = c(0, 1.25),
        axes = FALSE, xlab = "", ylab = ""
    )
    box()
    axis(1, cex.axis = 1.75)
    title(xlab = "n", font.lab = 3, cex.lab = 1.75)
    if(show_ylabel) {
        axis(2, cex.axis = 1.75)
        title(ylab = TeX(r"($p_1$, $p_2$)", italic = TRUE), cex.lab = 1.75)
    } else {
        axis(2, tick = TRUE, labels = FALSE)
    }
    matpoints(ns, t(p2), pch = 1, lwd = 2, cex = 2, col = 2)
    matlines(ns, t(p2), lty = ltys, lwd = 2, col = 2)
}

makeplot <- function(networks, labels, filename, singleletter = FALSE) {
    pdf(filename, height = 12, width = 12)
    par(mfrow = c(4, 4), mar = c(4, 4, 0, 0)+0.5)
    idx <- 1
    for(i in 1:2) {# net
        net <- networks[i]
        simgroup <- list(
            u = paste(net, c("doublewell", "mutualistic", "genereg"), sep = "-"),
            D = paste(net, c("doublewell", "mutualistic", "genereg", "SIS"), sep = "-")
        )
        for(j in 1:2) {# bparam
            for(k in 1:4) {
                if(k %in% 2:4) show_ylabel <- FALSE else show_ylabel <- TRUE
                
                if(k == 4 && j == 1 && i == 1) {
                    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
                    legend(
                        "left", bty = "n", cex = 1.5, pt.cex = 1.75, pt.lwd = 2, lwd = 2,
                        legend = c("p1", "p2", "Homogeneous", "Heterogeneous stress",
                                   "Heterogeneous stress\nand noise"),
                        pch = c(0, 1, NA, NA, NA), lty = c(NA, NA, 1:3), col = c(1, 2, 10, 10, 10)
                    )
                } else if(k == 4 && j == 1 && i == 2) {
                    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
                } else {
                    if(j == 1) {
                        plotdata <- pdata_u[grep(simgroup$u[k], names(pdata_u))]
                        plotit(plotdata, show_ylabel = show_ylabel)
                        labeler(idx, labels, singleletter)
                        idx <- idx + 1
                    } else {
                        plotdata <- pdata_D[grep(simgroup$D[k], names(pdata_D))]
                        plotit(plotdata, show_ylabel = show_ylabel)
                        labeler(idx, labels, singleletter)
                        idx <- idx + 1
                    }
                }
            }
        }
    }
    dev.off()
}

labeler <- function(idx, labels, singleletter = FALSE) {
    if(singleletter) {
        mtext(labels[idx], line = -2.2, adj = 0.02, font = 2, cex = 2)
    } else {
        mtext(labels[idx], line = -2.2, adj = 0.02, font = 1, cex = 1.5)
    }
}

formain <- c("ba", "chesapeake")

makeplot(
    formain,
    paste0("(", letters, ")"),
    "./overall-summary-alt.pdf",
    singleletter = TRUE
)
