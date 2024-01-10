df <- read.csv("./tau-tau-data-highlowinput.csv")
report <- FALSE # TRUE

printreport <- function(net) {
    input <- df[df$net == net & df$metric == "input", ]
    d <- df[df$net == net & df$metric == "d", ]
    summary(abs(input$tau) - abs(d$tau))
}

plotit <- function(net) {
    imagefile <- paste0("./tau-tau-highlowinput-", net, ".pdf")
    
    pdf(imagefile, height = 8, width = 12)
    par(mfrow = c(2, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
    palette("Tableau 10")

    labelsize <- 2.5
    ticksize <- 2

    ylim <- xlim <- range(abs(df[df$net == net, "tau"]))
    
    for(i in 1:5) {

        plotdf <- df[df$n == i & df$net == net, ]
        
        x <- abs(plotdf$tau[plotdf$metric == "input"])
        y <- abs(plotdf$tau[plotdf$metric == "d"])

        condition <- apply(plotdf[, c("varyv", "varys")], 1, function(row) {
            if(row[1] == FALSE && row[2] == FALSE) {
                return("const")
            } else if(row[1] == TRUE && row[2] == FALSE) {
                return("varyv")
            } else if(row[2] == TRUE) {
                return("varynoise")
            } else {
                return(NA)
            }
        })

        plotdf$cond <- factor(condition, levels = c("const", "varyv", "varynoise"))
        pchs <- as.integer(plotdf$cond[plotdf$metric == "d"]) + 20

        plotdf$bparam <- factor(plotdf$bparam, levels = c("u", "D"))
                
        colors <- as.integer(factor(plotdf$model[plotdf$metric == "d"]))
        bgs <- mapply(
            function(color, alpha) adjustcolor(color, alpha.f = alpha),
            colors, -(as.integer(plotdf$bparam[plotdf$metric == "d"]) - 2)*.5
        )
        
        plot(
            x, y, col = colors, pch = pchs, lwd = 2, cex = 3,
            bg = bgs,
            xlab = "High/Low Input", ylab = "",
            xlim = xlim, ylim = ylim,
            cex.lab = labelsize, cex.axis = ticksize
        )
        title(ylab = "d", font.lab = 3, cex.lab = labelsize)
        abline(a = 0, b = 1, lty = 1, col = "black", lwd = .5)
        
        mtext(paste0("(", letters[i], ")"), line = -2.5, adj = 0.02, cex = 2)
    }

    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    legendcolors <- c(1:4, rep("gray", 5))
    legendalphas <- c(rep(.5, 8), 0)
    legend(
        "left", bty = "n", cex = 1.75, pt.cex = 3, pt.lwd = 3,
        col = legendcolors, pch = c(rep(21, 4), 21:23, 21, 21),
        pt.bg = mapply(function(color, alpha) adjustcolor(color, alpha), legendcolors, legendalphas),
        legend = c("Double-well", "Gene regulatory", "Mutualistic species", "SIS",
                   "Homogeneous", "Heterogeneous stress", "Heterogeneous stress and noise",
                   "Bifurcation parameter is u", "Bifurcation parameter is D")
    )

    dev.off()
}

for(net in unique(df$net)) {
    plotit(net)

    if(report) printreport(net)
}
