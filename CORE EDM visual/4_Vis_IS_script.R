## ESU

par(mar = c(1, 0, 1, 5))

radarchart(rec4_spider, axistype=0,
           #custom polygon
           pcol=pal, pfcol=trans.pal,  plwd=2, plty=1, seg = 3,
           #custom the grid
           cglcol="grey", cglty=1, cglwd=0.8,
           #custom labels
           vlcex=.9, vlabels = c("ARC", "PDO", "NPGO", "UPW"), 
           title="What variables are found in top models?")

legend(x=.9, y=.8, legend = c("Ave Rank", "No. Models", "Highest Rank"), bty = "n", pch=20 , col=pal, text.col = "grey", cex=1, pt.cex=2)


par(op)

## IMN


radarchart(rec4_spider, axistype=0,
           #custom polygon
           pcol=pal, pfcol=trans.pal,  plwd=2, plty=1, seg = 3,
           #custom the grid
           cglcol="grey", cglty=1, cglwd=0.8,
           #custom labels
           vlcex=.9, vlabels = c("CSL", "NPGO", "PDO", "UPW",  "Harvest",
                                 "Hatchery", "ORCA"), 
           title="What variables are found in top models?")


legend(x=.9, y=.7, legend = c("Ave Rank", "No. Models", "Highest Rank"), bty = "n", pch=20 , col=pal, text.col = "grey", cex=1, pt.cex=2)

## MFS

radarchart(rec4_spider, axistype=0,
           #custom polygon
           pcol=pal, pfcol=trans.pal,  plwd=2, plty=1, seg = 3,
           #custom the grid
           cglcol="grey", cglty=1, cglwd=0.8,
           #custom labels
           vlcex=.9, vlabels = c("Hatchery", "NPGO", "PDO", "UPW", "ARC", "Harvest"), 
           title="What variables are found in top models?")
legend(x=.9, y=.5, legend = c("Ave Rank", "No. Models", "Highest Rank"), bty = "n", pch=20 , col=pal, text.col = "grey", cex=1, pt.cex=2)


## UPS
op <- par(mar = c(1, 1, 1, 1))
par(mar = c(1, 1, 1, 1))

radarchart(rec4_spider, axistype=0,
           #custom polygon
           pcol=pal, pfcol=trans.pal,  plwd=2, plty=1, seg = 3,
           #custom the grid
           cglcol="grey", cglty=1, cglwd=0.8,
           #custom labels
           vlcex=.9, vlabels = c("CSL", "NPGO", "PDO", "UPW",  "Harvest",
                                 "Hatchery", "ORCA", "SSL"), 
           title="What variables are found in top models?")

legend(x=.9, y=.8, legend = c("Ave Rank", "No. Models", "Highest Rank"), bty = "n", pch=20 , col=pal, text.col = "grey", cex=1, pt.cex=2)


par(op)


