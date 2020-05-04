
plot.umap = function(x, labels,
                     main="A UMAP visualization ",
                     colors=rainbow(19),
                     pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
                     cex.main=1, cex.legend=1) {
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  } 
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}
set.seed(123456)


iris.data = iris[, grep("Sepal|Petal", colnames(iris))]
iris.labels = iris[, "Species"]

iris.umap = umap(iris.data)
iris.umap
plot.iris(iris.umap, iris.labels)

#roots.data <- t(assay(ddsMatrix))
roots.data <- t(assay(vsd.QA))
roots.treatment.labels <- meta$Treatment
roots.week.labels <- meta$Week
roots.umap = umap(roots.data)

plot.umap(roots.umap, roots.treatment.labels, colors=c("#009444", "#BE9230"))
plot.umap(roots.umap, roots.week.labels, colors=rainbow(24))
