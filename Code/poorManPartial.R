poorManPartial <- function(mod = NULL,
                           dat = NULL,
                           plot_var = TRUE,
                           varimpFilename = NULL,
                           plot_plotmo = TRUE,
                           plotmoFilename = NULL,
                           plot_pdp = FALSE,
                           partialFilename = NULL) {
  
  require(pdp)
  require(ggplot2)
  require(ggalt)
  require(gridExtra)
  require(plotmo)
  
  vars <- varImp(mod)$importance
  vars$varname <- row.names(vars)
  vars <- vars[rev(order(vars$Overall)), ]
  vars$rank <- 1:nrow(vars)
  
  if (plot_var == TRUE) {
    cat("Plotting variable importance \n")
    p1 <- ggplot(data = vars, aes(x = Overall, y = reorder(varname, Overall))) +
      geom_lollipop(horizontal = TRUE) +
      labs(x = "Variable relative importance", y = "Variable name") +
      theme_bw() +
      theme(axis.text = element_text(colour = "black"))
    ggsave(paste0(varimpFilename),
           p1,
           device = "pdf",
           width = 100, height = 100, units = "mm")
  }
  
  if (plot_plotmo == TRUE) {
    cat("Plotting plotmo \n")
    pdf(file = paste0(plotmoFilename), paper = "a4r")
    plotmo(mod,
           type = "prob",
           all1 = TRUE,
           degree1 = vars$varname,
           ylim=c(0,1),
           caption = paste(this.species, this.stage),
           nrug = 10, rug.ticksize = 0.1, rug.lwd = 1, rug.col = "gray50")
    while (!is.null(dev.list()))  dev.off()
  }
  
  if (plot_pdp == TRUE) {
    cat("Plotting partial dependence \n")
    plts <- list()
    for (i in 1:length(vars$varname)) {
      print(vars$varname[i])
      plts[i] <- list(partial(mod,
                              pred.var = paste(vars$varname[i]),
                              plot = TRUE,
                              prob = TRUE,
                              parallel = FALSE,
                              rug = T,
                              smooth = FALSE,
                              train = dat))
    }
    #dev.off()
    while (!is.null(dev.list()))  dev.off()
    ggsave(paste(partialFilename),
           marrangeGrob(plts,
                        ncol=4,
                        nrow=3,
                        as.table = TRUE,
                        top = paste(this.species, this.stage)),
           device = "pdf",
           width = 297, height = 210, units = "mm")
  }
  
}