figure_root <- "new/figures"
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
pal <- c(Rel = "#0072B2", Pearson = "#D55E00", RMSE = "#009E73", Random = "#777777")
export_plot <- function(name, draw, width = 7.2, height = 4.2) {
  pdf(file.path(figure_root, paste0(name, ".pdf")), width = width, height = height, family = "Times", useDingbats = FALSE)
  draw(); dev.off()
  png(file.path(figure_root, paste0(name, ".png")), width = width, height = height, units = "in", res = 300, type = "cairo", family = "serif")
  draw(); dev.off()
}
panel_letter <- function(x) mtext(x, side = 3, adj = 0, line = .3, font = 2, cex = 1.1)
read_result <- function(path) read.csv(file.path("new/results", path), check.names = FALSE)
