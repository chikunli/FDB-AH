entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(entry)), "lib", "init.R"))
source(file.path(repro_dir, "lib", "real.R"))
suppressPackageStartupMessages(library(ggplot2))
library(grid)
if (!capabilities("cairo")) stop("Figure1 requires an R installation with Cairo PDF support.")

normal <- run_record("normal", figure = TRUE)
rotor <- run_record("rotor", figure = TRUE)
methods <- c("FDB-AH", "FDB-HF", "FDB-75")
colours <- c("FDB-AH" = "#D55E00", "FDB-HF" = "#444444", "FDB-75" = "#009E73")
shapes <- c("FDB-AH" = 16, "FDB-HF" = 15, "FDB-75" = 18)
lines <- c("FDB-AH" = "dashed", "FDB-HF" = "dashed", "FDB-75" = "longdash")
font <- "Times New Roman"
theme_set(theme_classic(base_size = 9, base_family = font) +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.title = element_text(size = 8.5),
        strip.background = element_blank(), strip.text = element_text(size = 9, hjust = 0),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 8.5, margin = margin(l = 2)),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(3, "pt"), legend.key.height = unit(10, "pt"),
        legend.key.width = unit(4, "mm"), legend.spacing.x = unit(3, "mm"),
        legend.key.spacing.x = unit(5, "mm"),
        panel.spacing = unit(10, "pt"), plot.margin = margin(6, 8, 4, 4)))

path_panel <- function(run, title, ticks) {
  path <- transform(run$path, panel = title)
  first <- subset(run$scores, seed == 1L, select = c(method, h))
  marks <- merge(first, path, by = "h", all.x = TRUE)
  marks$method <- factor(marks$method, levels = methods)
  stopifnot(all(is.finite(marks$d_next_sq)))
  ggplot(path, aes(h, d_next_sq)) +
    geom_vline(data = marks, aes(xintercept = h, colour = method, linetype = method),
               linewidth = 0.32, alpha = 0.7) +
    geom_line(linewidth = 0.4, colour = "#2B4C6F") +
    geom_point(data = marks, aes(colour = method, shape = method), size = 1.8) +
    geom_text(data = subset(marks, method == "FDB-AH"),
              aes(label = paste0("h = ", h), colour = method),
              family = font, size = 2.6, hjust = 1.1, vjust = -0.8) +
    facet_wrap(~panel) +
    scale_y_log10(expand = expansion(mult = c(0.08, 0.22))) +
    scale_x_continuous(breaks = ticks, expand = expansion(mult = c(0.025, 0.025))) +
    scale_colour_manual(values = colours) + scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = lines) +
    labs(x = expression("Initial " * italic(h)), y = "(h + 1)st squared\nMahalanobis distance") +
    theme(legend.position = "none")
}

rate_panel_data <- function(run, metric, title) {
  do.call(rbind, lapply(methods, function(m) {
    rates <- run$scores[[metric]][run$scores$method == m]
    data.frame(method = m, panel = title, value = mean(rates), low = min(rates), high = max(rates))
  }))
}
rates <- rbind(rate_panel_data(normal, "FPR", "c  Normal: FPR"),
               rate_panel_data(rotor, "FPR", "d  Anomalous: FPR"),
               rate_panel_data(rotor, "FNR", "e  Anomalous: FNR"))
rates$method <- factor(rates$method, levels = methods)
lower <- ggplot(rates, aes(method, value, colour = method, shape = method)) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.13, linewidth = 0.45, show.legend = FALSE) +
  geom_point(size = 2.3) +
  geom_text(aes(y = high + 0.85, label = sprintf("%.3f", value)),
            family = font, size = 2.65, show.legend = FALSE) +
  facet_wrap(~panel, ncol = 3) + scale_colour_manual(values = colours) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(limits = c(0, 18), breaks = c(0, 5, 10, 15), expand = expansion(mult = c(0.025, 0))) +
  scale_x_discrete(expand = expansion(add = 0.6)) + labs(x = NULL, y = "Error rate (%)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(0, 8, 4, 4),
        strip.text = element_text(size = 9, hjust = 0, margin = margin(0, 4.5, 2.5, 4.5)))

upper_a <- path_panel(normal, "a  SKAB: normal", c(750, 1000, 1250, 1500))
upper_b <- path_panel(rotor, "b  SKAB: anomalous", c(550, 700, 850, 1000))
grDevices::cairo_pdf(file.path(output_dir, "Figure1.pdf"), width = 6.5, height = 4.9, family = font)
grid.newpage()
upper_height <- 5.3 * 0.54
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(upper_height, 4.9 - upper_height), "in"))))
print(upper_a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(upper_b, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(lower, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
popViewport()
invisible(grDevices::dev.off())
