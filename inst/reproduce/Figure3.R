entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(entry)), "lib", "init.R"))
sim <- new.env(parent = baseenv())
sys.source(file.path(repro_dir, "lib", "simulation.R"), sim)
suppressPackageStartupMessages(library(ggplot2))
if (!capabilities("cairo")) stop("Figure 3 requires an R installation with Cairo PDF support.")

design <- subset(sim$simulation_design(), study == "curve")
raw <- sim$simulate(design, repro_map)
groups <- split(raw, raw$id)
design$mean_h <- vapply(as.character(design$id), function(id) {
  rows <- groups[[id]]
  stopifnot(identical(sort(rows$rep), 1:1000))
  mean(rows$h)
}, numeric(1L))
utils::write.csv(design[c("setting", "eps", "type", "r", "mean_h")],
                 file.path(output_dir, "Figure3.csv"), row.names = FALSE)

types <- names(sim$type_labels)
design$type <- factor(design$type, levels = types)
design$panel <- factor(design$eps, levels = c(0.1, 0.3),
                       labels = c("a   \u03b5 = 0.10", "b   \u03b5 = 0.30"))
reference <- unique(design[c("panel", "eps", "n")])
reference$regular <- reference$n * (1 - reference$eps)
font <- "Times New Roman"
plot <- ggplot(design, aes(r, mean_h, colour = type, shape = type, linetype = type)) +
  geom_hline(data = reference, aes(yintercept = regular),
             colour = "grey45", linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.4) + geom_point(size = 1.8) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c("#D55E00", "#CC79A7", "#009E73", "#0072B2"),
                      breaks = types, labels = sim$type_labels) +
  scale_shape_manual(values = c(16, 17, 15, 18), breaks = types, labels = sim$type_labels) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "dotdash"),
                        breaks = types, labels = sim$type_labels) +
  scale_x_continuous(breaks = 6:10, expand = expansion(mult = c(0.035, 0.04))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(x = "Signal parameter r", y = "Mean initial subset size") +
  guides(colour = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 9, base_family = font) +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(margin = margin(t = 2.25)),
        axis.text.y = element_text(margin = margin(r = 2.25)),
        axis.title = element_text(size = 8.5),
        axis.title.x = element_text(margin = margin(t = 2.25)),
        axis.title.y = element_text(margin = margin(r = 2.25)),
        axis.line = element_line(linewidth = 9/22, colour = "black"),
        axis.ticks = element_line(linewidth = 9/22, colour = "black"),
        axis.ticks.length = grid::unit(2.25, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, hjust = 0, margin = margin(0, 4.5, 2.5, 4.5)),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 8.5, margin = margin(l = 2)),
        legend.key.width = grid::unit(4, "mm"), legend.key.height = grid::unit(10, "pt"),
        legend.key.spacing.x = grid::unit(5, "mm"), legend.spacing.x = grid::unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(7, "pt"), panel.spacing = grid::unit(10, "pt"),
        panel.grid = element_blank(), plot.margin = margin(6, 8, 4, 4))
ggsave(file.path(output_dir, "Figure3.pdf"), plot = plot, device = grDevices::cairo_pdf,
       width = 6.5, height = 3.15, units = "in", family = font)
