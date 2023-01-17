
##
#
# This script builds a figure displaying the results of the benchmark script.
#
##


# Libraries

library(dplyr)
library(ggplot2)
library(ggforce)
library(scico)
library(grid)
library(patchwork)


# Parameters

theme_set(theme_bw(base_size = 14))


# Load data

threads_data <- read.table(
  file = "benchmark/benchmark_29.12.22_threads",
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
size_data <- read.table(
  file = "benchmark/benchmark_29.12.22_size",
  header = T,
  sep = "\t",
  stringsAsFactors = F
)


# build thread figure

median_thread <- median(unique(threads_data$threads))
max_thread <- max(threads_data$threads)
n_threads <- length(unique(threads_data$threads))

threads_plot_data <- threads_data %>% 
  mutate(
    x = log10(peptides),
    dx = 1/4*(threads - median_thread)/max_thread,
    x_bin = factor(x, levels = 1:6),
    time = ifelse(time == 0, 1, time),
    y = log10(time * threads),
    col = factor(threads)
  )

plot_threads <- ggplot() +
  geom_abline(
    slope = 1,
    intercept = -1,
    linetype = "dashed"
  ) +
  geom_point(
    data = threads_plot_data,
    mapping = aes(
      x = x + dx,
      y = y,
      col = col
    ),
    alpha = 0.8
  ) +
  scale_x_continuous(
    name = "# Peptides",
    breaks = 1:6,
    labels = 10^(1:6)
  ) +
  scale_y_continuous(
    name = "Time * threads",
    breaks = log10(c(10, 100, 1000, 10 * 1000, 60 * 1000)),
    labels = c("10 ms", "100 ms", "1 s", "10 s", "1 min")
  ) +
  scale_color_manual(
    name = "Threads",
    values = scico(
      n = n_threads,
      palette = "imola",
      begin = 0.2,
      end = 0.7
      )
  )+
  guides(
    color = guide_legend(
      override.aes = list(
        alpha = 1,
        size = 2
      )
    )
  )  + 
  theme(
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )

png(
  filename = "benchmark/benchmark_29.12.22_threads.png",
  width = 400,
  height = 300
)
grid.draw(plot_threads)
device <- dev.off()


# build size figure

size_plot_data <- size_data %>% 
  mutate(
    x = sites * occupancy,
    y = time / 1000,
    col = factor(modifications)
  )

size_plot <- ggplot() +
  geom_point(
    data = size_plot_data,
    mapping = aes(
      x = x,
      y = y,
      col = col
    ),
    alpha = 0.5
  ) +
  scale_x_continuous(
    name = "Possible * occupied sites"
  ) +
  scale_y_continuous(
    name = "Time [ms/peptide]",
  ) + 
  scale_color_scico_d(
    name = "# Mods",
    palette = "bilbao",
    begin = 0.5,
    end = 1,
    direction = -1
  ) +
  guides(
    color = guide_legend(
      reverse = T,
      override.aes = list(
        alpha = 1,
        size = 2
        )
      )
  ) +
  theme(
    panel.border = element_blank(),
    legend.position = "right"
  )

png(
  filename = "benchmark/benchmark_29.12.22_size.png",
  width = 400,
  height = 300
)
grid.draw(size_plot)
device <- dev.off()

