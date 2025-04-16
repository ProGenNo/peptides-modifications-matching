
##
#
# This script builds a figure displaying the modification matching for the histone example.
#
##


# Libraries

library(dplyr)
library(ggplot2)
library(ggforce)
library(scico)
library(grid)
library(janitor)
library(patchwork)
library(here)


# Parameters

theme_set(theme_minimal(base_size = 22))


# Load data

modification_occurrence <- read.table(
  file = here("histone/modification_occurrence.gz"),
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ''
)

modification_weight <- read.table(
  file = here("histone/modification_weight.gz"),
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ''
)

# Plot all PSMs

# psm_id <- 2059

for (psm_id in unique(modification_occurrence$psm)) {
  
  psm_occurrence <- modification_occurrence %>% 
    filter(
      psm == psm_id
    )
  
  psm_weight <- modification_weight %>% 
    filter(
      psm == psm_id
    )
  
  sites <- sort(unique(psm_weight$site))
  
  first_site <- psm_weight %>% 
    group_by(
      modification
    ) %>% 
    summarise(
      first_site = min(site)
    ) %>% 
    arrange(
      first_site
    )
  
  modifications <- first_site$modification
  
  n_sites <- length(sites)
  n_modifications <- length(modifications)
  
  psm_occurrence <- psm_occurrence %>% 
    mutate(
      modifiation_label = as.character(round(modification, digits = 2)),
      modifiation_label = ifelse(occurrence > 1, paste0(modifiation_label, " (", occurrence, ")"), modifiation_label),
      modification_factor = factor(modification, levels = modifications),
      x = as.numeric(modification_factor) / (n_modifications + 1)
    )
  
  psm_weight <- psm_weight %>% 
    mutate(
      site_factor = factor(site, levels = sites),
      x_site = as.numeric(site_factor)/(n_sites + 1),
      selected_factor = factor(selected, levels = c(1, 0))
    ) %>% 
    left_join(
      psm_occurrence %>% 
        select(
          modification,
          x_modification = x
        ),
      by = "modification"
    )
  
  psm_sites <- psm_weight %>% 
    select(
      site,
      x_site
    ) %>% 
    distinct()
  
  assinment_plot <- ggplot() +
    geom_text(
      data = psm_sites,
      mapping = aes(
        x = x_site,
        y = 2,
        label = site
      ),
      vjust = 0,
      size = 8
    ) +
    geom_text(
      data = psm_occurrence,
      mapping = aes(
        x = x,
        y = 1,
        label = modifiation_label
      ),
      vjust = 1,
      size = 8
    ) +
    geom_segment(
      data = psm_weight,
      mapping = aes(
        x = x_site,
        xend = x_modification,
        y = 1.95,
        yend = 1.05,
        col = weight,
        linetype = selected_factor
      )
    ) +
    scale_color_scico(
      name = "Score",
      palette = "grayC",
      limits = c(0, 1)
    ) + 
    scale_linetype_manual(
      name = NULL,
      values = c("solid", "dotted"),
      labels = c("Selected", "Rejected")
    ) +
    ggtitle(
      psm_occurrence$sequence[1]
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    )
  
  png(
    filename = here(paste0("histone/plots/", psm_id, ".png")),
    width = 900,
    height = 600
  )
  grid.draw(assinment_plot)
  device <- dev.off()
  
}


# Select complex assignment examples

score_per_psm <- modification_weight %>% 
  select(
    psm, weight
  ) %>% 
  distinct() %>% 
  group_by(
    psm
  ) %>% 
  summarize(
    n_scores = n()
  ) %>% 
  arrange(
    desc(n_scores)
  )

modification_overlap <- modification_weight %>% 
  select(
    psm, site, modification
  ) %>% 
  distinct() %>% 
  group_by(
    psm, site
  ) %>% 
  summarize(
    mods_per_site = n()
  ) %>% 
  group_by(
    psm
  ) %>% 
  summarize(
    max_mods_per_site = max(mods_per_site)
  ) %>% 
  filter(
    max_mods_per_site > 1
  )

# Plot examples

plot_psm <- function(psm_id, label) {
  
  psm_occurrence <- modification_occurrence %>% 
    filter(
      psm == psm_id
    )
  
  psm_weight <- modification_weight %>% 
    filter(
      psm == psm_id
    )
  
  sites <- sort(unique(psm_weight$site))
  
  first_site <- psm_weight %>% 
    group_by(
      modification
    ) %>% 
    summarise(
      first_site = min(site)
    ) %>% 
    arrange(
      first_site
    )
  
  modifications <- first_site$modification
  
  n_sites <- length(sites)
  n_modifications <- length(modifications)
  
  psm_occurrence <- psm_occurrence %>% 
    mutate(
      modifiation_label = as.character(round(modification, digits = 2)),
      modifiation_label = ifelse(occurrence > 1, paste0(modifiation_label, " (", occurrence, ")"), modifiation_label),
      modification_factor = factor(modification, levels = modifications),
      x = as.numeric(modification_factor) / (n_modifications + 1)
    )
  
  psm_weight <- psm_weight %>% 
    mutate(
      site_factor = factor(site, levels = sites),
      x_site = as.numeric(site_factor)/(n_sites + 1),
      selected_factor = factor(selected, levels = c(1, 0))
    ) %>% 
    left_join(
      psm_occurrence %>% 
        select(
          modification,
          x_modification = x
        ),
      by = "modification"
    )
  
  psm_sites <- psm_weight %>% 
    select(
      site,
      x_site
    ) %>% 
    distinct()
  
  assinment_plot <- ggplot() +
    geom_text(
      data = psm_sites,
      mapping = aes(
        x = x_site,
        y = 2,
        label = site
      ),
      vjust = 0,
      size = 8
    ) +
    geom_text(
      data = psm_occurrence,
      mapping = aes(
        x = x,
        y = 1,
        label = modifiation_label
      ),
      vjust = 1,
      size = 8
    ) +
    geom_segment(
      data = psm_weight,
      mapping = aes(
        x = x_site,
        xend = x_modification,
        y = 1.95,
        yend = 1.05,
        col = weight,
        linetype = selected_factor
      )
    ) +
    scale_color_scico(
      name = "Score",
      palette = "grayC",
      limits = c(0, 1)
    ) + 
    scale_linetype_manual(
      name = NULL,
      values = c("solid", "dotted"),
      labels = c("Selected", "Rejected"),
      drop = F
    ) +
    scale_y_continuous(
      expand = c(0, 0.15)
    ) +
    ggtitle(
      label = label,
      subtitle = psm_occurrence$sequence[1]
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    )
  
  return(assinment_plot)
  
}

plot_1 <- plot_psm(6353, "A")
plot_2 <- plot_psm(2059, "B")

png(
  filename = here("histone/fig_histones.png"),
  width = 900,
  height = 600
)
grid.draw(
  plot_1 / plot_2 + plot_layout(guides = 'collect')
)
device <- dev.off()



