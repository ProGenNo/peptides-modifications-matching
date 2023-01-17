
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


# Parameters

theme_set(theme_void(base_size = 14))
sequence <- "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE"

# Load data

modification_occurrence <- read.table(
  file = "histone/modification_occurrence.gz",
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ''
)

modification_weight <- read.table(
  file = "histone/modification_weight.gz",
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ''
)


# Plot one example

cpt <- 0

ids <- c(20, 8, 1537, 35)

unique_ids <- unique(modification_occurrence$psm)

# for (id in sample(unique_ids, length(unique_ids))) {
for (id in ids) {
  
  modifications_at_id <- as.character(unique(modification_occurrence$modification[modification_occurrence$psm == id]))
  
  if (length(modifications_at_id) == 2) {
    
    psm_mapping <- modification_occurrence %>% 
      filter(
        psm == id
      ) %>% 
      select(
        -psm
      ) %>% 
      left_join(
        modification_weight %>% 
          filter(
            psm == id
          ) %>% 
          select(
            -psm
          ),
        by = "modification"
      ) %>% 
      mutate(
        modification = as.character(modification),
        xend = NA,
        xend = ifelse(modification == modifications_at_id[1], -1, 1),
        modification_name = ifelse(modification == "28.03130012828", "Dimethylation", modification),
        modification_name = ifelse(modification == "42.0105646837", "Acetylation", modification_name),
        modification_name = ifelse(modification == "42.0105646837", "Acetylation", modification_name),
        modification_name = ifelse(modification == "79.96633052075", "Phosphorylation", modification_name),
        modification_name = ifelse(modification == "14.01565006414", "Mehtylation", modification_name),
        modification_name = ifelse(modification == "42.04695019242", "Trimehtylation", modification_name),
        modification_name = paste0(modification_name, " (", occurrence, ")"),
        weight_round = round(weight),
        site = ifelse(site == 0, 1, site),
        selected = factor(selected)
      )
    
    sequence_df <- data.frame(
      position = 1:nchar(sequence) - nchar(sequence) / 2
    )
    
    for (i in 1:nchar(sequence)) {
      
      sequence_df$letter[i] <- substr(sequence, i, i)
      
    }
    
    min_max <- psm_mapping %>% 
      select(xend, modification_name, site) %>% 
      group_by(
        modification_name
      ) %>% 
      summarize(
        site_min = min(site),
        site_max = max(site),
        xend = min(xend)
      )
    
    label_df <- psm_mapping %>% 
      select(
        modification_name, xend
      ) %>% 
      distinct()
    
    sequence_plot <- ggplot() +
      geom_text(
        data = sequence_df,
        mapping = aes(
          x = position,
          y = 0,
          label = letter
        )
      ) +
      geom_text(
        data = label_df,
        mapping = aes(
          x = 0,
          y = xend,
          label = modification_name
        )
      ) +
      geom_segment(
        data = psm_mapping,
        mapping = aes(
          x = site - nchar(sequence) / 2,
          xend = site - nchar(sequence) / 2,
          y = xend / 10,
          yend = 4 * xend / 10
        )
      ) +
      geom_text(
        data = psm_mapping,
        mapping = aes(
          x = site - nchar(sequence) / 2,
          y = xend / 2,
          label = weight_round,
          color = selected
        )
      ) +
      geom_segment(
        data = psm_mapping,
        mapping = aes(
          x = site - nchar(sequence) / 2,
          xend = site - nchar(sequence) / 2,
          y = 6 * xend / 10,
          yend = 9 * xend / 10
        )
      ) +
      geom_segment(
        data = min_max,
        mapping = aes(
          x = site_min - nchar(sequence) / 2,
          xend = site_max - nchar(sequence) / 2,
          y = 9 * xend / 10,
          yend = 9 * xend / 10
        )
      ) +
      scale_color_manual(
        values = c("grey", "black")
      ) +
      theme(
        legend.position = "none"
      )
    
    png(
      filename = paste0("histone/figure_histone_", id, ".png"),
      height = 150,
      width = 700
    )
    grid.draw(sequence_plot)
    device <- dev.off()
    
    cpt <- cpt + 1
    
    if (cpt >= 10) {
      
      break()
      
    }
  }
}




