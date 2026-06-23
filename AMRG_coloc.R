library(ggplot2)
library(cowplot)
library(dplyr)

AMRGbyP # col1=Plasmid ; col2= Genus; col3=ARG count ; col4= MRG count
set.seed(123) 
## List of unique plasmids for permutation reference
unique_plasmids <- unique(AMRGbyP$Plasmid)
n_plasmids <- length(unique_plasmids)

# Extract original numeric vectors (one per unique plasmid)
plasmid_ARG <- AMRGbyP %>% group_by(Plasmid) %>% summarize(ARG = unique(ARG), .groups='drop') %>% pull(ARG)
plasmid_MRG <- AMRGbyP %>% group_by(Plasmid) %>% summarize(MRG = unique(MRG), .groups='drop') %>% pull(MRG)

# Calculate total units
total_ARG_units <- sum(plasmid_ARG)
total_MRG_units <- sum(plasmid_MRG)
cat("Total ARG units:", total_ARG_units, "| Total MRG units:", total_MRG_units, "\n")

n_perm <- 10000
perm_probs_ARGifMRG <- numeric(n_perm)

for(i in 1:n_perm) {
  # Multinomial unit redistribution 
  perm_ARG <- as.numeric(rmultinom(1, size = total_ARG_units, prob = rep(1/n_plasmids, n_plasmids)))
  perm_MRG <- as.numeric(rmultinom(1, size = total_MRG_units, prob = rep(1/n_plasmids, n_plasmids)))
  
  # Build a mapping dataframe for this permutation
  perm_map <- data.frame(
    Plasmid = unique_plasmids,
    ARG = perm_ARG,
    MRG = perm_MRG
  )
  
  # Assign values back to all rows by Plasmid
  perm_AMRGbyP <- AMRGbyP %>%
    select(-ARG, -MRG) %>%
    left_join(perm_map, by = "Plasmid")
  
  # Aggregate at Genus level as before
  perm_genus_summary <- perm_AMRGbyP %>%
    group_by(Genus) %>%
    summarize(ARG = sum(ARG), MRG = sum(MRG), .groups='drop')
  
  perm_probs_ARGifMRG[i] <- sum(perm_genus_summary$ARG > 0 & perm_genus_summary$MRG > 0) / sum(perm_genus_summary$MRG > 0)
}

#Count Genera
obs_genus_summary <- AMRGbyP %>%
  group_by(Genus) %>%
  summarize(ARG = sum(ARG), MRG = sum(MRG), .groups='drop')
# Calculate observed values
# ARG & MRG if MRG
prob_ARGifMRG <- sum(obs_genus_summary$ARG > 0 & obs_genus_summary$MRG > 0) / sum(obs_genus_summary$MRG > 0)
prob_ARGifMRG
#ARG and MRG
prob_both_obs <-
  sum(obs_genus_summary$ARG > 0 & obs_genus_summary$MRG > 0) /
  nrow(obs_genus_summary)
prob_both_obs
# Empirical p-value 
p_val_ARGifMRG <- mean(perm_probs_ARGifMRG > prob_ARGifMRG)
cat("Observed prob_ARGifMRG:", prob_ARGifMRG, "\nP-value:", p_val_ARGifMRG, "\n")

# Plot 
AMRGbyP_ARGifMRG <- data.frame(perm_probs_ARGifMRG)
ARGifMRGprobplot <- ggplot(AMRGbyP_ARGifMRG, aes(x = perm_probs_ARGifMRG)) + 
  geom_histogram(bins = 50, alpha = 0.7) + 
  geom_vline(xintercept = prob_ARGifMRG,
             color = "blue", linetype = "dashed", size = 1.2) +
  geom_vline(xintercept = quantile(perm_probs_ARGifMRG, 0.95),
             color = "blue", linetype = "dotted") +
  theme_minimal() +
  theme(
    text       = element_text(size = 20),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(color = "black"),
    plot.margin = unit(c(2, 2, 1, 2), "cm")
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(
    x = "Prob. ARG on MRG-encoding plasmids\n(Unit redistribution; 10,000 perms)", 
    y = "Count",
    title = paste(
      "Pobs. =", round(prob_ARGifMRG, 3),
      "\np-value =", round(p_val_ARGifMRG, 3)
    )
  )

  print(ARGifMRGprobplot)
