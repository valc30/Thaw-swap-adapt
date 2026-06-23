##############################################################################
Probability of finding the same taxonomic class or sample type within modules.
##############################################################################
library(ggplot2)
plasmid_metadata # Module, Orignial_sample, Sample_type of original sample, taxonomy of plasmid contig
classplasmid_metadata # Module, Orignial_sample, Sample_type of original sample, taxonomy of plasmid contigs assigned to a class 

#1 Create contingency tables
df1 <- table(classplasmid_metadata$Module, classplasmid_metadata$iPHoP_Order)
df2 <- table(plasmid_metadata$Module, plasmid_metadata$Sample_type)
df3 <- table(plasmid_metadata$Module, plasmid_metadata$Original_sample)

df1 # rows = modules, columns = taxa ; number of plasmids in each module
df2 # rows = modules, columns = sample-type ; number of plasmids in each module
df3 # rows = modules, columns = sample-origin ; number of plasmids in each module

###Module and taxonomic class
#1 Compute observed chi-square statistic
chisq_stat <- function(df) {
  expected <- outer(rowSums(df), colSums(df)) / sum(df)
  sum((df - expected)^2 / expected)
}

obs_stat_order <- chisq_stat(df1)
obs_stat_order

#2 permutation test
set.seed(123)

nperm <- 10000
perm_stats_order <- numeric(nperm)

for (i in 1:nperm) {
  perm_Order <- sample(classplasmid_metadata$iPHoP_Order)   # shuffle labels
  perm_tab <- table(classplasmid_metadata$Module, perm_Order)
  perm_stats_order[i] <- chisq_stat(perm_tab)
}

p_value_order <- mean(perm_stats_order >= obs_stat_order)
p_value_order

perm_df_order <- data.frame(chi2 = perm_stats_order)
ggplot(perm_df_order, aes(x = chi2)) +
  geom_histogram(bins = 40, fill = "grey80", color = "white") +
  geom_vline(xintercept = obs_stat_order, color = "blue", linetype = "dashed", size = 1.2) +
  annotate("text",
           x = obs_stat_order,
           y = Inf,
           label = paste("Observed χ² =", round(obs_stat_order, 2),
                         "\np =", signif(p_value_order, 3)),
           vjust = 2,
           hjust = 1.1,
           color = "black") +
  labs(title = expression("iPHoP order permutation test for " * chi^2),
       x = expression(chi^2),
       y = "Frequency") +
  theme_minimal()


#3 proportion of modules containing only 1 class
mean(apply(df1, 1, function(x) sum(x > 0) == 1))

mod_class <- data.frame(
  module = rownames(df1),
  class  = apply(df1, 1, function(x) {
    cl <- colnames(df1)[x > 0]
    paste(cl, collapse = ",")
  }),
  row.names = NULL
)

head(mod_class)
write.csv(mod_class, "module_class.csv")

#5 size effect: Cramer
n_order <- sum(df1)
r_order <- nrow(df1)
c_order <- ncol(df1)

cramers_v_order <- sqrt(obs_stat_order / (n_order * (min(r_order - 1, c_order - 1))))
cramers_v_order

#Taxonomic class shows near-complete association with module structure (Cramér’s V = 0.98), indicating that modules are almost entirely taxonomically homogeneous.



###Module_sample_type 
#2 Compute observed chi-square statistic
obs_stat_type <- chisq_stat(df2)
obs_stat_type

#3 permutation test
set.seed(123)

nperm <- 10000
perm_stats_type <- numeric(nperm)

for (i in 1:nperm) {
  perm_sample_type <- sample(plasmid_metadata$Sample_type)   # shuffle labels
  perm_tab <- table(plasmid_metadata$Module, perm_sample_type)
  perm_stats_type[i] <- chisq_stat(perm_tab)
}

p_value_type <- mean(perm_stats_type >= obs_stat_type)
p_value_type

perm_df_type <- data.frame(chi2 = perm_stats_type)
ggplot(perm_df_type, aes(x = chi2)) +
  geom_histogram(bins = 40, fill = "grey80", color = "white") +
  geom_vline(xintercept = obs_stat_type, color = "blue", linetype = "dashed", size = 1.2) +
  annotate("text",
           x = obs_stat_type,
           y = Inf,
           label = paste("Observed χ² =", round(obs_stat_type, 2),
                         "\np =", signif(p_value_type, 3)),
           vjust = 2,
           hjust = 1.1,
           color = "black") +
  labs(title = expression("Sample-type permutation test for " * chi^2),
       x = expression(chi^2),
       y = "Frequency") +
  theme_minimal()


#4 proportion of modules containing only 1 sample-type
mean(apply(df2, 1, function(x) sum(x > 0) == 1))

mod_sample_type <- data.frame(
  module = rownames(df2),
  sample_type  = apply(df2, 1, function(x) {
    cl <- colnames(df2)[x > 0]
    paste(cl, collapse = ",")
  }),
  row.names = NULL
)

head(mod_sample_type)
write.csv(mod_sample_type, "module_sample_type.csv")

#5 size effect: Cramer
n_type <- sum(df2)
r_type <- nrow(df2)
c_type <- ncol(df2)

cramers_v_type <- sqrt(obs_stat_type / (n_type * (min(r_type - 1, c_type - 1))))
cramers_v_type

#Sample-type shows near-complete association with module structure (Cramér’s V = 0.72), indicating that modules are largely sample-type specific.


###Identify how many modules are pure with each Sample_type and shared modules
pure_modules <- apply(df2, 1, function(x) sum(x > 0) == 1)
pure_type <- apply(df2[pure_modules, , drop = FALSE], 1, function(x) {
  names(which(x > 0))
})
table(pure_type)

all(rownames(df3) == rownames(df2))  # should be TRUE
sum(df2[, "sample_type1"] > 0 & df2[, "sample_type2"] > 0)
sum(df3[, "sample1"])
sum((df3[, "sample1"] > 0) & (df2[, "sample_type2"] > 0))

# Count how many of the samples are present in each module
sample_count <- apply(df3[, c("sample1", "sample2", "sample3")], 1, function(x) sum(x > 0))

# Count modules shared by 2 or 3 selected samples
modules_shared <- sum(sample_count >= 2)
modules_shared
