# Module-loading-on-Linux

(base) [oxo242@pegasus PERU_NHW_AF]$ cat check_bims.R

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(writexl)  # optional Excel export

# 1. Read BIM files
nhw_bim <- fread("apoe_pm1mb_nhw.bim",
                 col.names = c("chr","rs","zero_pos","pos","allele1","allele2")) %>% mutate(nhw_bim_order = row_number())

peru_bim <- fread("apoe_pm1mb_peruvian.bim",
                  col.names = c("chr","rs","zero_pos","pos","per_allele1","per_allele2")) %>% mutate(per_bim_order = row_number())

# 2. Read FRQ files
nhw_frq <- fread("nhw_out.frq", skip = 1,
                 col.names = c("CHR","SNP","A1","A2","MAF_NHW","NCHROBS")) %>% mutate(nhw_bim_order = row_number())

per_frq <- fread("per_out.frq", skip = 1,
                 col.names = c("CHR","SNP","per_A1","per_A2","MAF_PER","NCHROBS" )) %>% mutate(per_bim_order = row_number())

# 3. Merge BIM + FRQ within each dataset
nhw_bim_check <- left_join(nhw_bim, nhw_frq, by = "nhw_bim_order")
per_bim_check <- left_join(peru_bim, per_frq, by = "per_bim_order") %>%
  filter(per_allele1 %in% c("A","C","G","T"),
         per_allele2 %in% c("A","C","G","T"))

# 4. Merge NHW and Peruvian datasets
all_check <- inner_join(nhw_bim_check, per_bim_check, by = c("chr","pos")) %>%
  mutate(maf_diff = MAF_NHW - MAF_PER)

# 5. Output results to CSV
write.csv(all_check, "all_check_results.csv", row.names = FALSE)
write.csv(subset(all_check, abs(maf_diff) > 0.2), "maf_diff_gt0.2.csv", row.names = FALSE)
write.csv(subset(all_check, allele1 != per_allele1), "allele_mismatches.csv", row.names = FALSE)

# Optional: Excel export (if writexl is installed)
write_xlsx(list(
  all_results = all_check,
  maf_diff_gt0.2 = subset(all_check, abs(maf_diff) > 0.2),
  allele_mismatches = subset(all_check, allele1 != per_allele1)
), "all_check_results.xlsx")

# 6. Clean data for plotting
all_check_clean <- all_check %>% filter(!is.na(pos), !is.na(maf_diff))

# 7. Scatter plot: MAF difference vs position
p <- ggplot(all_check_clean, aes(x = pos, y = maf_diff, color = factor(chr))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "MAF Differences: NHW vs Peruvian",
       x = "Genomic Position",
       y = "MAF Difference (NHW - Peruvian)") +
  theme(legend.position = "none")

ggsave("maf_diff_plot.jpg", plot = p, width = 8, height = 5, dpi = 300)

# 8. Manhattan-style plot with cumulative position
all_check_clean <- all_check_clean %>%
  group_by(chr) %>%
  mutate(chr_len = max(pos)) %>%
  ungroup() %>%
  arrange(chr, pos) %>%
  mutate(bp_cum = pos + cumsum(lag(chr_len, default = 0)))

p_manhattan <- ggplot(all_check_clean, aes(x = bp_cum, y = maf_diff, color = factor(chr))) +
  geom_point(alpha = 0.7, size = 1.2) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(title = "Genome-wide Manhattan-style Plot of MAF Differences",
       x = "Chromosome",
       y = "MAF Difference (NHW - Peruvian)") +
  theme(legend.position = "none")

ggsave("maf_diff_manhattan_genomewide.jpg", plot = p_manhattan, width = 12, height = 5, dpi = 300)
