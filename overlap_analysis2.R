#Script for overlap analysis for LOPIT datasets in single xlxs file
#make changes to file_path, replicate names, 

library(readxl)
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)

# ---- Load all four replicates ----
file_path <- "TMT_replicates.xlsx"
sheets <- excel_sheets(file_path)
cat("Sheets found:", sheets, "\n")

#repA <- read_xlsx(file_path, sheet = "ReplicateA")
repB <- read_xlsx(file_path, sheet = "ReplicateB")
repC <- read_xlsx(file_path, sheet = "ReplicateC")
pilot <- read_xlsx(file_path, sheet = "Pilot")

# ---- Extract unique Accession IDs (column D) ----
#acc_A <- unique(na.omit(repA$Accession))
acc_B <- unique(na.omit(repB$Accession))
acc_C <- unique(na.omit(repC$Accession))
acc_Pilot <- unique(na.omit(pilot$Accession))

cat("\nUnique Accessions per replicate:\n")
#cat("  ReplicateA:", length(acc_A), "\n")
cat("  ReplicateB:", length(acc_B), "\n")
cat("  ReplicateC:", length(acc_C), "\n")
cat("  Pilot:     ", length(acc_Pilot), "\n")

# ---- Accession ID sets ----
acc_list <- list(
  #ReplicateA = acc_A,
  ReplicateB = acc_B,
  ReplicateC = acc_C,
  Pilot = acc_Pilot
)

# ---- Overlap statistics ----
all_acc <- unique(c(acc_B, acc_C, acc_Pilot))
cat("\nTotal unique accessions across all replicates:", length(all_acc), "\n")

core_overlap <- Reduce(intersect, acc_list)
cat("Core overlap (found in ALL 3 replicates):", length(core_overlap), "\n")

# Pairwise overlaps
pairs <- combn(names(acc_list), 2)
cat("\nPairwise overlaps:\n")
for (i in 1:ncol(pairs)) {
  n <- length(intersect(acc_list[[pairs[1, i]]], acc_list[[pairs[2, i]]]))
  cat(sprintf("  %s & %s: %d\n", pairs[1, i], pairs[2, i], n))
}

# Triple overlaps
triples <- combn(names(acc_list), 3)
cat("\nTriple overlaps:\n")
for (i in 1:ncol(triples)) {
  n <- length(Reduce(intersect, acc_list[triples[, i]]))
  cat(sprintf("  %s & %s & %s: %d\n", triples[1, i], triples[2, i], triples[3, i], n))
}

# ---- Unique to each replicate ----
cat("\nUnique to each replicate (not found in any other):\n")
for (nm in names(acc_list)) {
  others <- setdiff(names(acc_list), nm)
  other_acc <- unique(unlist(acc_list[others]))
  unique_to <- setdiff(acc_list[[nm]], other_acc)
  cat(sprintf("  %s only: %d\n", nm, length(unique_to)))
}

# ---- Venn Diagram 1: VennDiagram package (publication-quality) ----
venn.plot <- venn.diagram(
  x = acc_list,
  filename = NULL,
  fill = c( "#377EB8", "#4DAF4A", "#984EA3"),
  alpha = 0.4,
  col = "grey30",
  lwd = 1.5,
  category.names = c("Rep B", "Rep C", "Pilot"),
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.dist = c(0.22, 0.22, 0.22),
  cex = 1.0,
  fontface = "bold",
  main = "Protein Accession Overlap\nAcross TMT Replicates",
  main.cex = 1.5,
  main.fontface = "bold",
  margin = 0.1
)

png("Venn_Accession_Overlap.png", width = 2400, height = 2400, res = 300)
grid::grid.draw(venn.plot)
dev.off()
cat("\nSaved: Venn_Accession_Overlap.png\n")

# ---- Venn Diagram 2: ggVennDiagram (ggplot-based) ----
p <- ggVennDiagram(acc_list,
                   label = "count",
                   label_alpha = 0,
                   set_size = 4.5) +
  scale_fill_gradient(low = "white", high = "#377EB8") +
  labs(title = "Protein Accession Overlap Across TMT Replicates",
       fill = "Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave("Venn_Accession_ggVenn.png", p, width = 10, height = 8, dpi = 300)
cat("Saved: Venn_Accession_ggVenn.png\n")

# ---- Build overlap summary table with descriptions ----
# Combine all accessions with descriptions from each replicate
#desc_A <- repA %>% select(Accession, Description) %>% distinct() %>% mutate(ReplicateA = TRUE)
desc_B <- repB %>% select(Accession, Description) %>% distinct() %>% mutate(ReplicateB = TRUE)
desc_C <- repC %>% select(Accession, Description) %>% distinct() %>% mutate(ReplicateC = TRUE)
desc_P <- pilot %>% select(Accession, Description) %>% distinct() %>% mutate(Pilot = TRUE)

# Merge all by Accession, keep first non-NA description
overlap_table <- tibble(Accession = all_acc) %>%
  #left_join(desc_A %>% select(Accession, ReplicateA), by = "Accession") %>%
  left_join(desc_B %>% select(Accession, ReplicateB), by = "Accession") %>%
  left_join(desc_C %>% select(Accession, ReplicateC), by = "Accession") %>%
  left_join(desc_P %>% select(Accession, Pilot), by = "Accession") %>%
  mutate(across(c(ReplicateA, ReplicateB, ReplicateC, Pilot), ~replace_na(., FALSE))) %>%
  mutate(n_replicates = ReplicateA + ReplicateB + ReplicateC + Pilot)

# Add description (take from first replicate that has it)
all_desc <- bind_rows(
  #desc_A %>% select(Accession, Description),
  desc_B %>% select(Accession, Description),
  desc_C %>% select(Accession, Description),
  desc_P %>% select(Accession, Description)
) %>%
  filter(!is.na(Description)) %>%
  distinct(Accession, .keep_all = TRUE)

overlap_table <- overlap_table %>%
  left_join(all_desc, by = "Accession") %>%
  arrange(desc(n_replicates), Accession)

# Save the full overlap table
write_csv(overlap_table, "overlap_summary_table.csv")
cat("Saved: overlap_summary_table.csv\n")

# Print summary
cat("\n---- Overlap Summary ----\n")
cat(sprintf("Found in 3/3 replicates: %d proteins\n", sum(overlap_table$n_replicates == 4)))
cat(sprintf("Found in 2/3 replicates: %d proteins\n", sum(overlap_table$n_replicates == 3)))
cat(sprintf("Found in 1/3 replicates: %d proteins\n", sum(overlap_table$n_replicates == 2)))
#cat(sprintf("Found in 0/3 replicates: %d proteins\n", sum(overlap_table$n_replicates == 1)))

# Print first few core overlap proteins with descriptions
cat("\n---- Top 20 Core Overlap Proteins (in all 4 replicates) ----\n")
core_table <- overlap_table %>% filter(n_replicates == 3) %>% head(20)
for (i in 1:nrow(core_table)) {
  cat(sprintf("  %s | %s\n", core_table$Accession[i],
              str_trunc(as.character(core_table$Description[i]), 80)))
}

cat("\nDone!\n")
