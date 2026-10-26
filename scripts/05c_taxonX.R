# 5c번: microbial read 수, Taxon_X 제외 후 남은 read 수, fold change
# (값은 모두 가상; Taxon_X 비율은 microbial read 중의 비율)
condition <- c("Untreated", "Treated")
total_reads <- c(20000000, 20000000)
host_fraction <- c(0.99, 0.90)
taxon_X_fraction <- c(0.10, 0.91)

# microbial read = 전체 read x (1 - host 비율)
microbial_reads <- total_reads * (1 - host_fraction)
# Taxon_X 제외 후 남은 read = microbial read x (1 - Taxon_X 비율)
non_X_reads <- microbial_reads * (1 - taxon_X_fraction)

res <- data.frame(condition, microbial_reads, non_X_reads)
print(res)

# fold change = Treated / Untreated (두 번째 값 / 첫 번째 값)
cat("Fold change (Treated/Untreated), microbial reads:",
    microbial_reads[2] / microbial_reads[1], "\n")      # 10
cat("Fold change (Treated/Untreated), non-Taxon_X reads:",
    non_X_reads[2] / non_X_reads[1], "\n")              # 1
