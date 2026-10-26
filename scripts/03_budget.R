# Problem 3
sequencing <- data.frame(
  platform = c("Illumina", "Oxford Nanopore", "PacBio HiFi"),
  yield_Gb = c(18, 30, 12),
  host_fraction = c(0.20, 0.70, 0.10),
  cost_USD = c(150, 200, 300)
)
budget <- 600
sequencing$num_samples <- floor(budget/sequencing$cost_USD)
sequencing$total_microbial_Gb <- sequencing$num_samples * sequencing$yield_Gb * (1 - sequencing$host_fraction)
print(sequencing)
