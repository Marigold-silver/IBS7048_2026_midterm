# 4c번: 상위 5개 ASV -> 정렬(DECIPHER) -> neighbor-joining 나무(phangorn)
suppressPackageStartupMessages({ library(dada2); library(DECIPHER); library(phangorn) })

# 04a_dada2.R에서 저장한 키메라 제거 ASV 표 (행 = 시료, 열 = ASV 서열)
seqtab.nochim <- readRDS("seqtab_nochim.rds")

# 두 시료의 read 합계로 ASV 순위를 매겨 상위 5개만 고르고,
# 많은 순서대로 ASV1..ASV5로 이름 붙임
top5 <- order(colSums(seqtab.nochim), decreasing = TRUE)[1:5]
seqs <- colnames(seqtab.nochim)[top5]      # 열 이름이 곧 DNA 서열
names(seqs) <- paste0("ASV", 1:5)
print(data.frame(ASV = names(seqs), reads = colSums(seqtab.nochim)[top5],
                 length = nchar(seqs), row.names = NULL))

# DECIPHER로 다중 서열 정렬
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = FALSE)

# phangorn용으로 변환 -> 서열 간 거리 계산 -> NJ 나무 생성
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm   <- dist.ml(phang.align)   # 서열 간 maximum-likelihood 거리
tree <- NJ(dm)                 # neighbor-joining 나무 (root 없음)

# 그림(PNG)과 Newick 텍스트 파일로 저장
png("asv_tree.png", width = 700, height = 450)
plot(tree, main = "Neighbor-joining tree of the top 5 ASVs")
add.scale.bar()                # 축척 막대 = 자리당 치환 수
dev.off()
write.tree(tree, file = "top5_ASV_NJ.nwk")   # Newick 파일
cat(write.tree(tree), "\n")
