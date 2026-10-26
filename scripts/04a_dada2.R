# 4a번: 예제 시료 2개(sam1, sam2)에 대한 DADA2 파이프라인
# 결과: 품질 그래프, read 추적표, 키메라 제거된 ASV 표
# 필요 패키지(README 참고): dada2, ggplot2

library(dada2)      # 노이즈 제거(denoising), 병합, 키메라 제거
library(ggplot2)    # 품질 그래프 저장용 ggsave()

# Bioconductor 병렬 처리를 한 프로세스에서만 실행.
# Windows 기본값은 R 작업 프로세스를 따로 띄워서 멈추거나 오류가 날 수 있음.
BiocParallel::register(BiocParallel::SerialParam())
# --- 입력 파일 ---------------------------------------------------------------
# dada2 패키지에 들어 있는 예제 read. F = forward(R1), R = reverse(R2)
rawpath <- system.file("extdata", package = "dada2")
fnFs <- file.path(rawpath, c("sam1F.fastq.gz", "sam2F.fastq.gz"))
fnRs <- file.path(rawpath, c("sam1R.fastq.gz", "sam2R.fastq.gz"))
sample.names <- c("sam1", "sam2")

# --- 1. 품질 프로파일 --------------------------------------------------------
# 사이클(x)별 quality score(y). 중앙값(초록 선)이 떨어지는 지점 앞에서
# read를 잘라야 함. 보통 reverse read가 더 일찍 떨어짐.
# 그래프 객체를 만든 뒤 ggsave로 PNG로 저장.
pF <- plotQualityProfile(fnFs)
ggsave("quality_forward.png", plot = pF, width = 9, height = 5, dpi = 100)
pR <- plotQualityProfile(fnRs)
ggsave("quality_reverse.png", plot = pR, width = 9, height = 5, dpi = 100)

# --- 2. 필터링 및 자르기 -----------------------------------------------------
# 필터링된 read를 저장할 경로("filtered" 폴더는 자동 생성됨)
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names; names(filtRs) <- sample.names
# truncLen = c(230, 160): forward는 230 bp, reverse는 160 bp에서 자름
#                         (이보다 짧은 read는 버려짐)
# maxEE = c(2, 2): expected error가 2를 넘는 read 제거
# maxN = 0: N이 하나라도 있는 read 제거;  truncQ = 2: Q <= 2인 첫 염기에서 자름
# rm.phix = TRUE: PhiX(대조군) read 제거
# out = 시료별 입력 read 수 / 출력 read 수 표
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(230, 160), maxEE = c(2, 2),
                     maxN = 0, truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = FALSE)
print(out)

# --- 3. 오류율 학습, 노이즈 제거 ---------------------------------------------
set.seed(100)                                    # 결과 재현을 위한 시드
# 데이터에서 염기 치환(A->C, A->G, ...)별 오류율을 학습
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
# 오류 모델로 오류가 낀 read를 실제 서열로 합쳐 ASV를 추정
# (forward와 reverse를 따로 수행)
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

# --- 4. paired read 병합 -----------------------------------------------------
# forward와 reverse가 겹치는 구간(최소 12 bp, 기본 설정에서 불일치 없음)으로
# 이어 붙임. 겹치지 않거나 불일치하는 쌍은 버려짐.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# --- 5. 서열 표 만들기, 키메라 제거 -------------------------------------------
# 행 = 시료, 열 = ASV 서열, 값 = read 수
seqtab <- makeSequenceTable(mergers)
cat("ASVs before chimera removal:", ncol(seqtab), "\n")
# 키메라 = PCR 중 두 모서열이 이어져 생긴 인공 서열
# method = "consensus": 여러 시료에서 일관되게 키메라로 판정된 ASV 제거
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = FALSE, verbose = TRUE)
cat("ASVs after chimera removal:", ncol(seqtab.nochim), "\n")
cat("Reads before/after chimera removal:", sum(seqtab), sum(seqtab.nochim), "\n")
cat("Chimeric ASVs removed:", ncol(seqtab) - ncol(seqtab.nochim), "\n")
cat("Chimeric reads removed:", sum(seqtab) - sum(seqtab.nochim),
    sprintf("(%.1f%% of merged reads)\n", 100 * (1 - sum(seqtab.nochim) / sum(seqtab))))

# --- 6. 단계별 read 추적표 ---------------------------------------------------
# getN: dada / merge 결과의 전체 read 수
getN <- function(x) sum(getUniques(x))
# 행 = 시료, 열 = 단계
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(track)
print(round(100 * track / track[, "input"], 1))   # 입력 대비 남은 read 비율(%)

# 키메라 제거된 표(04c_tree.R에서 사용)와 추적표 저장
saveRDS(seqtab.nochim, "seqtab_nochim.rds")
write.csv(track, "read_tracking.csv")
