# 4b번: truncLen = c(120, 120)로 하면 어떻게 되는가?
# amplicon은 252-253 bp인데 120 + 120 = 240 bp이므로 forward와 reverse가
# 겹치지 못함. mergePairs()는 12 bp 이상 겹쳐야 하므로 하나도 병합되지 않음.
suppressPackageStartupMessages(library(dada2))

# --- 입력 파일 (4a와 같은 예제 데이터) ---------------------------------------
rawpath <- system.file("extdata", package = "dada2")
fnFs <- file.path(rawpath, c("sam1F.fastq.gz", "sam2F.fastq.gz"))
fnRs <- file.path(rawpath, c("sam1R.fastq.gz", "sam2R.fastq.gz"))
sample.names <- c("sam1", "sam2")

# 필터링된 read 저장 경로 (4a와 다른 폴더)
filtFs <- file.path("filtered120", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered120", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names; names(filtRs) <- sample.names

# --- 학생이 바꾼 설정으로 필터링 및 자르기 -----------------------------------
# truncLen = c(120, 120): 두 read 모두 120 bp에서 자름
# maxEE = c(2, 2): expected error가 2를 넘는 read 제거
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(120, 120),
                     maxEE = c(2, 2), maxN = 0, truncQ = 2, rm.phix = TRUE,
                     multithread = FALSE)
print(out)   # c(230,160)보다 살아남는 read는 많지만 병합은 불가능

# --- 오류율 학습 및 노이즈 제거 ----------------------------------------------
set.seed(100)
errF <- learnErrors(filtFs, multithread = FALSE, verbose = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE, verbose = FALSE)
dadaFs <- dada(filtFs, err = errF, multithread = FALSE, verbose = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE, verbose = FALSE)

# --- 병합 (기본값 minOverlap = 12) -------------------------------------------
cat("--- mergePairs with default minOverlap = 12 ---\n")
m <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
print(sapply(m, nrow))   # 시료별 병합된 서열 수: 0 과 0
