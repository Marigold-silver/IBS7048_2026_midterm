# 6b번: score가 0.20에서 0.60으로 증가할 때의 효과
score <- c(0.20, 0.60)
OR   <- c(estimate = 0.47, lower = 0.22, upper = 0.97)       # score 1 단위 증가당
CASI <- c(estimate = -0.92, lower = -1.51, upper = -0.33)    # score 1 단위 증가당
p0 <- 0.40                                                   # 가정한 초기 확률

delta <- diff(score)     # score 변화량 = 0.6 - 0.2 = 0.4

# log-odds가 score에 선형이므로 log(OR)도 delta배가 됨 -> OR^delta
# (신뢰구간 하한/상한도 같은 방식으로 변환)
OR_delta <- OR^delta
# CASI는 score에 선형 -> 계수 x delta (신뢰구간도 동일)
CASI_delta <- CASI * delta
print(round(OR_delta, 3)); print(round(CASI_delta, 3))

# 초기 확률이 p0인 아이의 새 확률:
# 확률 -> odds, odds ratio를 곱함, odds -> 확률
odds0 <- p0 / (1 - p0)
odds1 <- odds0 * OR_delta
p1 <- odds1 / (1 + odds1)
print(round(p1, 3))
cat("Change in probability:", round(p1["estimate"] - p0, 3), "\n")
