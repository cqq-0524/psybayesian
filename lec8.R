# lec8.R
# 与 lec8.qmd 对齐的 R 脚本（qmd 是幻灯片代码的唯一真源）。
# 每个代码块按 qmd 中出现的顺序排列；学生练习的骨架块以注释形式保留。

# ---- qmd chunk: setup (line 27) ----


# 数值展示：避免 knitr 把 inline 数值自动转成 LaTeX 科学计数法（会丢掉小数位）
fmt_bf <- function(x) {
  if (abs(x) >= 1e-2 && abs(x) < 1e6) {
    formatC(x, format = "f", digits = 2)
  } else {
    formatC(x, format = "e", digits = 2)
  }
}

# 区间文本，如 [0.76, 0.89]
hdi_txt <- function(x) paste0("[", fmt_bf(x$CI_low), ", ", fmt_bf(x$CI_high), "]")

# ---- qmd chunk: packages (line 209) ----

# 安装和加载包
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace('pacman', quietly = TRUE)) {
    install.packages('pacman')
}
pacman::p_load("tidyverse","ggplot2", "dplyr","gridExtra","papaja", "patchwork","bayesplot","rstan","bridgesampling","bayestestR") # "bayesrules"

# ---- qmd chunk: cache-helper (line 222) ----

# 采样结果缓存：tmpdata/ 下已有同名 .rds 就直接读取；换个文件名或删掉即可重跑
cache_rds <- function(path, code) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(path)) {
    readRDS(path)
  } else {
    value <- code
    saveRDS(value, path)
    value
  }
}

# ---- qmd chunk: bb-fit (line 241) ----

# 准备数据
stan_data <- list(N = 118,Y = 99)

# 定义模型
stan_model_code <- "
data {
  int<lower=0> N;       
  int<lower=0, upper=N> Y; 
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  theta ~ beta(2, 2);
  Y ~ binomial(N, theta);
}
"

# 进行采样（结果缓存到 tmpdata/，重复渲染直接读取；改模型后删掉该 .rds 即可重跑）
fit_bb <- cache_rds("tmpdata/lec8_fit_bb.rds",
  stan(
    model_code = stan_model_code,      # 定义的模型或模型文件路径
    data = stan_data,                  # 输入数据
    chains = 4,                        # 马尔可夫链数量
    iter = 5000,                       # 总迭代次数（每个链）
    warmup = 1000,                     # 热身迭代次数（不保存）
    seed = 84735
  )
)

# ---- qmd chunk: plot-posterior-areas (line 278) ----

pp <- bayesplot::mcmc_areas(fit_bb, regex_pars = "theta",  prob = 0.95) +
 labs(
   title = "Posterior distributions",
   subtitle = "with medians and 95% intervals"
 ) +
 papaja::theme_apa()

 pp

# ---- qmd chunk: posterior-samples (line 338) ----

# 提取后验样本（后续各节统一复用 theta_post）
theta_post <- as.numeric(rstan::extract(fit_bb)$theta)
theta_post_df <- data.frame(theta = theta_post)

# 95% HDI（后续文字与摘要页复用）
hdi_post <- bayestestR::hdi(theta_post_df, ci = 0.95)
hdi_post  # 可自行更改不同范围的 HDI 区间

# ---- qmd chunk: rope-call (line 404) ----

# get ROPE
bayestestR::rope(theta_post_df, range = c(0.4, 0.6))

# ---- qmd chunk: plot-rope (line 425) ----

result <- bayestestR::rope(theta_post_df, range = c(0.4, 0.6))

plot(result, rope_color = "grey70") + papaja::theme_apa()

# ---- qmd chunk: plot-six-tests (line 474) ----

# ============================================
# 1. 后验样本（直接复用前面提取的 theta_post）
# ============================================

# ============================================
# 2. 定义假设检验函数（基于 theta，ROPE = [0.4, 0.6]）
# ============================================

# 区间假设：ROPE 检验与等价性检验共用（null = 区间外，interest = 区间内）
region_test <- function(samples, rope_lower = 0.4, rope_upper = 0.6,
                        name = "ROPE Test (Practical Equivalence)") {
  interest_region <- samples >= rope_lower & samples <= rope_upper  # inside ROPE
  null_region <- !interest_region                                    # outside ROPE
  
  list(
    name = name,
    null = paste0("θ < ", rope_lower, " or θ > ", rope_upper, " (outside ROPE)"),
    interest = paste0(rope_lower, " ≤ θ ≤ ", rope_upper, " (inside ROPE)"),
    p_null = sum(null_region) / length(samples),
    p_interest = sum(interest_region) / length(samples),
    rope_lower = rope_lower,
    rope_upper = rope_upper
  )
}

# 单侧 NHST: H0: θ ≤ 0.5 vs H1: θ > 0.5
one_sided_nhst <- function(samples, null_val = 0.5) {
  null_region <- samples <= null_val
  interest_region <- samples > null_val
  
  list(
    name = "One-sided NHST",
    null = paste0("θ ≤ ", null_val),
    interest = paste0("θ > ", null_val),
    p_null = sum(null_region) / length(samples),
    p_interest = sum(interest_region) / length(samples),
    null_val = null_val
  )
}

# 双侧最小效应检验: H0: θ ∈ [0.4, 0.6] (no meaningful effect)
two_sided_met <- function(samples, rope_lower = 0.4, rope_upper = 0.6) {
  null_region <- samples >= rope_lower & samples <= rope_upper
  interest_region <- !null_region
  
  list(
    name = "Two-sided Minimal Effect Test",
    null = paste0(rope_lower, " ≤ θ ≤ ", rope_upper, " (trivial effect)"),
    interest = paste0("θ < ", rope_lower, " or θ > ", rope_upper, " (meaningful effect)"),
    p_null = sum(null_region) / length(samples),
    p_interest = sum(interest_region) / length(samples),
    rope_lower = rope_lower,
    rope_upper = rope_upper
  )
}

# 单侧最小效应检验: H0: θ ≤ 0.6 vs H1: θ > 0.6
one_sided_met <- function(samples, threshold = 0.6) {
  null_region <- samples <= threshold
  interest_region <- samples > threshold
  
  list(
    name = "One-sided Minimal Effect Test",
    null = paste0("θ ≤ ", threshold),
    interest = paste0("θ > ", threshold),
    p_null = sum(null_region) / length(samples),
    p_interest = sum(interest_region) / length(samples),
    threshold = threshold
  )
}

# 等价性检验：与 ROPE 检验共用同一套区间逻辑，仅名称不同
equivalence_test <- function(samples, rope_lower = 0.4, rope_upper = 0.6) {
  region_test(samples, rope_lower, rope_upper, name = "Equivalence Test")
}

# 非劣性检验: H0: θ < 0.4 vs H1: θ ≥ 0.4
non_inferiority_test <- function(samples, threshold = 0.4) {
  null_region <- samples < threshold
  interest_region <- samples >= threshold
  
  list(
    name = "Non-inferiority Test",
    null = paste0("θ < ", threshold),
    interest = paste0("θ ≥ ", threshold),
    p_null = sum(null_region) / length(samples),
    p_interest = sum(interest_region) / length(samples),
    threshold = threshold
  )
}

# ============================================
# 3. 执行所有检验 
# ============================================

results <- list(
  region_test(theta_post, rope_lower = 0.4, rope_upper = 0.6),
  one_sided_nhst(theta_post, null_val = 0.5),
  two_sided_met(theta_post, rope_lower = 0.4, rope_upper = 0.6),
  one_sided_met(theta_post, threshold = 0.6),
  equivalence_test(theta_post, rope_lower = 0.4, rope_upper = 0.6),
  non_inferiority_test(theta_post, threshold = 0.4)
)

# ============================================
# 4. 可视化函数
# ============================================

visualize_test_theta <- function(samples, test_result, x_range = c(0, 1)) {
  df <- data.frame(theta = samples)

  # 所有面板统一：阴影 = 零假设区间（null region）
  shade_null <- function(p, lower, upper) {
    p + annotate("rect", xmin = lower, xmax = upper, ymin = 0, ymax = Inf,
                 alpha = 0.18, fill = "steelblue")
  }

  p <- ggplot(df, aes(x = theta)) +
    geom_density(fill = "gray80", alpha = 0.5) +
    labs(
      title = test_result$name,
      subtitle = paste0(
        "Null: ", test_result$null, "\n",
        "Interest: ", test_result$interest, "\n",
        "p(null) = ", sprintf("%.3f", test_result$p_null),
        " | p(interest) = ", sprintf("%.3f", test_result$p_interest)
      ),
      x = "Parameter θ",
      y = "Posterior Density"
    ) +
    papaja::theme_apa() +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7),
      axis.title = element_text(size = 9)
    ) +
    xlim(x_range)

  # 阴影与虚线都标出零假设区间
  if (test_result$name %in% c("ROPE Test (Practical Equivalence)", "Equivalence Test")) {
    p <- shade_null(p, -Inf, test_result$rope_lower)
    p <- shade_null(p, test_result$rope_upper, Inf)
    p <- p + geom_vline(xintercept = c(test_result$rope_lower, test_result$rope_upper),
                        linetype = "dotted", color = "steelblue", linewidth = 0.8)
  } else if (test_result$name == "Two-sided Minimal Effect Test") {
    p <- shade_null(p, test_result$rope_lower, test_result$rope_upper)
    p <- p + geom_vline(xintercept = c(test_result$rope_lower, test_result$rope_upper),
                        linetype = "dotted", color = "steelblue", linewidth = 0.8)
  } else if (test_result$name == "One-sided NHST") {
    p <- shade_null(p, -Inf, test_result$null_val)
    p <- p + geom_vline(xintercept = test_result$null_val, linetype = "dotted",
                        color = "steelblue", linewidth = 0.8)
  } else {
    # One-sided Minimal Effect Test / Non-inferiority Test
    p <- shade_null(p, -Inf, test_result$threshold)
    p <- p + geom_vline(xintercept = test_result$threshold, linetype = "dotted",
                        color = "steelblue", linewidth = 0.8)
  }

  p
}

# 生成图形
plots <- lapply(results, function(r) visualize_test_theta(theta_post, r))

# 组合图形（2列）
do.call(grid.arrange, c(plots, ncol = 2))

# ---- qmd chunk: rope-print (line 649) ----

# ============================================
# 5. ROPE、HDI、后验摘要
# ============================================

# 使用 bayestestR 的 rope() 函数验证（ROPE = [0.4, 0.6]）
rope_bayestestR <- bayestestR::rope(theta_post, range = c(0.4, 0.6), ci = 0.95)
cat("\n=== bayestestR ROPE Analysis (95% CI) ===\n")
print(rope_bayestestR)

# ---- qmd chunk: hdi-summary (line 665) ----

# 95% HDI（复用前文已算好的 hdi_post）
cat("\n95% Highest Density Interval (HDI):", 
    round(hdi_post$CI_low, 4), "to", round(hdi_post$CI_high, 4), "\n")

# 后验均值与中位数
cat("Posterior Mean of θ:", round(mean(theta_post), 4), "\n")
cat("Posterior Median of θ:", round(median(theta_post), 4), "\n")

# 判断是否等效：如果 HDI 完全在 [0.4, 0.6] 内，则 strong evidence for equivalence
hdi_in_rope <- hdi_post$CI_low > 0.4 & hdi_post$CI_high < 0.6
cat("HDI fully within ROPE [0.4, 0.6]?", ifelse(hdi_in_rope, "YES", "NO"), "\n")

# ---- qmd chunk: print-fit (line 816) ----

print(fit_bb, pars = c("theta"))

# ---- qmd chunk: prior-fit (line 846) ----

# =====================================
# 从先验采样（无数据似然）
# =====================================
# 定义模型
stan_prior_code <- "
data {
  int<lower=0> N;       
  int<lower=0, upper=N> Y; 
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  theta ~ beta(2, 2);
  
  // 注释掉似然
  // Y ~ binomial(N, theta);
}
"

# 编译先验模型
prior_model <- stan_model(model_code = stan_prior_code)

# 从先验采样（使用虚拟数据；结果同样缓存到 tmpdata/）
fit_prior <- cache_rds("tmpdata/lec8_fit_prior.rds",
  sampling(
    prior_model,
    data = stan_data,  # 虚拟
    iter = 5000,
    warmup = 1000,
    chains = 4,
    seed = 84735
  )
)

# 提取先验样本（后续 3.1 节复用）
theta_prior <- as.numeric(rstan::extract(fit_prior)$theta)

# ---- qmd chunk: bf-kde (line 909) ----

# Savage-Dickey：BF01 = 后验在 θ = 0.5 处的密度 / 先验在 θ = 0.5 处的密度
# 本例为共轭模型：先验 Beta(2, 2)，后验 Beta(2 + 99, 2 + 118 - 99) = Beta(101, 21)
# 直接用解析密度；核密度估计在 1e-12 量级的密度处分辨不出（会得到 0）
prior_at_null     <- dbeta(0.5, shape1 = 2,   shape2 = 2)
posterior_at_null <- dbeta(0.5, shape1 = 101, shape2 = 21)

# ---- qmd chunk: str-posterior (line 919) ----

str(theta_post)

# ---- qmd chunk: str-prior (line 925) ----

str(theta_prior)

# ---- qmd chunk: bf-savage-dickey (line 934) ----

# 计算 Savage-Dickey 比
BF_01 <- posterior_at_null / prior_at_null  # 支持 H0 相对于 H1
BF_10 <- 1 / BF_01  # 支持 H1（θ ≠ 0.5）相对于 H0（θ = 0.5）

cat("The BF_10 is", fmt_bf(BF_10), "\n")
cat("The BF_01 is", fmt_bf(BF_01), "\n")

# ---- qmd chunk: plot-prior-posterior (line 948) ----

# 绘图用的核密度：先验与后验在同一网格上各估一次
grid_x <- range(c(theta_prior, theta_post))
dens_prior <- density(theta_prior, from = grid_x[1], to = grid_x[2], n = 1024)
dens_post  <- density(theta_post,  from = grid_x[1], to = grid_x[2], n = 1024)

# 准备绘图数据
df_prior <- data.frame(
  theta = dens_prior$x,
  density = dens_prior$y,
  type = "Prior"
)

df_posterior <- data.frame(
  theta = dens_post$x,
  density = dens_post$y,
  type = "Posterior"
)

df_combined <- bind_rows(df_prior, df_posterior)

# 主密度图
main_plot <- ggplot(df_combined, aes(x = theta, y = density, color = type)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("Prior" = "steelblue", "Posterior" = "darkorange")) +
  labs(
    title = "Prior and Posterior Distribution Comparison",
    x = expression(theta),
    y = "Density",
    color = "Distribution"
  ) +
  papaja::theme_apa() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

# 分面版本：拆成两个独立面板，保证每个面板都有自己的 x 轴
panel_prior <- ggplot(df_prior, aes(x = theta, y = density)) +
  geom_area(fill = "steelblue", alpha = 0.6) +
  labs(title = "Prior", x = expression(theta), y = "Density") +
  xlim(0, 1) +
  papaja::theme_apa()

panel_posterior <- ggplot(df_posterior, aes(x = theta, y = density)) +
  geom_area(fill = "darkorange", alpha = 0.6) +
  labs(title = "Posterior", x = expression(theta), y = "Density") +
  xlim(0, 1) +
  papaja::theme_apa()

# 组合图形
combined_plot <- main_plot / (panel_prior | panel_posterior) +
  plot_annotation(title = "Bayesian Analysis: Prior vs Posterior Distributions",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

print(combined_plot)

# ---- qmd chunk: plot-jasp-style (line 1038) ----

jasp_samples <- bind_rows(
  data.frame(theta = theta_prior, type = "Prior"),
  data.frame(theta = theta_post,  type = "Posterior")
)

ggplot(jasp_samples, aes(x = theta, fill = type)) +
  geom_density(alpha = 0.6, colour = NA) +
  scale_fill_manual(values = c(Prior = "grey70", Posterior = "steelblue")) +
  annotate("text", x = 0.02, y = Inf, hjust = 0, vjust = 1.4,
           label = paste0("BF10 = ", fmt_bf(BF_10), "\nBF01 = ", fmt_bf(BF_01)),
           size = 4) +
  labs(x = expression(theta), y = "Density", fill = NULL) +
  xlim(0, 1) +
  papaja::theme_apa() +
  theme(legend.position = "top")

# ---- qmd chunk: bf-bridge (line 1064) ----

# H1：用后验样本 + 解析 log 后验密度做桥采样
log_post <- function(theta, data) {
  dbeta(theta, 2, 2, log = TRUE) + dbinom(data$Y, data$N, theta, log = TRUE)
}

samples_H1 <- matrix(theta_post, ncol = 1, dimnames = list(NULL, "theta"))

bridge_H1 <- bridgesampling::bridge_sampler(
  samples_H1,
  log_posterior = log_post,
  data = list(N = 118, Y = 99),
  lb = c(theta = 0),
  ub = c(theta = 1),
  silent = TRUE
)

# H0：θ 固定在 0.5，边际似然就是该点的似然
logml_H0 <- dbinom(99, 118, 0.5, log = TRUE)

BF01_bridge <- exp(logml_H0 - bridge_H1$logml)
BF10_bridge <- 1 / BF01_bridge

cat("BF_01 (bridge sampling) =", fmt_bf(BF01_bridge), "\n")
cat("BF_10 (bridge sampling) =", fmt_bf(BF10_bridge), "\n")

# ---- qmd chunk: bf-rope (line 1153) ----

# 区间假设：H0: θ ∈ [0.4, 0.6]（实用等效），H1: θ ∉ [0.4, 0.6]
p_rope_post  <- pbeta(0.6, 101, 21) - pbeta(0.4, 101, 21)   # 后验落在 ROPE 内
p_rope_prior <- pbeta(0.6, 2, 2)    - pbeta(0.4, 2, 2)      # 先验落在 ROPE 内

BF01_rope <- (p_rope_post / (1 - p_rope_post)) *
             ((1 - p_rope_prior) / p_rope_prior)
BF10_rope <- 1 / BF01_rope

cat("P(θ ∈ ROPE | Y)      =", fmt_bf(p_rope_post), "\n")
cat("P(θ ∈ ROPE)  (prior) =", fmt_bf(p_rope_prior), "\n")
cat("BF_01 =", fmt_bf(BF01_rope), " | BF_10 =", fmt_bf(BF10_rope), "\n")

# ---- qmd chunk: exercise-step1 (line 1200) ----

# [练习骨架，qmd 中为 eval: false]
#第一步：后验样本提取及对应转换

# 用 as.vector() 统一格式
# theta_vec <- as.vector(...)  
# 查看总样本数（确保所有有效迭代都被包含）
# total_samples <- length(theta_vec)  
# cat("theta的后验总样本数：", total_samples, "\n")

# ---- qmd chunk: exercise-step2 (line 1214) ----

# [练习骨架，qmd 中为 eval: false]
#第二步： 计算ROPE区间内的后验概率比

# 定义 ROPE 区间
# rope_low <- ...
# rope_high <- ...

# 筛选：判断每个样本是否落在 (0.4, 0.6) 内
# theta_in_rope <- (... > rope_low) & (... < rope_high)  

# 计算 ROPE 内的后验概率比例
# theta_rope_prob <- mean(theta_in_rope)  

# 输出结果
# cat("参数 theta 落在 ROPE (", rope_low, ",", rope_high, ") 内的后验概率：", round(theta_rope_prob, 3), "\n")
# 可选：输出落在 ROPE 内的样本数
#cat("落在 ROPE 内的样本数：", sum(theta_in_rope), "\n")

# ---- qmd chunk: exercise-step3 (line 1237) ----

# [练习骨架，qmd 中为 eval: false]
# 第三步： 计算贝叶斯因子：BF₀₁ = 支持 H₀（theta在ROPE内）的概率 / 支持 H₁（theta在ROPE外）的概率
# bf01 <- (...) / (...)  
# cat("贝叶斯因子 BF₀₁（支持 theta 在 (0.4,0.6) 内）：", round(bf01, 3), "\n")

# ---- qmd chunk: helper-odds (line 1259) ----

# 请直接运行，无需修改
# 计算「区间内 / 区间外」的概率比（练习2用；region 默认取 ROPE = [0.4, 0.6]）
calculate_odds <- function(theta_samples, region = c(0.4, 0.6)) {
  # 区间内的样本
  in_range <- theta_samples[theta_samples >= region[1] & theta_samples <= region[2]]

  # 区间外的样本
  out_of_range <- theta_samples[theta_samples < region[1] | theta_samples > region[2]]

  # 区间内外的比例
  P_in_range <- length(in_range) / length(theta_samples)
  P_out_of_range <- length(out_of_range) / length(theta_samples)

  # 概率比（区间外 / 区间内）
  P_out_of_range / P_in_range
}

# ---- qmd chunk: exercise-directional-samples (line 1283) ----

# [练习骨架，qmd 中为 eval: false]
#=====================================
#      基于方向的贝叶斯因子计算
#      自行练习
#=====================================

# 获取theta的先验分布的采样
#theta_prior =  


# 获取theta的后验分布的采样
# theta_post = ...

# ---- qmd chunk: exercise-directional-bf (line 1300) ----

# [练习骨架，qmd 中为 eval: false]
# 定义区间
# region = ...

# 计算先验比
# prior_odds = ...

# 计算后验比
# posterior_odds = ...

# 计算贝叶斯因子
# BF_10 = ...

# cat("贝叶斯因子（BF_10）:", BF_10,"\n")

# ---- qmd chunk: exercise-directional-plot (line 1320) ----

# [练习骨架，qmd 中为 eval: false]
# df_plot <- tibble(
#   prior = ...,
#   posterior = ...
# ) %>%
#   pivot_longer(
#     cols = everything(),
#     names_to = "distribution",
#     values_to = "theta"
#   )

# === 方向性假设：H1: theta < 0.5 ===
# h0_boundary <- 0.5  # H₀: theta >= 0.5

# 计算最大密度用于标签定位
# dens_prior <- density(theta_prior, from = 0, to = 1)
# dens_post <- density(theta_post, from = 0, to = 1)
# max_dens <- max(dens_prior$y, dens_post$y)

# 绘图
# ggplot(df_plot, aes(x = theta, fill = distribution, color = distribution)) +
#   geom_density(alpha = 0.3, linewidth = 1) +


#   annotate("rect", 
#            xmin = h0_boundary, xmax = Inf,
#            ymin = -Inf, ymax = Inf,
#            fill = "blue", alpha = 0.2, inherit.aes = FALSE) +

  # 垂直线在 0 处（假设边界）
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +

  # 标签：H₀ 在右侧，H₁ 在左侧
#   annotate("text",
#            x = 0.75,  # H₀ 区域（右侧）
#            y = max_dens * 0.8,
#            label = "H0: θ ≥ 0.5",
#            color = "blue",
#            fontface = "bold",
#            size = 4,
#            parse = FALSE) +

#   annotate("text",
#            x = 0.25,  # H₁ 区域（左侧）
#            y = max_dens * 0.6,
#            label = "H1: θ < 0.5",
#            color = "darkorange",
#            fontface = "bold",
#            size = 4,
#            parse = FALSE) +

#   xlim(0, 1) +

#   labs(
#     title = "Directional Hypothesis: Testing for Negative Effect",
#     x = expression(theta),
#     y = "Density",
#     fill = "Distribution",
#     color = "Distribution"
#   ) +

#   scale_fill_manual(values = c("prior" = "cyan", "posterior" = "coral")) +
#   scale_color_manual(values = c("prior" = "cyan", "posterior" = "coral")) +

#   theme_minimal() +
#   theme(
#     panel.border = element_rect(fill = NA, colour = "black"),
#     panel.grid.minor = element_blank(),
#     legend.position = "upper right"
#   )
