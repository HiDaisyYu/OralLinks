############################################################
# full_mr_loop.R
# Industrial-scale Mendelian Randomization (MR) pipeline
# --------------------------------------------------------
# Steps:
# 1. 设置 API token (服务器环境用 Sys.setenv)
# 2. 定义口腔暴露关键词 & 系统性疾病关键词字典
# 3. 从 OpenGWAS 拉取 GWAS 元数据，筛选暴露和结局集合
# 4. 对每个暴露提取工具变量 (p < 5e-8, clumping)
# 5. 对每个结局提取 SNP 数据，harmonise，跑 MR
# 6. 保存所有结果，筛选显著结果 (FDR < 0.05)
# 7. 生成下载清单 (暴露+结局的 ID 和 Trait 名称)
############################################################

# -----------------------------
# 0. 设置 API token (服务器环境)
# -----------------------------
Sys.setenv(IEU_TOKEN = "your_api_token_here")

# -----------------------------
# 1. 加载必要的 R 包
# -----------------------------
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(stringr)

# -----------------------------
# 2. 定义关键词
# -----------------------------
# Oral (暴露)关键词
oral_keys <- c(
  "periodont","gingiv","gum","tooth","teeth","dental","caries","decay",
  "mouth","oral","tongue","saliv","xerostomia","denture","edentul",
  "leukoplakia","canker","stomatitis"
)

# Systemic (结局)关键词字典
systemic_dict <- list(
  Cognitive = c("alzheimer","dementia","parkinson","cognitive","memory","schizophrenia",
                "depress","anxiety","bipolar","insomnia","neuroticism"),
  Addiction = c("smoking","tobacco","cigarette","nicotine","alcohol","drink","liquor",
                "cannabis","marijuana","weed","dependence"),
  Autoimmune = c("sjogren","lupus","sle","rheumatoid","arthritis","psoriasis",
                 "autoimmune","type 1 diabet","multiple sclerosis","celiac"),
  Cardio = c("coronary","myocardial","stroke","heart failure","atrial fib",
             "hypertension","blood pressure"),
  Metabolic = c("type 2 diabet","bmi","obesity","lipid","cholesterol","hdl","ldl","triglyceride"),
  Digestive = c("bowel","colitis","crohn","ibd","gastric","stomach","esophag",
                "reflux","gerd","peptic","pancrea","liver","gut"),
  Head_Neck_Cancer = c("head and neck","oral cancer","mouth cancer","tongue cancer",
                       "esophageal cancer","nasopharyngeal","larynx","pharynx")
)

# -----------------------------
# 3. 拉取 GWAS 元数据并筛选
# -----------------------------
message("Fetching GWAS metadata from OpenGWAS...")
all_gwas <- as.data.table(available_outcomes())

# 只保留欧洲人群 + 样本量 > 10000
all_gwas <- all_gwas[population == "European" & sample_size > 10000]

# 暴露集合：口腔关键词
is_oral <- function(trait) {
  t <- tolower(trait)
  any(str_detect(t, str_c(oral_keys, collapse="|")))
}
exposures <- all_gwas[sapply(trait, is_oral)]
exposures <- unique(exposures[, .(id, trait, sample_size, ncase, author, year)])

# 结局集合：系统性关键词
label_system <- function(trait, dict) {
  t <- tolower(trait)
  for (sys in names(dict)) {
    if (any(str_detect(t, str_c(dict[[sys]], collapse="|")))) return(sys)
  }
  return(NA_character_)
}
outcomes <- all_gwas[sapply(trait, function(x) !is.na(label_system(x, systemic_dict)))]
outcomes$System <- sapply(outcomes$trait, label_system, dict = systemic_dict)
outcomes <- unique(outcomes[, .(id, trait, System, sample_size, ncase, author, year)])

message(sprintf("Exposures (oral): %d | Outcomes (systemic): %d",
                nrow(exposures), nrow(outcomes)))

# -----------------------------
# 4. MR 参数
# -----------------------------
instrument_p_threshold <- 5e-8   # 工具变量筛选阈值
mr_result_p_threshold  <- 0.05   # MR 结果显著性阈值
clump_r2    <- 0.001
clump_kb    <- 10000
pop_ref     <- "EUR"
min_instruments <- 3  # 最少工具变量数

# -----------------------------
# 5. MR 循环
# -----------------------------
results_all <- list()

for (i in 1:nrow(exposures)) {
  exp_id   <- exposures$id[i]
  exp_name <- exposures$trait[i]
  message(sprintf("\n[Exposure] %s (%s)", exp_name, exp_id))

  # 提取暴露工具变量
  exp_iv <- tryCatch({
    extract_instruments(outcomes = exp_id, p1 = instrument_p_threshold,
                        clump = TRUE, r2 = clump_r2, kb = clump_kb, pop = pop_ref)
  }, error = function(e) NULL)

  if (is.null(exp_iv) || nrow(exp_iv) < min_instruments) {
    message("  Instruments insufficient — skip.")
    next
  }

  # 循环结局
  for (j in 1:nrow(outcomes)) {
    out_id   <- outcomes$id[j]
    out_name <- outcomes$trait[j]
    out_sys  <- outcomes$System[j]

    message(sprintf("  [Outcome] %s (%s) | System: %s", out_name, out_id, out_sys))

    out_dat <- tryCatch({
      extract_outcome_data(snps = exp_iv$SNP, outcomes = out_id, proxies = TRUE)
    }, error = function(e) NULL)

    if (is.null(out_dat) || nrow(out_dat) == 0) next

    dat <- tryCatch({
      harmonise_data(exposure_dat = exp_iv, outcome_dat = out_dat)
    }, error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0) next

    mr_res <- tryCatch({
      mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
    }, error = function(e) NULL)
    if (is.null(mr_res) || nrow(mr_res) == 0) next

    mr_res$exposure_id <- exp_id
    mr_res$exposure    <- exp_name
    mr_res$outcome_id  <- out_id
    mr_res$outcome     <- out_name
    mr_res$System      <- out_sys

    results_all[[length(results_all)+1]] <- mr_res
  }
}

# -----------------------------
# 6. 保存结果 & 筛选显著
# -----------------------------
if (length(results_all) == 0) {
  stop("MR loop produced no results.")
}

res <- bind_rows(results_all)

# FDR 校正
res <- res %>%
  group_by(method) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ungroup()

# 保存所有结果
fwrite(res, "mr_results_all.tsv", sep="\t")

# 筛选显著结果
sig <- res %>% filter(FDR < 0.05)
fwrite(sig, "mr_results_significant.tsv", sep="\t")
message(sprintf("\nSignificant pairs (FDR<0.05): %d", nrow(sig)))

# -----------------------------
# 7. 生成下载清单
# -----------------------------
if (nrow(sig) > 0) {
  dl <- sig %>%
    select(exposure_id, exposure, outcome_id, outcome, System) %>%
    distinct()
  fwrite(dl, "sumstats_download_list.tsv", sep="\t")
  message("Saved: sumstats_download_list.tsv")
}
