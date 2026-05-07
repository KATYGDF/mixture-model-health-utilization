# =============================================================================
# construcao_progressiva.R
# Evolução metodológica: da Poisson à Mistura NB customizada
# Variáveis: V3 = consultas + ps + exames + internacoes + terapias
#
# Partes:
#   I   — Sobredispersão: Poisson vs Binomial Negativa (GLM univariado)
#   II  — Heterogeneidade latente: NB simples não captura subpopulações
#   III — Primeiros modelos de mistura: mixtools (Gaussiana em log1p)
#   IV  — Mixturas com famílias de contagem: flexmix (Poisson e NB)
#   V   — Curvas de BIC e probabilidades posteriores
#
# Pacotes requeridos:
#   MASS      — glm.nb (Binomial Negativa)
#   mixtools  — normalmixEM (mistura Gaussiana, via EM)
#   flexmix   — flexmix + FLXMRglm + FLXMRnegbin (mixturas GLM)
#
# USO: Source (Ctrl+Shift+S) no RStudio
#      Resultados salvos em rds/construcao_progressiva.rds
# =============================================================================


# ── 1. Parâmetros ─────────────────────────────────────────────────────────────

CAMINHO_CSV <- "../../dados/dados_simulados.csv"
SEED        <- 42L
G_MAX       <- 6L
VARS_V3     <- c("consultas", "ps", "exames", "internacoes", "terapias")
DIR_RDS     <- "rds"


# ── 2. Pacotes ────────────────────────────────────────────────────────────────

pkgs  <- c("MASS", "mixtools", "flexmix")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Instalando: ", paste(novos, collapse = ", "))
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(MASS); library(mixtools); library(flexmix)
})


# ── 3. Dados ──────────────────────────────────────────────────────────────────

df_raw <- read.csv(CAMINHO_CSV, stringsAsFactors = FALSE, encoding = "UTF-8")

candidatos <- c("grupo_latente", "grupo", "grupo_real", "grupo_idx")
col_grupo  <- candidatos[candidatos %in% names(df_raw)][1]
if (is.na(col_grupo)) stop("Coluna de grupo não encontrada.")

df <- df_raw
names(df)[names(df) == col_grupo] <- "grupo_real"
df$grupo_real <- factor(df$grupo_real)

cat(sprintf("\nDados: %d obs · %d grupos reais\n", nrow(df), nlevels(df$grupo_real)))
cat("Variáveis V3:", paste(VARS_V3, collapse = ", "), "\n\n")


# ══ PARTE I: Sobredispersão ═══════════════════════════════════════════════════

cat("══ PARTE I: Sobredispersão ══════════════════════════════════════\n")

# Estatísticas descritivas por variável
stats_mv <- data.frame(
  variavel  = VARS_V3,
  media     = sapply(VARS_V3, function(v) round(mean(df[[v]]), 2)),
  variancia = sapply(VARS_V3, function(v) round(var(df[[v]]),  2)),
  zeros_pct = sapply(VARS_V3, function(v) round(100 * mean(df[[v]] == 0), 1))
)
stats_mv$razao_v_m <- round(stats_mv$variancia / stats_mv$media, 1)

cat("\nRazão Var/Média por variável (Poisson assume razão = 1):\n")
print(stats_mv, row.names = FALSE)

# GLM Poisson vs Binomial Negativa — variável exames (mais overdispersa)
y_exames <- df$exames
fit_pois  <- glm(y_exames ~ 1, family = poisson)
fit_nb    <- MASS::glm.nb(y_exames ~ 1)

cat(sprintf("\nexames — AIC Poisson: %.1f | AIC NB: %.1f (ganho: %.1f)\n",
            AIC(fit_pois), AIC(fit_nb), AIC(fit_pois) - AIC(fit_nb)))
cat(sprintf("         theta NB: %.4f | sobredispersão 1/theta = %.4f\n",
            fit_nb$theta, 1 / fit_nb$theta))


# ══ PARTE II: Heterogeneidade latente ════════════════════════════════════════

cat("\n══ PARTE II: Heterogeneidade latente ═══════════════════════════\n")

# Resíduos de Pearson do modelo NB univariado
resid_nb <- residuals(fit_nb, type = "pearson")
cat(sprintf("Resíduos de Pearson NB: média=%.3f · sd=%.3f\n",
            mean(resid_nb), sd(resid_nb)))

# Curtose empírica dos resíduos (> 3 indica cauda pesada ou bimodalidade)
kurt_nb <- mean((resid_nb - mean(resid_nb))^4) / var(resid_nb)^2
cat(sprintf("Curtose dos resíduos: %.2f (Normal = 3; > 3 ⟹ estrutura residual)\n",
            kurt_nb))
cat("Uma NB global homogênea não captura subpopulações distintas.\n")


# ══ PARTE III: mixtools — mistura Gaussiana em log1p ═════════════════════════

cat("\n══ PARTE III: mixtools — mistura Gaussiana em log1p ════════════\n")

y_log <- log1p(y_exames)

# BIC manual para normalmixEM (g = 2..G_MAX)
# Nota: mixtools não calcula BIC diretamente — calculado a partir do log-lik
bic_gauss <- numeric(G_MAX - 1)
for (g in 2:G_MAX) {
  set.seed(SEED)
  fit_tmp    <- suppressWarnings(
    mixtools::normalmixEM(y_log, k = g, maxit = 500, epsilon = 1e-6)
  )
  k_par          <- 3L * g - 1L   # g médias + g dp + (g-1) pesos
  bic_gauss[g-1] <- -2 * fit_tmp$loglik + k_par * log(length(y_log))
  cat(sprintf("  [g=%d: BIC=%.1f]", g, bic_gauss[g-1]))
}
cat("\n")

g_gauss <- which.min(bic_gauss) + 1L
cat(sprintf("  g selecionado: %d\n", g_gauss))

set.seed(SEED)
fit_gauss <- suppressWarnings(
  mixtools::normalmixEM(y_log, k = g_gauss, maxit = 500, epsilon = 1e-6)
)
cat(sprintf("  Pesos:  %s\n", paste(round(fit_gauss$lambda, 3), collapse = " | ")))
cat(sprintf("  Médias: %s  (escala log1p)\n",
            paste(round(fit_gauss$mu, 3), collapse = " | ")))
cat(sprintf("  DPs:    %s\n", paste(round(fit_gauss$sigma, 3), collapse = " | ")))
cat("\n  Limitação: suporte contínuo ℝ ≠ suporte discreto ℕ₀.\n")
cat("  P(Y < 0) > 0 em escala original — incompatível com contagens.\n")


# ══ PARTE IV: flexmix — Poisson e NB ═════════════════════════════════════════

cat("\n══ PARTE IV: flexmix — mistura Poisson e NB ════════════════════\n")

# ── IV.a: Mistura Poisson (univariado, exames) ────────────────────────────────
cat("  4a. Mistura Poisson · exames · g = 2..6\n")

bic_pois_fx <- sapply(2:G_MAX, function(g) {
  set.seed(SEED)
  mod <- flexmix::flexmix(
    y_exames ~ 1, k = g,
    model   = flexmix::FLXMRglm(family = "poisson"),
    control = list(iter.max = 300, minprior = 0.02)
  )
  cat(sprintf("     [g=%d: BIC=%.1f]", g, BIC(mod)))
  BIC(mod)
})
cat("\n")
g_pois_fx <- which.min(bic_pois_fx) + 1L
cat(sprintf("     g selecionado: %d\n", g_pois_fx))
cat("     Limitação: equidispersão forçada DENTRO de cada componente.\n")

# ── IV.b: Mistura NB (univariado, exames) ─────────────────────────────────────
cat("\n  4b. Mistura NB · exames · g = 2..6\n")

bic_nb_fx  <- numeric(G_MAX - 1)
fits_nb_fx <- vector("list", G_MAX - 1)

for (g in 2:G_MAX) {
  set.seed(SEED)
  mod <- tryCatch(
    flexmix::flexmix(
      y_exames ~ 1, k = g,
      model   = flexmix::FLXMRnegbin(),
      control = list(iter.max = 300, minprior = 0.02)
    ),
    error = function(e) NULL
  )
  if (!is.null(mod)) {
    bic_nb_fx[g-1]    <- BIC(mod)
    fits_nb_fx[[g-1]] <- mod
    cat(sprintf("     [g=%d: BIC=%.1f]", g, BIC(mod)))
  }
}
cat("\n")

bic_nb_validos <- bic_nb_fx[bic_nb_fx > 0]
g_nb_fx        <- which.min(bic_nb_validos) + 1L
mod_nb_best    <- fits_nb_fx[[g_nb_fx - 1L]]

cat(sprintf("     g selecionado: %d\n", g_nb_fx))
cat("     Família NB adequada — mas apenas para exames isolado.\n")


# ══ PARTE V: Probabilidades posteriores e resumo final ═══════════════════════

cat("\n══ PARTE V: Probabilidades posteriores ═════════════════════════\n")

post_nb  <- flexmix::posterior(mod_nb_best)
prob_max <- apply(post_nb, 1, max)

cat(sprintf("  Prob. máxima — média: %.3f · mediana: %.3f · pct>0,9: %.1f%%\n",
            mean(prob_max), median(prob_max), 100 * mean(prob_max > 0.9)))

cat("\n── Resumo da evolução (AIC/BIC): exames ─────────────────────────\n")
cat(sprintf("  GLM Poisson (1 comp):   AIC = %.1f\n", AIC(fit_pois)))
cat(sprintf("  GLM NB      (1 comp):   AIC = %.1f\n", AIC(fit_nb)))
cat(sprintf("  Mistura Poisson:        BIC = %.1f  (g=%d)\n",
            min(bic_pois_fx), g_pois_fx))
cat(sprintf("  Mistura NB (flexmix):   BIC = %.1f  (g=%d)\n",
            min(bic_nb_validos), g_nb_fx))
cat("\n")
cat("Nota: avaliação acima é UNIVARIADA (apenas exames).\n")
cat("O EM customizado (analises.R) opera sobre Y multivariado (5 vars).\n")
cat("Mistura NB multivariada com EM customizado: ARI=0,729 · g=4 em V3.\n")


# ── 6. Salvar resultados ──────────────────────────────────────────────────────

if (!dir.exists(DIR_RDS)) dir.create(DIR_RDS, recursive = TRUE)

saveRDS(
  list(
    stats_mv    = stats_mv,
    fit_pois    = fit_pois,
    fit_nb      = fit_nb,
    fit_gauss   = fit_gauss,
    bic_gauss   = bic_gauss,
    bic_pois_fx = bic_pois_fx,
    bic_nb_fx   = bic_nb_fx,
    mod_nb_best = mod_nb_best,
    post_nb     = post_nb
  ),
  file.path(DIR_RDS, "construcao_progressiva.rds")
)

cat(sprintf("\nSalvo em '%s/construcao_progressiva.rds'\n", DIR_RDS))
cat("Próximo passo: analises.html — EM customizado multivariado (V1–V5)\n")
