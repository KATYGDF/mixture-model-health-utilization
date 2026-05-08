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
#   V   — Curvas de BIC e probabilidades posteriores (exames, univariado)
#   VI  — flexmix NB univariado: as 5 variáveis de V3 em profundidade
#         Funções exploradas: stepFlexmix · getModel · parameters · prior
#                             posterior · clusters · ICL · refit · summary
#
# Pacotes requeridos:
#   MASS      — glm.nb (Binomial Negativa)
#   mixtools  — normalmixEM (mistura Gaussiana, via EM)
#   flexmix   — flexmix + stepFlexmix + FLXMRglm + FLXMRnegbin
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
N_REP_FLEX  <- 5L   # reinícios por g no stepFlexmix


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

stats_mv <- data.frame(
  variavel  = VARS_V3,
  media     = sapply(VARS_V3, function(v) round(mean(df[[v]]), 2)),
  variancia = sapply(VARS_V3, function(v) round(var(df[[v]]),  2)),
  zeros_pct = sapply(VARS_V3, function(v) round(100 * mean(df[[v]] == 0), 1))
)
stats_mv$razao_v_m <- round(stats_mv$variancia / stats_mv$media, 1)

cat("\nRazão Var/Média por variável (Poisson assume razão = 1):\n")
print(stats_mv, row.names = FALSE)

y_exames <- df$exames
fit_pois  <- glm(y_exames ~ 1, family = poisson)
fit_nb    <- MASS::glm.nb(y_exames ~ 1)

cat(sprintf("\nexames — AIC Poisson: %.1f | AIC NB: %.1f (ganho: %.1f)\n",
            AIC(fit_pois), AIC(fit_nb), AIC(fit_pois) - AIC(fit_nb)))
cat(sprintf("         theta NB: %.4f | sobredispersão 1/theta = %.4f\n",
            fit_nb$theta, 1 / fit_nb$theta))


# ══ PARTE II: Heterogeneidade latente ════════════════════════════════════════

cat("\n══ PARTE II: Heterogeneidade latente ═══════════════════════════\n")

resid_nb <- residuals(fit_nb, type = "pearson")
cat(sprintf("Resíduos de Pearson NB: média=%.3f · sd=%.3f\n",
            mean(resid_nb), sd(resid_nb)))

kurt_nb <- mean((resid_nb - mean(resid_nb))^4) / var(resid_nb)^2
cat(sprintf("Curtose dos resíduos: %.2f (Normal = 3; > 3 ⟹ estrutura residual)\n",
            kurt_nb))
cat("Uma NB global homogênea não captura subpopulações distintas.\n")


# ══ PARTE III: mixtools — mistura Gaussiana em log1p ═════════════════════════

cat("\n══ PARTE III: mixtools — mistura Gaussiana em log1p ════════════\n")

y_log <- log1p(y_exames)

bic_gauss <- numeric(G_MAX - 1)
for (g in 2:G_MAX) {
  set.seed(SEED)
  fit_tmp    <- suppressWarnings(
    mixtools::normalmixEM(y_log, k = g, maxit = 500, epsilon = 1e-6)
  )
  k_par          <- 3L * g - 1L
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


# ══ PARTE IV: flexmix — Poisson e NB ═════════════════════════════════════════

cat("\n══ PARTE IV: flexmix — mistura Poisson e NB ════════════════════\n")

# ── IV.a: Mistura Poisson ─────────────────────────────────────────────────────
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

# ── IV.b: Mistura NB ──────────────────────────────────────────────────────────
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


# ══ PARTE V: Probabilidades posteriores e resumo (exames, univariado) ════════

cat("\n══ PARTE V: Probabilidades posteriores ═════════════════════════\n")

post_nb  <- flexmix::posterior(mod_nb_best)
prob_max <- apply(post_nb, 1, max)

cat(sprintf("  Prob. máxima — média: %.3f · mediana: %.3f · pct>0,9: %.1f%%\n",
            mean(prob_max), median(prob_max), 100 * mean(prob_max > 0.9)))

cat("\n── Resumo evolução (exames, univariado) ─────────────────────────\n")
cat(sprintf("  GLM Poisson (g=1):    AIC = %.1f\n", AIC(fit_pois)))
cat(sprintf("  GLM NB      (g=1):    AIC = %.1f\n", AIC(fit_nb)))
cat(sprintf("  Mistura Poisson:      BIC = %.1f  (g=%d)\n",
            min(bic_pois_fx), g_pois_fx))
cat(sprintf("  Mistura NB (flexmix): BIC = %.1f  (g=%d)\n",
            min(bic_nb_validos), g_nb_fx))


# ══ PARTE VI: flexmix NB univariado — as 5 variáveis de V3 ══════════════════
# Explora as funcionalidades do flexmix em profundidade para cada variável:
#
#   stepFlexmix(formula, data, k, nrep, model, control)
#     Ajusta flexmix para cada valor de k com nrep reinícios aleatórios cada.
#     Mantém o melhor resultado (maior log-verossimilhança) por k.
#     Retorna objeto S4 class "stepFlexmix" com slot @models.
#
#   getModel(objeto_stepFlexmix, which)
#     which = inteiro  → modelo com k componentes específico
#     which = "BIC"    → modelo com menor BIC
#     which = "AIC"    → modelo com menor AIC
#     which = "ICL"    → modelo com menor ICL
#
#   parameters(modelo)   → matriz de parâmetros por componente
#   prior(modelo)        → vetor de pesos π̂_h (proporções de mistura)
#   posterior(modelo)    → matriz n × g de probabilidades posteriores τ̂_ih
#   clusters(modelo)     → vetor n de classificações MAP (argmax da posterior)
#   ICL(modelo)          → integrated completed likelihood (BIC + penalidade de entropia)
#   refit(modelo)        → re-estima cada componente para obter erros-padrão
#   summary(refit(...))  → coeficientes, SE, z-valor e p-valor por componente

cat("\n\n══ PARTE VI: flexmix NB univariado — as 5 variáveis ═══════════\n")

resultados_univ <- list()

for (var in VARS_V3) {

  y <- df[[var]]
  n_zero <- sum(y == 0)
  cat(sprintf("\n%s\n", strrep("─", 60)))
  cat(sprintf("  Variável: %s  (n_zeros=%d · %.1f%% · média=%.2f · Var/Média=%.1f)\n",
              var, n_zero, 100 * mean(y == 0), mean(y), var(y) / mean(y)))
  cat(strrep("─", 60), "\n")

  # ── stepFlexmix: g = 1..G_MAX, N_REP_FLEX reinícios cada ─────────────────
  # Vantagem sobre loop manual: guarda automaticamente o melhor resultado
  # (por log-lik) para cada k, permite getModel() por critério.
  set.seed(SEED)
  fx_steps <- suppressWarnings(
    flexmix::stepFlexmix(
      as.formula(paste0(var, " ~ 1")),
      data    = df,
      k       = 1:G_MAX,
      nrep    = N_REP_FLEX,
      model   = flexmix::FLXMRnegbin(),
      control = list(iter.max = 300, minprior = 0.02)
    )
  )

  # ── Tabela de critérios por g: logLik · BIC · AIC · ICL ──────────────────
  crit_tab <- do.call(rbind, lapply(seq_len(G_MAX), function(g) {
    mod <- tryCatch(flexmix::getModel(fx_steps, which = g), error = function(e) NULL)
    if (is.null(mod)) return(data.frame(g=g, k_livre=NA, logLik=NA,
                                        BIC=NA, AIC=NA, ICL=NA))
    ll  <- as.numeric(logLik(mod))
    df_ll <- attr(logLik(mod), "df")   # número de parâmetros livres
    data.frame(
      g      = g,
      k_livre = df_ll,
      logLik = round(ll,           1),
      BIC    = round(BIC(mod),     1),
      AIC    = round(AIC(mod),     1),
      ICL    = round(flexmix::ICL(mod), 1)
    )
  }))

  cat("\n  Critérios por g (BIC = –2ℓ + d·log n  ·  ICL = BIC + 2·Entropia(τ)):\n")
  print(crit_tab, row.names = FALSE)

  # ── Melhor g por BIC: getModel(fx_steps, "BIC") ──────────────────────────
  best      <- flexmix::getModel(fx_steps, which = "BIC")
  g_best    <- best@k
  cat(sprintf("\n  g* (BIC) = %d  ·  BIC = %.1f  ·  logLik = %.1f\n",
              g_best, BIC(best), as.numeric(logLik(best))))

  # ── prior(): pesos dos componentes π̂_h ───────────────────────────────────
  # Interpretação: proporção estimada de beneficiários em cada componente
  pesos <- flexmix::prior(best)
  cat(sprintf("  Pesos π̂_h:  %s\n",
              paste(sprintf("G%d=%.3f", seq_along(pesos), pesos), collapse = "  ")))

  # ── parameters(): parâmetros estimados por componente ────────────────────
  # FLXMRnegbin usa log-link → linha 1 = log(µ), linha 2 = θ (size)
  pars   <- flexmix::parameters(best)
  mu     <- exp(pars[1, ])          # µ em escala original
  theta  <- pars[2, ]               # θ (parâmetro de dispersão NB₂)
  var_teo <- mu + mu^2 / theta      # Var[Y] = µ + µ²/θ

  pars_tab <- data.frame(
    comp       = paste0("G", seq_len(g_best)),
    peso       = round(pesos, 4),
    mu         = round(mu,    3),
    theta      = round(theta, 3),
    var_teorica = round(var_teo, 2),
    var_mu      = round(var_teo / mu, 2)   # razão Var/µ (= 1 seria Poisson)
  )
  cat("\n  Parâmetros por componente (µ = escala original, θ = size NB₂):\n")
  print(pars_tab, row.names = FALSE)

  # ── posterior(): probabilidades posteriores τ̂_ih ─────────────────────────
  # Matriz n × g. posterior(model)[i, h] = P(Z_i = h | Y = y_i)
  post     <- flexmix::posterior(best)
  max_post <- apply(post, 1, max)

  # Entropia por observação: H_i = -Σ_h τ_ih · log(τ_ih)
  # H = 0 → certeza total; H = log(g) → incerteza máxima (uniforme)
  ent     <- -rowSums(post * log(post + 1e-15))
  ent_max <- log(g_best)    # entropia máxima possível para g componentes

  cat(sprintf(
    "\n  Certeza (max τ̂_ih): média=%.3f · mediana=%.3f · %%>0.9=%.1f%%\n",
    mean(max_post), median(max_post), 100 * mean(max_post > 0.9)
  ))
  cat(sprintf(
    "  Entropia média: %.4f  (máx teórico = log(%d) = %.4f  ·  relativa = %.1f%%)\n",
    mean(ent), g_best, ent_max, 100 * mean(ent) / ent_max
  ))

  # ── clusters(): classificação MAP (hard assignment) ───────────────────────
  # clusters(model) = apply(posterior(model), 1, which.max)
  cl <- flexmix::clusters(best)
  tab_cl <- table(real = df$grupo_real, estimado = cl)
  cat("\n  Tabela real × estimado (hard assignment via MAP):\n")
  print(tab_cl)

  # Tamanho de cada componente estimado
  tab_comp <- table(estimado = cl)
  cat(sprintf("  Tamanho dos componentes: %s\n",
              paste(sprintf("G%s=%d", names(tab_comp), as.integer(tab_comp)),
                    collapse = "  ")))

  # ── refit() + summary(): coeficientes com erros-padrão ───────────────────
  # refit() re-estima os parâmetros de cada componente separadamente,
  # calculando a Hessiana para obter SE via Fisher information.
  # summary(refit(model)) exibe coef + SE + z + p por componente.
  ref <- tryCatch(flexmix::refit(best), error = function(e) NULL)
  if (!is.null(ref)) {
    cat("\n  summary(refit) — coeficientes com erros-padrão por componente:\n")
    s <- summary(ref)
    # Extrai tabela de cada componente
    for (h in seq_len(g_best)) {
      comp_nm <- paste0("Comp.", h)
      if (!is.null(s@coef[[comp_nm]])) {
        cat(sprintf("    Componente %d (peso=%.3f):\n", h, pesos[h]))
        ct <- s@coef[[comp_nm]]
        # log(µ): exp para mostrar µ em escala original
        if ("coef.(Intercept)" %in% rownames(ct)) {
          ct_show <- ct
          ct_show["coef.(Intercept)", "Estimate"] <-
            exp(ct["coef.(Intercept)", "Estimate"])
          rownames(ct_show)[rownames(ct_show) == "coef.(Intercept)"] <- "mu (orig.)"
        } else {
          ct_show <- ct
        }
        print(ct_show)
      }
    }
  }

  # Armazenar resultado completo
  resultados_univ[[var]] <- list(
    fx_steps  = fx_steps,
    best      = best,
    g_best    = g_best,
    crit_tab  = crit_tab,
    pars_tab  = pars_tab,
    post      = post,
    max_post  = max_post,
    ent       = ent,
    clusters  = cl
  )
}


# ── Tabela resumo comparativa das 5 variáveis ─────────────────────────────────

cat(sprintf("\n\n%s\n", strrep("═", 60)))
cat("RESUMO COMPARATIVO — flexmix NB univariado — V3\n")
cat(strrep("═", 60), "\n\n")

resumo_univ <- do.call(rbind, lapply(VARS_V3, function(var) {
  r   <- resultados_univ[[var]]
  bst <- r$best
  data.frame(
    variavel      = var,
    g_star        = r$g_best,
    BIC_star      = round(r$crit_tab$BIC[r$g_best], 1),
    pi_min        = round(min(flexmix::prior(bst)), 3),
    pi_max        = round(max(flexmix::prior(bst)), 3),
    mu_min        = round(min(exp(flexmix::parameters(bst)[1, ])), 2),
    mu_max        = round(max(exp(flexmix::parameters(bst)[1, ])), 2),
    post_mediana  = round(median(r$max_post), 3),
    pct_cert_90   = round(100 * mean(r$max_post > 0.9), 1),
    ent_rel_pct   = round(100 * mean(r$ent) / log(max(r$g_best, 2)), 1)
  )
}))

print(resumo_univ, row.names = FALSE)
cat("\nLegenda:\n")
cat("  g_star       = g selecionado pelo BIC\n")
cat("  pi_min/max   = menor e maior peso de componente\n")
cat("  mu_min/max   = menor e maior média estimada (escala original)\n")
cat("  post_mediana = mediana da max(τ̂_ih) — certeza mediana de classificação\n")
cat("  pct_cert_90  = % de beneficiários com max(τ̂_ih) > 0,90\n")
cat("  ent_rel_pct  = entropia média como % da entropia máxima (menor = melhor separação)\n")

cat("\nNota: avaliação UNIVARIADA — cada variável independentemente.\n")
cat("O EM customizado (comparar_variaveis.R) opera sobre Y multivariado (V3: 5 vars).\n")


# ── 7. Salvar resultados ──────────────────────────────────────────────────────

if (!dir.exists(DIR_RDS)) dir.create(DIR_RDS, recursive = TRUE)

saveRDS(
  list(
    stats_mv       = stats_mv,
    fit_pois       = fit_pois,
    fit_nb         = fit_nb,
    fit_gauss      = fit_gauss,
    bic_gauss      = bic_gauss,
    bic_pois_fx    = bic_pois_fx,
    bic_nb_fx      = bic_nb_fx,
    mod_nb_best    = mod_nb_best,
    post_nb        = post_nb,
    resultados_univ = resultados_univ,
    resumo_univ    = resumo_univ
  ),
  file.path(DIR_RDS, "construcao_progressiva.rds")
)

cat(sprintf("\nSalvo em '%s/construcao_progressiva.rds'\n", DIR_RDS))
cat("Próximo passo: analises.html — EM customizado multivariado (V1–V3)\n")
