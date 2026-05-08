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

CAMINHO_CSV  <- "../../dados/dados_simulados.csv"
SEED         <- 42L
G_MAX        <- 6L
VARS_V3      <- c("consultas", "ps", "exames", "internacoes", "terapias")
DIR_RDS      <- "rds"
N_REP_FLEX   <- 5L    # reinícios por g no stepFlexmix
DIR_FIGS_HTML <- "../../docs/assets/img/construcao"   # figuras para a página web


# ── 2. Pacotes ────────────────────────────────────────────────────────────────

pkgs  <- c("MASS", "mixtools", "flexmix", "ggplot2", "patchwork", "scales")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Instalando: ", paste(novos, collapse = ", "))
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(MASS); library(mixtools); library(flexmix)
  library(ggplot2); library(patchwork); library(scales)
})


# ── 2b. FLXMRnegbin (compatibilidade) ────────────────────────────────────────
# FLXMRnegbin foi adicionado ao flexmix a partir da versão 2.3-14.
# Se a versão instalada não a contém, definimos uma versão local com API idêntica:
#   parameters(model) → matriz 2×g: linha 1 = log(µ_h), linha 2 = θ_h

if (!tryCatch({ flexmix::FLXMRnegbin; TRUE }, error = function(e) FALSE)) {
  FLXMRnegbin <- function(formula = . ~ .) {
    z <- methods::new("FLXMR", formula = formula, name = "FLXMRnegbin",
                      weighted = TRUE)

    z@defineComponent <- function(para) {
      mu_log <- para$mu_log
      theta  <- para$theta
      methods::new("FLXcomponent",
        # lista de comprimento 1 → drop extrai vetor nomeado → cbind → matriz 2×g
        parameters = list(coef = c("(Intercept)" = mu_log, "theta" = theta)),
        predict    = function(x, ...) matrix(exp(mu_log), nrow(x), 1L),
        logLik     = function(x, y)
          dnbinom(as.integer(y), mu = exp(mu_log), size = theta, log = TRUE),
        df = 2L
      )
    }

    z@fit <- function(x, y, w, ...) {
      y_int <- as.integer(y)
      w_n   <- w / sum(w)                          # pesos normalizados
      mu0   <- pmax(sum(w_n * y_int), 1e-4)        # MLE fechado para µ
      opt   <- tryCatch(
        optimize(
          function(lt) -sum(w_n * dnbinom(y_int, mu = mu0,
                                          size = exp(lt), log = TRUE)),
          interval = c(-6, 8)
        ),
        error = function(e) NULL
      )
      theta <- if (!is.null(opt) && is.finite(opt$minimum)) exp(opt$minimum) else 1.0
      z@defineComponent(list(mu_log = log(mu0), theta = theta))
    }

    z
  }
  message("FLXMRnegbin definido localmente (flexmix < 2.3-14).")
}


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
      model   = FLXMRnegbin(),
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

# Usa NA em vez de filtrar: preserva o mapeamento índice→g (1 ↔ g=2, 2 ↔ g=3, …)
bic_nb_na <- replace(bic_nb_fx, bic_nb_fx == 0, NA_real_)

if (all(is.na(bic_nb_na))) {
  warning("Nenhum ajuste flexmix NB convergiu para exames (FLXMRnegbin falhou em todos os g).")
  mod_nb_best <- NULL
  g_nb_fx     <- NA_integer_
} else {
  idx_best_nb <- which.min(bic_nb_na)    # posição em bic_nb_fx (1 = g=2, …)
  g_nb_fx     <- idx_best_nb + 1L        # g real
  mod_nb_best <- fits_nb_fx[[idx_best_nb]]
  cat(sprintf("     g selecionado: %d\n", g_nb_fx))
}

bic_nb_validos <- bic_nb_na[!is.na(bic_nb_na)]   # mantém para prints abaixo


# ══ PARTE V: Probabilidades posteriores e resumo (exames, univariado) ════════

cat("\n══ PARTE V: Probabilidades posteriores ═════════════════════════\n")

if (!is.null(mod_nb_best)) {
  post_nb  <- flexmix::posterior(mod_nb_best)
  prob_max <- apply(post_nb, 1, max)
  cat(sprintf("  Prob. máxima — média: %.3f · mediana: %.3f · pct>0,9: %.1f%%\n",
              mean(prob_max), median(prob_max), 100 * mean(prob_max > 0.9)))
} else {
  post_nb  <- NULL
  prob_max <- NULL
  cat("  (ajuste NB não disponível — Parte V ignorada)\n")
}

cat("\n── Resumo evolução (exames, univariado) ─────────────────────────\n")
cat(sprintf("  GLM Poisson (g=1):    AIC = %.1f\n", AIC(fit_pois)))
cat(sprintf("  GLM NB      (g=1):    AIC = %.1f\n", AIC(fit_nb)))
cat(sprintf("  Mistura Poisson:      BIC = %.1f  (g=%d)\n",
            min(bic_pois_fx), g_pois_fx))
if (length(bic_nb_validos) > 0 && !is.na(g_nb_fx)) {
  cat(sprintf("  Mistura NB (flexmix): BIC = %.1f  (g=%d)\n",
              min(bic_nb_validos), g_nb_fx))
} else {
  cat("  Mistura NB (flexmix): sem convergência\n")
}


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
      model   = FLXMRnegbin(),
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


# ══ PARTE VII: Visualizações ═════════════════════════════════════════════════

cat("\n══ PARTE VII: Gerando figuras ══════════════════════════════════\n")

for (d in c(DIR_RDS, DIR_FIGS_HTML))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Paleta e tema consistentes com comparar_variaveis.R
COR_COMP <- setNames(
  c("#4a7c59","#2980b9","#e6a817","#c0392b","#8e44ad","#16a085"),
  as.character(1:6)
)
TEMA <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 9),
    strip.background  = element_rect(fill = "#f0f0f0", colour = NA),
    legend.position   = "bottom",
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(colour = "grey40", size = 9),
    plot.caption      = element_text(colour = "grey55", size = 8, hjust = 0)
  )

salvar <- function(p, nome, w = 9, h = 5) {
  caminho <- file.path(DIR_FIGS_HTML, nome)
  ggsave(caminho, p, width = w, height = h, dpi = 150, bg = "white")
  cat(sprintf("  %s\n", nome))
}


# ── FIG 01 · Sobredispersão: razão Var/Média por variável ────────────────────

df_sd <- stats_mv
df_sd$cor <- ifelse(df_sd$razao_v_m > 10, "alta",
             ifelse(df_sd$razao_v_m > 5,  "moderada", "baixa"))
df_sd$variavel <- factor(df_sd$variavel, levels = rev(VARS_V3))

p01 <- ggplot(df_sd, aes(razao_v_m, variavel, fill = cor)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f×", razao_v_m)),
            hjust = -0.15, size = 3.5) +
  geom_vline(xintercept = 1, colour = "grey30", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 1.3, y = 0.55, label = "Poisson\n(razão = 1)",
           hjust = 0, size = 3, colour = "grey30") +
  scale_fill_manual(
    values = c(alta = "#c0392b", moderada = "#e6a817", baixa = "#4a7c59"),
    name = "Sobredispersão", guide = guide_legend(reverse = TRUE)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = "Razão Variância / Média",
    y = NULL,
    title = "Sobredispersão por variável — V3",
    subtitle = "Razão > 1 viola a hipótese de equidispersão da Poisson",
    caption = sprintf("n = %d beneficiários  ·  grupos simulados = %d",
                      nrow(df), nlevels(df$grupo_real))
  ) +
  TEMA + theme(legend.position = "right")

salvar(p01, "fig01_sobredispersao.png", w = 8, h = 4.5)


# ── FIG 02 · Poisson vs NB: PMF empírica × ajustada (exames) ─────────────────

x_max02 <- as.integer(quantile(y_exames, 0.97))
x_seq02 <- 0:x_max02
mu_pois  <- mean(y_exames)
mu_nb    <- exp(coef(fit_nb)[1])
th_nb    <- fit_nb$theta

emp02 <- data.frame(
  y    = as.integer(names(table(y_exames[y_exames <= x_max02]))),
  prop = as.numeric(table(y_exames[y_exames <= x_max02])) / length(y_exames)
)
mod02 <- rbind(
  data.frame(y = x_seq02, prob = dpois(x_seq02, lambda = mu_pois),
             modelo = sprintf("Poisson (λ=%.1f)", mu_pois)),
  data.frame(y = x_seq02, prob = dnbinom(x_seq02, mu = mu_nb, size = th_nb),
             modelo = sprintf("NB (μ=%.1f, θ=%.2f)", mu_nb, th_nb))
)

p02 <- ggplot() +
  geom_col(data = emp02, aes(y, prop),
           fill = "grey82", colour = "grey70", width = 0.9) +
  geom_line(data = mod02, aes(y, prob, colour = modelo), linewidth = 1.1) +
  geom_point(data = mod02, aes(y, prob, colour = modelo), size = 1.5) +
  scale_colour_manual(
    values = c("#c0392b", "#2980b9"), name = NULL,
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  scale_x_continuous(breaks = seq(0, x_max02, by = 5)) +
  labs(
    x = "exames (contagem)",
    y = "Proporção / Probabilidade",
    title = "PMF empírica de exames vs. ajuste Poisson e NB",
    subtitle = "Barras = frequência observada  ·  Linhas = PMF teórica",
    caption = sprintf(
      "AIC Poisson = %.0f  ·  AIC NB = %.0f  ·  ganho = %.0f  ·  θ̂ = %.3f",
      AIC(fit_pois), AIC(fit_nb), AIC(fit_pois) - AIC(fit_nb), th_nb
    )
  ) +
  TEMA + theme(legend.position = c(0.72, 0.78),
               legend.background = element_rect(fill = "white", colour = NA))

salvar(p02, "fig02_pois_vs_nb.png", w = 9, h = 5)


# ── FIG 03 · Resíduos de Pearson do GLM NB global ───────────────────────────

df_res <- data.frame(resid = resid_nb)
p03 <- ggplot(df_res, aes(resid)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = "#2980b9", colour = "white", alpha = 0.7) +
  geom_density(colour = "#1a5276", linewidth = 0.9) +
  stat_function(fun = dnorm, colour = "#c0392b", linewidth = 0.9, linetype = "dashed",
                args = list(mean = mean(resid_nb), sd = sd(resid_nb))) +
  annotate("text", x = max(resid_nb) * 0.6, y = Inf,
           label = sprintf("Curtose = %.2f\n(Normal = 3.00)", kurt_nb),
           vjust = 1.4, hjust = 0.5, size = 3.3, colour = "#c0392b") +
  labs(
    x = "Resíduo de Pearson",
    y = "Densidade",
    title = "Resíduos de Pearson — GLM Binomial Negativa global (exames)",
    subtitle = "Azul = densidade empírica  ·  Vermelho tracejado = Normal com mesma média e DP",
    caption = "Curtose >> 3 e bimodalidade nos resíduos indicam mistura de subpopulações"
  ) +
  TEMA

salvar(p03, "fig03_residuos_nb.png", w = 8.5, h = 5)


# ── FIG 04 · Mistura Gaussiana em log1p(exames) — mixtools ───────────────────

x_grid04 <- seq(min(y_log) - 0.2, max(y_log) + 0.2, length.out = 400)

comp04 <- do.call(rbind, lapply(seq_len(g_gauss), function(h) {
  data.frame(
    x    = x_grid04,
    dens = fit_gauss$lambda[h] * dnorm(x_grid04,
                                       mean = fit_gauss$mu[h],
                                       sd   = fit_gauss$sigma[h]),
    comp = paste0("G", h)
  )
}))
mix04 <- data.frame(
  x    = x_grid04,
  dens = rowSums(sapply(seq_len(g_gauss), function(h)
    fit_gauss$lambda[h] * dnorm(x_grid04, fit_gauss$mu[h], fit_gauss$sigma[h])))
)

p04 <- ggplot() +
  geom_histogram(data = data.frame(x = y_log),
                 aes(x, y = after_stat(density)),
                 bins = 40, fill = "grey85", colour = "white") +
  geom_line(data = comp04, aes(x, dens, colour = comp), linewidth = 1, linetype = "dashed") +
  geom_line(data = mix04,  aes(x, dens), colour = "black", linewidth = 1.2) +
  scale_colour_manual(values = unname(COR_COMP[1:g_gauss]), name = "Componente") +
  labs(
    x = "log(1 + exames)",
    y = "Densidade",
    title = sprintf("Mistura Gaussiana em log₁₊(exames) — %d componentes (mixtools)", g_gauss),
    subtitle = "Linhas coloridas = componentes ponderados  ·  Linha preta = densidade da mistura",
    caption = "Limitação: suporte ℝ da Normal inclui valores negativos impossíveis para contagens"
  ) +
  TEMA

salvar(p04, "fig04_gaussiana_mix.png", w = 8.5, h = 5)


# ── FIG 05 · Curvas BIC: Gaussiana / Poisson / NB — exames, g = 2..6 ────────

bic05 <- data.frame(
  g      = rep(2:G_MAX, 3),
  BIC    = c(bic_gauss, bic_pois_fx, bic_nb_na),   # bic_nb_na tem mesma
  modelo = rep(c("Gaussiana (log1p)", "Poisson", "Binomial Negativa"), each = G_MAX - 1)
  # comprimento que bic_gauss (G_MAX - 1), NAs excluídos pelo ggplot
)
bic05$modelo <- factor(bic05$modelo,
                        levels = c("Gaussiana (log1p)", "Poisson", "Binomial Negativa"))

# marcar mínimo por modelo
bic05_min <- do.call(rbind, lapply(levels(bic05$modelo), function(m) {
  sub <- bic05[bic05$modelo == m, ]
  sub[which.min(sub$BIC), ]
}))

p05 <- ggplot(bic05, aes(g, BIC, colour = modelo, group = modelo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_point(data = bic05_min, shape = 21, size = 5.5,
             stroke = 1.8, fill = NA, colour = "grey20") +
  scale_colour_manual(
    values = c("Gaussiana (log1p)" = "#8e44ad",
               "Poisson"           = "#e6a817",
               "Binomial Negativa" = "#2980b9"),
    name = NULL
  ) +
  scale_x_continuous(breaks = 2:G_MAX) +
  labs(
    x = "Número de componentes (g)",
    y = "BIC",
    title = "Curvas BIC por família de mistura — exames (univariado)",
    subtitle = "Círculo destacado = g* selecionado (mínimo BIC por família)",
    caption = "Nota: BIC NB < BIC Poisson < BIC Gaussiana — confirma que NB é a família correta para dados de contagem sobredispersos"
  ) +
  TEMA

salvar(p05, "fig05_bic_familias.png", w = 9, h = 5)


# ── FIG 06 · Histograma + componentes NB ajustados — as 5 variáveis ──────────
# Um painel por variável; barras = PMF empírica; linhas = componentes + mistura

plots_comp <- lapply(VARS_V3, function(var) {
  y_v    <- df[[var]]
  res_v  <- resultados_univ[[var]]
  pars_v <- flexmix::parameters(res_v$best)
  mu_v   <- exp(pars_v[1, ])
  th_v   <- pars_v[2, ]
  pi_v   <- flexmix::prior(res_v$best)
  g_v    <- res_v$g_best

  x_max_v <- as.integer(quantile(y_v, 0.97))
  x_seq_v <- 0:x_max_v

  # PMF empírica (truncada no p97)
  emp_v <- data.frame(
    y    = x_seq_v,
    prop = sapply(x_seq_v, function(yk) mean(y_v == yk))
  )

  # PMF da mistura
  mix_v <- data.frame(
    y   = x_seq_v,
    pmf = sapply(x_seq_v, function(yk)
      sum(pi_v * dnbinom(yk, mu = mu_v, size = th_v)))
  )

  # PMF por componente (contribuição ponderada)
  comp_v <- do.call(rbind, lapply(seq_len(g_v), function(h) {
    data.frame(
      y    = x_seq_v,
      pmf  = pi_v[h] * dnbinom(x_seq_v, mu = mu_v[h], size = th_v[h]),
      comp = paste0("G", h)
    )
  }))

  pesos_lbl <- paste(sprintf("G%d: μ=%.1f π=%.2f", seq_len(g_v), mu_v, pi_v),
                     collapse = "\n")

  ggplot() +
    geom_col(data = emp_v, aes(y, prop),
             fill = "grey83", colour = "grey70", width = 0.85) +
    geom_line(data = comp_v, aes(y, pmf, colour = comp),
              linewidth = 0.85, linetype = "dashed") +
    geom_line(data = mix_v,  aes(y, pmf),
              colour = "black", linewidth = 1.1) +
    scale_colour_manual(values = unname(COR_COMP[1:g_v]), name = NULL,
                        guide = guide_legend(nrow = 1)) +
    annotate("text", x = x_max_v * 0.55, y = Inf,
             label = pesos_lbl, vjust = 1.3, hjust = 0,
             size = 2.4, colour = "grey20", lineheight = 1.3) +
    labs(x = var, y = "Probabilidade",
         title = sprintf("%s  (g*=%d)", var, g_v),
         subtitle = sprintf("BIC* = %.0f  ·  mediana max(τ̂) = %.3f",
                            min(res_v$crit_tab$BIC, na.rm = TRUE),
                            median(res_v$max_post))) +
    TEMA +
    theme(legend.position  = "none",
          plot.title       = element_text(size = 10),
          plot.subtitle    = element_text(size = 8),
          axis.title.x     = element_text(face = "bold"))
})

fig06 <- wrap_plots(plots_comp, ncol = 3) +
  plot_annotation(
    title    = "Componentes NB ajustados por variável — flexmix NB univariado (V3)",
    subtitle = "Barras = PMF empírica  ·  Linhas coloridas = componentes ponderados  ·  Linha preta = mistura",
    caption  = "g* selecionado pelo BIC  ·  x truncado no percentil 97",
    theme    = TEMA
  )

salvar(fig06, "fig06_componentes_nb_v3.png", w = 15, h = 9)


# ── FIG 07 · Curvas BIC e ICL — as 5 variáveis, g = 1..6 ────────────────────

crit_long <- do.call(rbind, lapply(VARS_V3, function(var) {
  ct <- resultados_univ[[var]]$crit_tab
  rbind(
    data.frame(variavel = var, g = ct$g, valor = ct$BIC, criterio = "BIC"),
    data.frame(variavel = var, g = ct$g, valor = ct$ICL, criterio = "ICL")
  )
}))
crit_long$variavel <- factor(crit_long$variavel, levels = VARS_V3)
crit_long$criterio <- factor(crit_long$criterio, levels = c("BIC", "ICL"))

# marcar g* por variável e critério
g_star_df <- do.call(rbind, lapply(VARS_V3, function(var) {
  ct <- resultados_univ[[var]]$crit_tab
  rbind(
    data.frame(variavel = var, g = ct$g[which.min(ct$BIC)],
               valor = min(ct$BIC, na.rm=TRUE), criterio = "BIC"),
    data.frame(variavel = var, g = ct$g[which.min(ct$ICL)],
               valor = min(ct$ICL, na.rm=TRUE), criterio = "ICL")
  )
}))
g_star_df$variavel <- factor(g_star_df$variavel, levels = VARS_V3)
g_star_df$criterio <- factor(g_star_df$criterio, levels = c("BIC", "ICL"))

p07 <- ggplot(crit_long[!is.na(crit_long$valor), ],
              aes(g, valor, colour = variavel, group = variavel)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_point(data = g_star_df,
             shape = 21, size = 5, stroke = 1.6, fill = NA, colour = "grey20") +
  facet_wrap(~ criterio, scales = "free_y") +
  scale_x_continuous(breaks = 1:G_MAX) +
  scale_colour_manual(values = unname(COR_COMP[1:5]), name = "Variável") +
  labs(
    x = "Número de componentes (g)",
    y = "Valor do critério",
    title = "Seleção de g por variável — BIC e ICL (flexmix NB univariado)",
    subtitle = "Círculo destacado = g* selecionado  ·  ICL penaliza sobreposição de componentes",
    caption = "ICL = BIC + 2·H(τ)  onde H(τ) = entropia das probabilidades posteriores"
  ) +
  TEMA

salvar(p07, "fig07_bic_icl_v3.png", w = 12, h = 5.5)


# ── FIG 08 · Distribuição de max(τ̂_ih) — certeza de classificação ──────────

post_long <- do.call(rbind, lapply(VARS_V3, function(var) {
  res_v <- resultados_univ[[var]]
  data.frame(
    variavel = factor(var, levels = VARS_V3),
    max_tau  = res_v$max_post,
    limiar   = 1 / res_v$g_best,
    g        = res_v$g_best
  )
}))

# mediana e % > 0.9 por variável
post_sum <- do.call(rbind, lapply(VARS_V3, function(var) {
  sub <- post_long[post_long$variavel == var, ]
  data.frame(
    variavel = factor(var, levels = VARS_V3),
    mediana  = median(sub$max_tau),
    pct90    = mean(sub$max_tau > 0.9) * 100,
    limiar   = unique(sub$limiar)
  )
}))
post_sum$label <- sprintf("med=%.2f\n%.0f%%>0.9", post_sum$mediana, post_sum$pct90)

p08 <- ggplot(post_long, aes(max_tau)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "#2980b9", colour = "white", alpha = 0.65) +
  geom_density(colour = "#1a5276", linewidth = 0.9) +
  geom_vline(data = post_sum, aes(xintercept = mediana),
             colour = "#c0392b", linetype = "dashed", linewidth = 0.8) +
  geom_vline(data = post_sum, aes(xintercept = limiar),
             colour = "#e6a817", linetype = "dotted", linewidth = 0.8) +
  geom_text(data = post_sum,
            aes(x = 0.04, y = Inf, label = label),
            hjust = 0, vjust = 1.4, size = 2.9, colour = "grey25",
            lineheight = 1.2) +
  facet_wrap(~ variavel, ncol = 3, scales = "free_y") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  labs(
    x       = "max τ̂_ih  (certeza de pertencimento ao melhor componente)",
    y       = "Densidade",
    title   = "Certeza de classificação por variável — flexmix NB univariado",
    subtitle = "Vermelho tracejado = mediana  ·  Amarelo pontilhado = chance aleatória (1/g*)",
    caption  = "Valores próximos de 1 = beneficiário bem alocado  ·  próximos de 1/g = sobreposição entre grupos"
  ) +
  TEMA

salvar(p08, "fig08_posteriores_v3.png", w = 13, h = 7)


cat("\nFiguras salvas em:", DIR_FIGS_HTML, "\n")
cat("  fig01_sobredispersao.png\n")
cat("  fig02_pois_vs_nb.png\n")
cat("  fig03_residuos_nb.png\n")
cat("  fig04_gaussiana_mix.png\n")
cat("  fig05_bic_familias.png\n")
cat("  fig06_componentes_nb_v3.png\n")
cat("  fig07_bic_icl_v3.png\n")
cat("  fig08_posteriores_v3.png\n")


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
