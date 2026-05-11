# ============================================================================
# explorar_modelos.R
# Exploração interativa — Misturas NB, Gaussiana e K-means
# Katy Garcia de Freitas · 2026
#
# USO:
#   - Abra no RStudio e rode com Source (Ctrl+Shift+S) ou por partes
#   - Os melhores fits são salvos em rds/ para uso no Rmd
#   - Ajuste os parâmetros na Seção 1 conforme necessário
# ============================================================================

# ── Seção 1 · Parâmetros ─────────────────────────────────────────────────────

CAMINHO_CSV <- "../../dados/dados_simulados.csv"
VARS_MODELO <- c("consultas", "ps", "exames")
VARS_EXTRA  <- c("consultas", "ps", "exames", "internacoes")           # + G3-discriminador
VARS_TUDO   <- c("consultas", "ps", "exames", "internacoes", "terapias") # + G2-discriminador
G_MAX       <- 6L      # número máximo de componentes a testar
SEED        <- 42L
N_INIT_NB   <- 8L      # reinicializações (1ª = K-means, demais = aleatórias)
MAX_ITER_NB <- 200L    # iterações máximas do EM
TOL_NB      <- 1e-5    # critério de convergência
DIR_RDS     <- "rds"   # pasta para salvar os objetos ajustados

# ── Seção 2 · Pacotes ─────────────────────────────────────────────────────────

pkgs <- c("dplyr", "tidyr", "ggplot2", "scales",
          "mclust", "cluster", "clue", "aricode", "kableExtra")

novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Instalando: ", paste(novos, collapse = ", "))
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(mclust); library(cluster); library(clue); library(aricode)
})

# ── Seção 3 · Dados ───────────────────────────────────────────────────────────

df_raw <- read.csv(CAMINHO_CSV, stringsAsFactors = FALSE, encoding = "UTF-8")

candidatos <- c("grupo_latente", "grupo", "grupo_real", "grupo_idx")
col_grupo  <- candidatos[candidatos %in% names(df_raw)][1]
if (is.na(col_grupo)) stop("Coluna de grupo não encontrada.")

df <- df_raw %>%
  rename_at(col_grupo, ~ "grupo_real") %>%
  mutate(grupo_real = factor(grupo_real))

df <- df[complete.cases(df[, VARS_MODELO]), ]

G_REAL         <- nlevels(df$grupo_real)
NIVEIS_REAL    <- levels(df$grupo_real)
grupo_real_int <- as.integer(df$grupo_real)

X_raw <- df %>% select(one_of(VARS_MODELO)) %>% as.matrix()
X_log <- log1p(X_raw)
X_std <- scale(X_log)

cat(sprintf("\n── Dados carregados: %d obs · %d vars · %d grupos ──\n",
            nrow(df), length(VARS_MODELO), G_REAL))
cat("Grupos:", paste(NIVEIS_REAL, collapse = " | "), "\n\n")

cat("Média por grupo e variável (referência):\n")
for (v in VARS_MODELO) {
  cat(sprintf("  %-12s", v))
  tapply(df[[v]], df$grupo_real, mean) %>%
    round(2) %>%
    { cat(paste(names(.), "=", ., collapse = "  "), "\n") }
}

# ── Seção 4 · Funções ─────────────────────────────────────────────────────────

# Atualiza theta_hj por otimização escalar no log-espaço
.update_theta <- function(y_j, tau_h, mu_hj) {
  neg_ll <- function(lt) {
    th <- exp(lt)
    -sum(tau_h * dnbinom(y_j, size = th, mu = mu_hj, log = TRUE))
  }
  lt_opt <- tryCatch(
    optimize(neg_ll, interval = c(-4, 8))$minimum,
    error = function(e) 0
  )
  exp(lt_opt)
}

# EM para mistura de Binomial Negativa (independência condicional)
#
# Inicialização:
#   init = 1 → K-means em log(1+Y): ponto de partida próximo da solução real
#   init > 1 → τ aleatório suavizado: explora outros máximos locais
nb_mix_em <- function(Y, g,
                      max_iter = MAX_ITER_NB,
                      tol      = TOL_NB,
                      n_init   = N_INIT_NB,
                      seed     = SEED,
                      verbose  = FALSE) {
  n <- nrow(Y); p <- ncol(Y)
  Y <- as.matrix(Y)
  best_loglik <- -Inf
  best_fit    <- NULL

  for (init in seq_len(n_init)) {

    # ── Inicialização ──────────────────────────────────────────────────────
    if (init == 1) {
      # K-means no espaço log: boa aproximação inicial sem custo extra
      set.seed(seed)
      km  <- kmeans(log1p(Y), centers = g, nstart = 10, iter.max = 100)
      tau <- matrix(0.02 / (g - 1 + 1e-9), n, g)   # massa pequena nos outros
      for (i in seq_len(n)) tau[i, km$cluster[i]] <- 0.98
    } else {
      set.seed(seed + init)
      tau <- matrix(runif(n * g) + 0.01, n, g)
    }
    tau   <- tau / rowSums(tau)
    theta <- matrix(2, g, p)
    loglik_old <- -Inf

    for (iter in seq_len(max_iter)) {
      # M-step
      pi_h <- colMeans(tau)
      n_h  <- colSums(tau)
      mu   <- t(tau) %*% Y / n_h
      mu   <- pmax(mu, 1e-6)
      for (h in seq_len(g))
        for (j in seq_len(p))
          theta[h, j] <- .update_theta(Y[, j], tau[, h], mu[h, j])
      theta <- pmax(theta, 0.05)

      # E-step (log-sum-exp)
      log_f <- matrix(0, n, g)
      for (h in seq_len(g)) {
        mu_mat    <- matrix(mu[h, ],    n, p, byrow = TRUE)
        theta_mat <- matrix(theta[h, ], n, p, byrow = TRUE)
        log_f[, h] <- log(pi_h[h]) +
          rowSums(dnbinom(Y, size = theta_mat, mu = mu_mat, log = TRUE))
      }
      log_max  <- apply(log_f, 1, max)
      log_norm <- log(rowSums(exp(log_f - log_max))) + log_max
      loglik   <- sum(log_norm)
      tau_new  <- exp(log_f - log_norm)
      tau_new  <- tau_new / rowSums(tau_new)

      if (iter > 1 && abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik
      tau        <- tau_new
    }

    if (verbose)
      cat(sprintf("    init %d/%d · g=%d · iter=%d · loglik=%.1f\n",
                  init, n_init, g, iter, loglik))

    if (loglik > best_loglik) {
      best_loglik <- loglik
      k_params    <- 2L * g * p + (g - 1L)
      best_fit    <- list(
        g              = g,
        pi             = pi_h,
        mu             = mu,
        theta          = theta,
        tau            = tau_new,
        loglik         = loglik,
        n_iter         = iter,
        classification = apply(tau_new, 1, which.max),
        BIC            = -2 * loglik + k_params * log(n),
        AIC            = -2 * loglik + 2 * k_params
      )
    }
  }
  best_fit
}

# Pareamento ótimo de rótulos (algoritmo húngaro)
#
# solve_LSAP retorna σ onde σ[i]=j significa "grupo real i → cluster j".
# Precisamos do INVERSO: "cluster j → grupo real i" = order(σ).
parear_otimo <- function(real, estimado) {
  g_r <- max(real); g_e <- max(estimado); g <- max(g_r, g_e)
  custo <- matrix(0, g, g)
  for (i in seq_len(g_r))
    for (j in seq_len(g_e))
      custo[i, j] <- -sum(real == i & estimado == j)
  sol        <- clue::solve_LSAP(custo + abs(min(custo)))
  mapeamento <- order(as.integer(sol))   # ← inverso da permutação
  mapeamento[estimado]
}

# Métricas de recuperação
metricas <- function(real_int, estimado_int, nome) {
  est_par <- parear_otimo(real_int, estimado_int)
  pur <- mean(sapply(unique(estimado_int), function(k) {
    idx <- estimado_int == k
    max(table(real_int[idx])) / sum(idx)
  }))
  list(
    metodo     = nome,
    n_clusters = length(unique(estimado_int)),
    ari        = mclust::adjustedRandIndex(real_int, estimado_int),
    nmi        = aricode::NMI(real_int, estimado_int),
    acuracia   = mean(est_par == real_int),
    pureza     = pur,
    est_par    = est_par
  )
}

# ── Seção 5 · Modelo 1 — Mistura Gaussiana ────────────────────────────────────

cat("\n══════════════════════════════════════════\n")
cat("MODELO 1 · Mistura Gaussiana (mclust)\n")
cat("══════════════════════════════════════════\n")

set.seed(SEED)
m1_bic <- mclust::mclustBIC(X_std, G = 1:G_MAX, verbose = FALSE)
m1_fit <- mclust::Mclust(X_std, G = 1:G_MAX, x = m1_bic, verbose = FALSE)
M1_G   <- m1_fit$G

cat(sprintf("Selecionado: %s · g = %d\n", m1_fit$modelName, M1_G))
cat("BIC por g:\n")
print(round(apply(m1_bic, 1, max, na.rm = TRUE), 1))

r1 <- metricas(grupo_real_int, m1_fit$classification, "Gaussiana")
cat(sprintf("\nARI = %.3f · NMI = %.3f · Acurácia = %.3f · Pureza = %.3f\n",
            r1$ari, r1$nmi, r1$acuracia, r1$pureza))

cat("\nDistribuição dos clusters:\n")
print(table(Cluster = m1_fit$classification, Real = df$grupo_real))

# ── Seção 6 · Modelo 2 — Mistura NB ──────────────────────────────────────────

cat("\n══════════════════════════════════════════\n")
cat("MODELO 2 · Mistura Binomial Negativa (EM)\n")
cat(sprintf("n_init=%d · max_iter=%d · tol=%.0e\n", N_INIT_NB, MAX_ITER_NB, TOL_NB))
cat("══════════════════════════════════════════\n")

Y_counts <- X_raw
m2_fits  <- vector("list", G_MAX)
m2_bic   <- numeric(G_MAX)
m2_aic   <- numeric(G_MAX)

for (g in seq_len(G_MAX)) {
  cat(sprintf("  g = %d ...", g))
  t0 <- proc.time()
  m2_fits[[g]] <- nb_mix_em(Y_counts, g = g, verbose = FALSE)
  dt <- (proc.time() - t0)["elapsed"]
  m2_bic[g]    <- m2_fits[[g]]$BIC
  m2_aic[g]    <- m2_fits[[g]]$AIC
  cat(sprintf(" BIC = %8.1f · AIC = %8.1f · iter = %d · %.1fs\n",
              m2_bic[g], m2_aic[g], m2_fits[[g]]$n_iter, dt))
}

M2_G   <- which.min(m2_bic)
m2_opt <- m2_fits[[M2_G]]

cat(sprintf("\nSelecionado: g = %d · loglik = %.1f · BIC = %.1f\n",
            M2_G, m2_opt$loglik, m2_opt$BIC))

cat("\nParâmetros estimados (médias μ):\n")
mu_tab <- as.data.frame(round(m2_opt$mu, 2))
colnames(mu_tab) <- VARS_MODELO
rownames(mu_tab) <- paste0("G", seq_len(M2_G))
print(mu_tab)

cat("\nParâmetros estimados (dispersão θ):\n")
th_tab <- as.data.frame(round(m2_opt$theta, 3))
colnames(th_tab) <- VARS_MODELO
rownames(th_tab) <- paste0("G", seq_len(M2_G))
print(th_tab)

cat("\nProporções das componentes:\n")
print(round(setNames(m2_opt$pi, paste0("G", seq_len(M2_G))), 3))

cat("\nIncerteza de classificação:\n")
incert <- 1 - apply(m2_opt$tau, 1, max)
cat(sprintf("  Média = %.3f · Mediana = %.3f · P95 = %.3f\n",
            mean(incert), median(incert), quantile(incert, 0.95)))

r2 <- metricas(grupo_real_int, m2_opt$classification, "Mistura NB")
cat(sprintf("\nARI = %.3f · NMI = %.3f · Acurácia = %.3f · Pureza = %.3f\n",
            r2$ari, r2$nmi, r2$acuracia, r2$pureza))

cat("\nDistribuição dos clusters:\n")
print(table(Cluster = m2_opt$classification, Real = df$grupo_real))

# ── g = G_REAL forçado (comparação direta) ───────────────────────────────────
if (M2_G != G_REAL) {
  cat(sprintf("\n── NB com g = %d forçado (número real de grupos) ──\n", G_REAL))
  t0 <- proc.time()
  m2_g4 <- nb_mix_em(Y_counts, g = G_REAL,
                     n_init   = N_INIT_NB,
                     max_iter = MAX_ITER_NB,
                     tol      = TOL_NB,
                     verbose  = FALSE)
  dt <- (proc.time() - t0)["elapsed"]
  cat(sprintf("loglik = %.1f · BIC = %.1f · iter = %d · %.1fs\n",
              m2_g4$loglik, m2_g4$BIC, m2_g4$n_iter, dt))

  mu_g4 <- as.data.frame(round(m2_g4$mu, 2))
  colnames(mu_g4) <- VARS_MODELO; rownames(mu_g4) <- paste0("G", seq_len(G_REAL))
  cat("\nMédias μ (g=4 forçado):\n"); print(mu_g4)

  r2_g4 <- metricas(grupo_real_int, m2_g4$classification, sprintf("Mistura NB (g=%d forçado)", G_REAL))
  cat(sprintf("\nARI = %.3f · NMI = %.3f · Acurácia = %.3f · Pureza = %.3f\n",
              r2_g4$ari, r2_g4$nmi, r2_g4$acuracia, r2_g4$pureza))
  cat("\nDistribuição dos clusters (g=4 forçado):\n")
  print(table(Cluster = m2_g4$classification, Real = df$grupo_real))

  # Guarda para o Rmd usar o melhor dos dois
  saveRDS(m2_g4, file.path(DIR_RDS, "m2_nb_g4.rds"))
} else {
  cat(sprintf("\n(BIC selecionou g = %d = número real de grupos)\n", M2_G))
}

# ── Seção 6b · NB com variáveis ampliadas (+ internacoes) ────────────────────
#
# internacoes separa fortemente G3 (agudo/hospitalar): ~80% de prevalência vs
# ~5% em G1. Incluir mesmo com muitos zeros: NB com mu pequeno e theta pequeno
# acomoda zeros naturalmente.

cat("\n══════════════════════════════════════════\n")
cat("MODELO 2b · Mistura NB com internacoes\n")
cat(sprintf("variáveis: %s\n", paste(VARS_EXTRA, collapse = ", ")))
cat("══════════════════════════════════════════\n")

vars_extra_ausentes <- VARS_EXTRA[!VARS_EXTRA %in% names(df)]
if (length(vars_extra_ausentes) > 0) {
  cat("ATENÇÃO: variáveis ausentes:", paste(vars_extra_ausentes, collapse = ", "), "\n")
  cat("Seção 6b ignorada.\n")
  m2b_ok <- FALSE
} else {
  m2b_ok    <- TRUE
  Y_extra   <- df %>% select(one_of(VARS_EXTRA)) %>% as.matrix()

  m2b_fits  <- vector("list", G_MAX)
  m2b_bic   <- numeric(G_MAX)

  for (g in seq_len(G_MAX)) {
    cat(sprintf("  g = %d ...", g))
    t0 <- proc.time()
    m2b_fits[[g]] <- nb_mix_em(Y_extra, g = g, verbose = FALSE)
    dt <- (proc.time() - t0)["elapsed"]
    m2b_bic[g]    <- m2b_fits[[g]]$BIC
    cat(sprintf(" BIC = %8.1f · iter = %d · %.1fs\n",
                m2b_bic[g], m2b_fits[[g]]$n_iter, dt))
  }

  M2b_G    <- which.min(m2b_bic)
  m2b_opt  <- m2b_fits[[M2b_G]]

  cat(sprintf("\nSelecionado: g = %d · BIC = %.1f\n", M2b_G, m2b_opt$BIC))

  mu_b <- as.data.frame(round(m2b_opt$mu, 2))
  colnames(mu_b) <- VARS_EXTRA; rownames(mu_b) <- paste0("G", seq_len(M2b_G))
  cat("\nMédias μ:\n"); print(mu_b)

  r2b <- metricas(grupo_real_int, m2b_opt$classification,
                  sprintf("Mistura NB (+internacoes, g=%d)", M2b_G))
  cat(sprintf("\nARI = %.3f · NMI = %.3f · Acurácia = %.3f · Pureza = %.3f\n",
              r2b$ari, r2b$nmi, r2b$acuracia, r2b$pureza))

  cat("\nDistribuição dos clusters:\n")
  print(table(Cluster = m2b_opt$classification, Real = df$grupo_real))

  saveRDS(m2b_opt,  file.path(DIR_RDS, "m2b_nb_extra.rds"))
  saveRDS(list(bic = m2b_bic), file.path(DIR_RDS, "m2b_criterios.rds"))
}

# ── Seção 6c · NB com todas as contagens (+ internacoes + terapias) ──────────
#
# terapias discrimina G2 (ambulatorial coordenado): alta prevalência vs G1/G3.
# internacoes discrimina G3 (agudo/hospitalar).
# Usando todas as 5 variáveis de contagem relevantes.

cat("\n══════════════════════════════════════════\n")
cat("MODELO 2c · Mistura NB com internacoes + terapias\n")
cat(sprintf("variáveis: %s\n", paste(VARS_TUDO, collapse = ", ")))
cat("══════════════════════════════════════════\n")

vars_tudo_ausentes <- VARS_TUDO[!VARS_TUDO %in% names(df)]
if (length(vars_tudo_ausentes) > 0) {
  cat("ATENÇÃO: variáveis ausentes:", paste(vars_tudo_ausentes, collapse = ", "), "\n")
  cat("Seção 6c ignorada.\n")
  m2c_ok <- FALSE
} else {
  m2c_ok   <- TRUE
  Y_tudo   <- df %>% select(one_of(VARS_TUDO)) %>% as.matrix()

  m2c_fits <- vector("list", G_MAX)
  m2c_bic  <- numeric(G_MAX)

  for (g in seq_len(G_MAX)) {
    cat(sprintf("  g = %d ...", g))
    t0 <- proc.time()
    m2c_fits[[g]] <- nb_mix_em(Y_tudo, g = g, verbose = FALSE)
    dt <- (proc.time() - t0)["elapsed"]
    m2c_bic[g]    <- m2c_fits[[g]]$BIC
    cat(sprintf(" BIC = %8.1f · iter = %d · %.1fs\n",
                m2c_bic[g], m2c_fits[[g]]$n_iter, dt))
  }

  M2c_G   <- which.min(m2c_bic)
  m2c_opt <- m2c_fits[[M2c_G]]

  cat(sprintf("\nSelecionado: g = %d · BIC = %.1f\n", M2c_G, m2c_opt$BIC))

  mu_c <- as.data.frame(round(m2c_opt$mu, 2))
  colnames(mu_c) <- VARS_TUDO; rownames(mu_c) <- paste0("G", seq_len(M2c_G))
  cat("\nMédias μ:\n"); print(mu_c)

  r2c <- metricas(grupo_real_int, m2c_opt$classification,
                  sprintf("Mistura NB (5 vars, g=%d)", M2c_G))
  cat(sprintf("\nARI = %.3f · NMI = %.3f · Acurácia = %.3f · Pureza = %.3f\n",
              r2c$ari, r2c$nmi, r2c$acuracia, r2c$pureza))

  cat("\nDistribuição dos clusters:\n")
  print(table(Cluster = m2c_opt$classification, Real = df$grupo_real))

  saveRDS(m2c_opt, file.path(DIR_RDS, "m2c_nb_tudo.rds"))
  saveRDS(list(bic = m2c_bic), file.path(DIR_RDS, "m2c_criterios.rds"))
}

# ── Seção 7 · Modelo 3 — K-means ─────────────────────────────────────────────

cat("\n══════════════════════════════════════════\n")
cat("MODELO 3 · K-means\n")
cat("══════════════════════════════════════════\n")

set.seed(SEED)

# Curva do cotovelo
cat("  Curva do cotovelo...")
wcss <- sapply(1:G_MAX, function(k) {
  kmeans(X_std, centers = k, nstart = 20, iter.max = 200)$tot.withinss
})
cat(" OK\n")

# Silhouette (dist() é pesado — salvo para reusar no Rmd)
cat("  Silhouette (dist)...")
dist_std <- dist(X_std)
sil_vals <- sapply(2:G_MAX, function(k) {
  km  <- kmeans(X_std, centers = k, nstart = 20, iter.max = 200)
  sil <- cluster::silhouette(km$cluster, dist_std)
  mean(sil[, "sil_width"])
})
cat(" OK\n")

# Ajuste final com k = G_REAL
M3_K  <- G_REAL
m3_km <- kmeans(X_std, centers = M3_K, nstart = 50, iter.max = 300)

r3 <- metricas(grupo_real_int, m3_km$cluster, "K-means")
cat(sprintf("k = %d · WCSS = %.1f\n", M3_K, m3_km$tot.withinss))
cat(sprintf("ARI = %.3f · NMI = %.3f · Acurácia = %.3f · Pureza = %.3f\n",
            r3$ari, r3$nmi, r3$acuracia, r3$pureza))

cat("\nDistribuição dos clusters:\n")
print(table(Cluster = m3_km$cluster, Real = df$grupo_real))

# ── Seção 8 · Comparativo final ───────────────────────────────────────────────

cat("\n══════════════════════════════════════════\n")
cat("COMPARATIVO FINAL\n")
cat("══════════════════════════════════════════\n")

res_lista <- list(r1, r2, r3)
if (exists("r2b")) res_lista <- c(res_lista, list(r2b))
if (exists("r2c")) res_lista <- c(res_lista, list(r2c))

tab <- data.frame(
  Método   = sapply(res_lista, `[[`, "metodo"),
  g        = sapply(res_lista, `[[`, "n_clusters"),
  ARI      = round(sapply(res_lista, `[[`, "ari"),      3),
  NMI      = round(sapply(res_lista, `[[`, "nmi"),      3),
  Acuracia = round(sapply(res_lista, `[[`, "acuracia"), 3),
  Pureza   = round(sapply(res_lista, `[[`, "pureza"),   3)
)
print(tab)

melhor <- tab$Método[which.max(tab$ARI)]
cat(sprintf("\nMelhor ARI: %s (%.3f)\n", melhor, max(tab$ARI)))

# Recomendação: qual variante do modelo NB usar no Rmd?
cat("\n── Recomendação para o Rmd ──────────────────────────────────\n")
nb_rows <- grep("NB", tab$Método)
if (length(nb_rows) > 0) {
  melhor_nb <- nb_rows[which.max(tab$ARI[nb_rows])]
  cat(sprintf("Melhor modelo NB: %s\n", tab$Método[melhor_nb]))
  cat(sprintf("  ARI = %.3f · Acurácia = %.3f\n",
              tab$ARI[melhor_nb], tab$Acuracia[melhor_nb]))
}

# ── Seção 9 · Salvar RDS para o Rmd ──────────────────────────────────────────

if (!dir.exists(DIR_RDS)) dir.create(DIR_RDS, recursive = TRUE)

saveRDS(m1_fit,  file.path(DIR_RDS, "m1_gaussiana.rds"))
saveRDS(m2_opt,  file.path(DIR_RDS, "m2_nb_otimo.rds"))
saveRDS(m2_fits, file.path(DIR_RDS, "m2_nb_todos.rds"))
saveRDS(list(bic = m2_bic, aic = m2_aic),
        file.path(DIR_RDS, "m2_criterios.rds"))
saveRDS(m3_km,   file.path(DIR_RDS, "m3_kmeans.rds"))
saveRDS(list(wcss = wcss, sil_vals = sil_vals),
        file.path(DIR_RDS, "m3_selecao.rds"))

cat(sprintf("\nObjetos salvos em '%s/':\n", DIR_RDS))
cat("  m1_gaussiana.rds\n")
cat("  m2_nb_otimo.rds · m2_nb_todos.rds · m2_criterios.rds\n")
if (exists("m2b_ok") && m2b_ok)
  cat("  m2b_nb_extra.rds · m2b_criterios.rds\n")
if (exists("m2c_ok") && m2c_ok)
  cat("  m2c_nb_tudo.rds  · m2c_criterios.rds\n")
cat("  m3_kmeans.rds   · m3_selecao.rds\n")
cat("\nPróximo passo:\n")
cat("  1. Veja qual NB teve melhor resultado (seção 'Recomendação' acima)\n")
cat("  2. Ajuste analises.Rmd para usar o melhor RDS\n")
cat("  3. knit analises.Rmd\n")
