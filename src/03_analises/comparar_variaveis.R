# =============================================================================
# comparar_variaveis.R
# Comparação sistemática: 5 conjuntos de variáveis × 2 métodos
#
# Conjuntos testados (adição progressiva):
#   V1  consultas + ps + exames                          (baseline)
#   V2  + internacoes                                    (discrimina G3)
#   V3  + terapias                                       (discrimina G4)
#   V4  + prestadores_distintos                          (sinal de fragmentação)
#   V5  + meses_utilizacao                               (padrão agudo vs. crônico)
#
# Métodos (dados de contagem → escolha motivada pela natureza discreta):
#   M2  Mistura NB  — EM do zero, distribuição NB por componente e variável
#   M3  K-means     — baseline geométrico livre de distribuição, log1p + scale
#
# Critério de seleção de g: BIC (M2) · Silhouette (M3)
# Métricas: ARI · NMI · Acurácia (pareamento húngaro) · Pureza · g · tempo(s)
#
# Nota sobre meses_utilizacao (V5):
#   Gerada via Beta-Binomial (limite 1–12), distribuição bimodal em U.
#   O NB aproxima de forma imperfeita — resultado deve ser interpretado
#   com cautela para M2 nesse conjunto.
#
# USO: Source (Ctrl+Shift+S) no RStudio
#      Resultados salvos em rds/comparacao_variaveis.rds
# =============================================================================


# ── 1. Parâmetros ─────────────────────────────────────────────────────────────

CAMINHO_CSV <- "../../dados/dados_simulados.csv"
G_MAX       <- 6L
SEED        <- 42L
N_INIT_NB   <- 8L
MAX_ITER_NB <- 200L
TOL_NB      <- 1e-5
DIR_RDS     <- "rds"

CONJUNTOS <- list(
  V1 = c("consultas", "ps", "exames"),
  V2 = c("consultas", "ps", "exames", "internacoes"),
  V3 = c("consultas", "ps", "exames", "internacoes", "terapias"),
  V4 = c("consultas", "ps", "exames", "internacoes", "terapias",
          "prestadores_distintos"),
  V5 = c("consultas", "ps", "exames", "internacoes", "terapias",
          "prestadores_distintos", "meses_utilizacao")
)

LABELS_CONJUNTOS <- c(
  V1 = "Baseline (3 vars)",
  V2 = "+ internacoes",
  V3 = "+ terapias",
  V4 = "+ prestadores",
  V5 = "+ meses_utilizacao *"
)


# ── 2. Pacotes ────────────────────────────────────────────────────────────────

pkgs <- c("dplyr", "cluster", "clue", "aricode")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Instalando: ", paste(novos, collapse = ", "))
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(dplyr); library(cluster); library(clue); library(aricode)
})


# ── 3. Dados ──────────────────────────────────────────────────────────────────

df_raw <- read.csv(CAMINHO_CSV, stringsAsFactors = FALSE, encoding = "UTF-8")

candidatos <- c("grupo_latente", "grupo", "grupo_real", "grupo_idx")
col_grupo  <- candidatos[candidatos %in% names(df_raw)][1]
if (is.na(col_grupo)) stop("Coluna de grupo não encontrada.")

df <- df_raw
names(df)[names(df) == col_grupo] <- "grupo_real"
df$grupo_real   <- factor(df$grupo_real)
G_REAL          <- nlevels(df$grupo_real)
NIVEIS_REAL     <- levels(df$grupo_real)
grupo_real_int  <- as.integer(df$grupo_real)

cat(sprintf("\nDados: %d obs · %d grupos reais\n", nrow(df), G_REAL))
cat("Grupos:", paste(NIVEIS_REAL, collapse = " | "), "\n\n")


# ── 4. Funções ────────────────────────────────────────────────────────────────

# Theta update (NB EM, log-espaço)
.update_theta <- function(y_j, tau_h, mu_hj) {
  neg_ll <- function(lt) {
    th <- exp(lt)
    -sum(tau_h * dnbinom(y_j, size = th, mu = mu_hj, log = TRUE))
  }
  exp(tryCatch(optimize(neg_ll, interval = c(-4, 8))$minimum, error = function(e) 0))
}

# EM para mistura NB (independência condicional)
# Init 1 = K-means em log1p(Y); demais = aleatório
nb_mix_em <- function(Y, g,
                      max_iter = MAX_ITER_NB, tol = TOL_NB,
                      n_init = N_INIT_NB, seed = SEED) {
  n <- nrow(Y); p <- ncol(Y)
  Y <- as.matrix(Y)
  best_loglik <- -Inf; best_fit <- NULL

  for (init in seq_len(n_init)) {
    if (init == 1) {
      set.seed(seed)
      km  <- kmeans(log1p(Y), centers = g, nstart = 10, iter.max = 100)
      tau <- matrix(0.02 / max(g - 1, 1), n, g)
      for (i in seq_len(n)) tau[i, km$cluster[i]] <- 0.98
    } else {
      set.seed(seed + init)
      tau <- matrix(runif(n * g) + 0.01, n, g)
    }
    tau   <- tau / rowSums(tau)
    theta <- matrix(2, g, p)
    loglik_old <- -Inf

    for (iter in seq_len(max_iter)) {
      pi_h <- colMeans(tau); n_h <- colSums(tau)
      mu   <- pmax(t(tau) %*% Y / n_h, 1e-6)
      for (h in seq_len(g))
        for (j in seq_len(p))
          theta[h, j] <- .update_theta(Y[, j], tau[, h], mu[h, j])
      theta <- pmax(theta, 0.05)

      log_f <- matrix(0, n, g)
      for (h in seq_len(g)) {
        log_f[, h] <- log(pi_h[h]) +
          rowSums(dnbinom(Y,
                          size = matrix(theta[h, ], n, p, byrow = TRUE),
                          mu   = matrix(mu[h, ],    n, p, byrow = TRUE),
                          log  = TRUE))
      }
      log_max  <- apply(log_f, 1, max)
      log_norm <- log(rowSums(exp(log_f - log_max))) + log_max
      loglik   <- sum(log_norm)
      tau_new  <- exp(log_f - log_norm); tau_new <- tau_new / rowSums(tau_new)

      if (iter > 1 && abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik; tau <- tau_new
    }

    if (loglik > best_loglik) {
      best_loglik <- loglik
      k_params    <- 2L * g * p + (g - 1L)
      best_fit    <- list(
        g = g, pi = pi_h, mu = mu, theta = theta, tau = tau_new,
        loglik = loglik, n_iter = iter,
        classification = apply(tau_new, 1, which.max),
        BIC = -2 * loglik + k_params * log(n),
        AIC = -2 * loglik + 2 * k_params
      )
    }
  }
  best_fit
}

# Pareamento ótimo (algoritmo húngaro) — usa inverso da permutação
parear_otimo <- function(real, estimado) {
  g <- max(max(real), max(estimado))
  custo <- matrix(0, g, g)
  for (i in seq_len(max(real)))
    for (j in seq_len(max(estimado)))
      custo[i, j] <- -sum(real == i & estimado == j)
  sol <- clue::solve_LSAP(custo + abs(min(custo)))
  order(as.integer(sol))[estimado]
}

# Métricas de recuperação
calc_metricas <- function(real_int, estimado_int) {
  est_par <- parear_otimo(real_int, estimado_int)
  pur <- mean(sapply(unique(estimado_int), function(k) {
    idx <- estimado_int == k
    max(table(real_int[idx])) / sum(idx)
  }))
  list(
    g        = length(unique(estimado_int)),
    ari      = round(aricode::ARI(real_int, estimado_int),              3),
    nmi      = round(aricode::NMI(real_int, estimado_int),              3),
    acuracia = round(mean(est_par == real_int),                         3),
    pureza   = round(pur,                                               3)
  )
}


# ── 5. Loop principal ─────────────────────────────────────────────────────────

cat("Iniciando comparação: ", length(CONJUNTOS), "conjuntos × 2 métodos\n")
cat(strrep("═", 60), "\n\n")

resultados <- list()   # acumula uma linha por (conjunto × método)
fits_nb    <- list()   # guarda os fits NB selecionados para reuso

for (nome_v in names(CONJUNTOS)) {

  vars <- CONJUNTOS[[nome_v]]
  label <- LABELS_CONJUNTOS[nome_v]

  # Verifica disponibilidade
  ausentes <- vars[!vars %in% names(df)]
  if (length(ausentes) > 0) {
    cat(sprintf("  [%s] Variáveis ausentes: %s — ignorado\n",
                nome_v, paste(ausentes, collapse = ", ")))
    next
  }

  cat(sprintf("── %s · %s ──\n", nome_v, label))

  Y_raw <- as.matrix(df[, vars])
  Y_std <- scale(log1p(Y_raw))

  # ── M2 · Mistura NB ───────────────────────────────────────────────────────
  cat("  M2 NB       ... ")
  m2_fits  <- vector("list", G_MAX)
  m2_bic   <- numeric(G_MAX)
  t0       <- proc.time()

  for (g in seq_len(G_MAX)) {
    m2_fits[[g]] <- nb_mix_em(Y_raw, g = g)
    m2_bic[g]    <- m2_fits[[g]]$BIC
    cat(sprintf("[g=%d:%.0f]", g, m2_bic[g]))
  }

  M2_G    <- which.min(m2_bic)
  m2_opt  <- m2_fits[[M2_G]]
  dt      <- (proc.time() - t0)["elapsed"]
  m2_met  <- calc_metricas(grupo_real_int, m2_opt$classification)

  cat(sprintf("\n         g=%d · ARI=%.3f · Acurácia=%.3f · %.1fs\n",
              m2_met$g, m2_met$ari, m2_met$acuracia, dt))

  fits_nb[[nome_v]] <- m2_opt   # guarda para RDS

  resultados[[length(resultados) + 1]] <- data.frame(
    Conjunto = label, Método = "Mistura NB",
    g = m2_met$g, ARI = m2_met$ari, NMI = m2_met$nmi,
    Acuracia = m2_met$acuracia, Pureza = m2_met$pureza,
    Tempo_s = round(dt, 1), stringsAsFactors = FALSE
  )

  # ── M3 · K-means ─────────────────────────────────────────────────────────
  # Seleção por silhouette (g=2..G_MAX); ajuste final com g selecionado
  cat("  M3 K-means  ... ")
  t0      <- proc.time()
  set.seed(SEED)
  dist_std <- dist(Y_std)
  sil_vals <- sapply(2:G_MAX, function(k) {
    km  <- kmeans(Y_std, centers = k, nstart = 15, iter.max = 200)
    mean(cluster::silhouette(km$cluster, dist_std)[, "sil_width"])
  })
  M3_K    <- (2:G_MAX)[which.max(sil_vals)]
  m3_km   <- kmeans(Y_std, centers = M3_K, nstart = 50, iter.max = 300)
  dt      <- (proc.time() - t0)["elapsed"]
  m3_met  <- calc_metricas(grupo_real_int, m3_km$cluster)

  cat(sprintf("g=%d · ARI=%.3f · Acurácia=%.3f · %.1fs\n",
              m3_met$g, m3_met$ari, m3_met$acuracia, dt))

  resultados[[length(resultados) + 1]] <- data.frame(
    Conjunto = label, Método = "K-means",
    g = m3_met$g, ARI = m3_met$ari, NMI = m3_met$nmi,
    Acuracia = m3_met$acuracia, Pureza = m3_met$pureza,
    Tempo_s = round(dt, 1), stringsAsFactors = FALSE
  )

  cat("\n")
}


# ── 6. Tabela comparativa ─────────────────────────────────────────────────────

tab <- do.call(rbind, resultados)
rownames(tab) <- NULL

cat("\n", strrep("═", 60), "\n")
cat("TABELA COMPARATIVA FINAL\n")
cat(strrep("═", 60), "\n\n")
print(tab, row.names = FALSE)

# Melhor por método (ARI)
cat("\n── Melhor conjunto por método (maior ARI) ───────────────\n")
for (met in unique(tab$Método)) {
  sub   <- tab[tab$Método == met, ]
  best  <- sub[which.max(sub$ARI), ]
  cat(sprintf("  %-14s → %s  (ARI=%.3f · Acurácia=%.3f · g=%d)\n",
              met, best$Conjunto, best$ARI, best$Acuracia, best$g))
}

# Melhor geral (ARI)
cat("\n── Melhor resultado geral ───────────────────────────────\n")
best_all <- tab[which.max(tab$ARI), ]
cat(sprintf("  %s + %s  →  ARI=%.3f · NMI=%.3f · Acurácia=%.3f · Pureza=%.3f\n\n",
            best_all$Método, best_all$Conjunto,
            best_all$ARI, best_all$NMI, best_all$Acuracia, best_all$Pureza))

cat("* V5 inclui meses_utilizacao (Beta-Binomial): NB aproxima com ressalva.\n")


# ── 7. Salvar RDS ─────────────────────────────────────────────────────────────

if (!dir.exists(DIR_RDS)) dir.create(DIR_RDS, recursive = TRUE)

saveRDS(tab,     file.path(DIR_RDS, "comparacao_variaveis.rds"))
saveRDS(fits_nb, file.path(DIR_RDS, "nb_fits_por_conjunto.rds"))

cat(sprintf("\nSalvo em '%s/':\n", DIR_RDS))
cat("  comparacao_variaveis.rds  — tabela completa\n")
cat("  nb_fits_por_conjunto.rds  — fits NB de cada conjunto\n")
cat("\nPróximo passo: usar comparacao_variaveis.rds no analises.Rmd\n")
