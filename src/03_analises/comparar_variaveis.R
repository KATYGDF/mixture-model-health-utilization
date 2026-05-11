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
#   M1  Mistura NB  — EM do zero, distribuição NB por componente e variável
#   M2  K-means     — baseline geométrico livre de distribuição, log1p + scale
#
# Critério de seleção de g: BIC (M1) · Silhouette (M2)
# Métricas: ARI · NMI · Acurácia (pareamento húngaro) · Pureza · g · tempo(s)
#
# Nota sobre meses_utilizacao (V5):
#   Gerada via Beta-Binomial (limite 1–12), distribuição bimodal em U.
#   O NB aproxima de forma imperfeita — resultado deve ser interpretado
#   com cautela para M1 nesse conjunto.
#
# USO: Source (Ctrl+Shift+S) no RStudio
#      Resultados salvos em rds/  ·  Figuras salvas em figs/
# =============================================================================


# ── 1. Parâmetros ─────────────────────────────────────────────────────────────

CAMINHO_CSV <- "../../dados/dados_simulados.csv"
G_MAX       <- 6L
SEED        <- 42L
N_INIT_NB   <- 8L
MAX_ITER_NB <- 200L
TOL_NB      <- 1e-5
DIR_RDS     <- "rds"
DIR_FIGS    <- "figs"

CONJUNTOS <- list(
  V1 = c("consultas", "ps", "exames"),
  V2 = c("consultas", "ps", "exames", "internacoes"),
  V3 = c("consultas", "ps", "exames", "internacoes", "terapias")
  # V4 e V5 desativados: não melhoram o modelo (confirmado empiricamente)
  # V4 = c("consultas", "ps", "exames", "internacoes", "terapias",
  #         "prestadores_distintos"),
  # V5 = c("consultas", "ps", "exames", "internacoes", "terapias",
  #         "prestadores_distintos", "meses_utilizacao")
)

LABELS_CONJUNTOS <- c(
  V1 = "Baseline (3 vars)",
  V2 = "+ internacoes",
  V3 = "+ terapias"
  # V4 = "+ prestadores",
  # V5 = "+ meses_utilizacao *"
)


# ── 2. Pacotes ────────────────────────────────────────────────────────────────

pkgs <- c("dplyr", "cluster", "clue", "aricode",
          "ggplot2", "patchwork", "scales", "tidyr")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Instalando: ", paste(novos, collapse = ", "))
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(dplyr); library(cluster); library(clue); library(aricode)
  library(ggplot2); library(patchwork); library(scales); library(tidyr)
})


# ── 3. Dados ──────────────────────────────────────────────────────────────────

df_raw <- read.csv(CAMINHO_CSV, stringsAsFactors = FALSE, encoding = "UTF-8")

candidatos <- c("grupo_latente", "grupo", "grupo_real", "grupo_idx")
col_grupo  <- candidatos[candidatos %in% names(df_raw)][1]
if (is.na(col_grupo)) stop("Coluna de grupo não encontrada.")

df <- df_raw
names(df)[names(df) == col_grupo] <- "grupo_real"
df$grupo_real   <- factor(df$grupo_real,
                          levels = c("baixo_uso", "ambulatorial_coordenado",
                                     "agudo_hospitalar", "atipico"))
G_REAL          <- nlevels(df$grupo_real)
NIVEIS_REAL     <- levels(df$grupo_real)
grupo_real_int  <- as.integer(df$grupo_real)

# Legenda numérica: "1 = agudo_hospitalar  ·  2 = ambulatorial_coordenado  ·  ..."
LEGENDA_GRUPOS <- paste(
  sprintf("%d = %s", seq_along(NIVEIS_REAL), NIVEIS_REAL),
  collapse = "  ·  "
)

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
fits_nb    <- list()   # fits NB selecionados (melhor g por conjunto)
km_results <- list()   # fits K-means selecionados (melhor k por conjunto)
bic_curves <- list()   # BIC(g=1..G_MAX) por conjunto
sil_curves <- list()   # silhouette(k=2..G_MAX) por conjunto

for (nome_v in names(CONJUNTOS)) {

  vars  <- CONJUNTOS[[nome_v]]
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

  # ── M1 · Mistura NB ───────────────────────────────────────────────────────
  cat("  M1 NB       ... ")
  M1_fits  <- vector("list", G_MAX)
  M1_bic   <- numeric(G_MAX)
  t0       <- proc.time()

  for (g in seq_len(G_MAX)) {
    M1_fits[[g]] <- nb_mix_em(Y_raw, g = g)
    M1_bic[g]    <- M1_fits[[g]]$BIC
    cat(sprintf("[g=%d:%.0f]", g, M1_bic[g]))
  }

  M1_G    <- which.min(M1_bic)
  M1_opt  <- M1_fits[[M1_G]]
  dt      <- (proc.time() - t0)["elapsed"]
  M1_met  <- calc_metricas(grupo_real_int, M1_opt$classification)

  cat(sprintf("\n         g=%d · ARI=%.3f · Acurácia=%.3f · %.1fs\n",
              M1_met$g, M1_met$ari, M1_met$acuracia, dt))

  fits_nb[[nome_v]]    <- M1_opt
  bic_curves[[nome_v]] <- M1_bic

  resultados[[length(resultados) + 1]] <- data.frame(
    ConjID = nome_v, Conjunto = label, Método = "Mistura NB",
    g = M1_met$g, ARI = M1_met$ari, NMI = M1_met$nmi,
    Acuracia = M1_met$acuracia, Pureza = M1_met$pureza,
    Tempo_s = round(dt, 1), stringsAsFactors = FALSE
  )

  # ── M2 · K-means ─────────────────────────────────────────────────────────
  cat("  M2 K-means  ... ")
  t0      <- proc.time()
  set.seed(SEED)
  dist_std <- dist(Y_std)
  sil_vals <- sapply(2:G_MAX, function(k) {
    km  <- kmeans(Y_std, centers = k, nstart = 15, iter.max = 200)
    mean(cluster::silhouette(km$cluster, dist_std)[, "sil_width"])
  })
  M2_K    <- (2:G_MAX)[which.max(sil_vals)]
  M2_km   <- kmeans(Y_std, centers = M2_K, nstart = 50, iter.max = 300)
  dt      <- (proc.time() - t0)["elapsed"]
  M2_met  <- calc_metricas(grupo_real_int, M2_km$cluster)

  cat(sprintf("g=%d · ARI=%.3f · Acurácia=%.3f · %.1fs\n",
              M2_met$g, M2_met$ari, M2_met$acuracia, dt))

  sil_curves[[nome_v]] <- sil_vals
  km_results[[nome_v]] <- list(
    cluster = M2_km$cluster,
    k       = M2_K,
    Y_raw   = Y_raw,
    Y_std   = Y_std,
    vars    = vars
  )

  resultados[[length(resultados) + 1]] <- data.frame(
    ConjID = nome_v, Conjunto = label, Método = "K-means",
    g = M2_met$g, ARI = M2_met$ari, NMI = M2_met$nmi,
    Acuracia = M2_met$acuracia, Pureza = M2_met$pureza,
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
print(tab[, -1], row.names = FALSE)   # exclui ConjID da impressão

# Melhor por método (ARI)
cat("\n── Melhor conjunto por método (maior ARI) ───────────────\n")
for (met in unique(tab$Método)) {
  sub  <- tab[tab$Método == met, ]
  best <- sub[which.max(sub$ARI), ]
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

for (d in c(DIR_RDS, DIR_FIGS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

saveRDS(tab,        file.path(DIR_RDS, "comparacao_variaveis.rds"))
saveRDS(fits_nb,    file.path(DIR_RDS, "nb_fits_por_conjunto.rds"))
saveRDS(km_results, file.path(DIR_RDS, "km_results_por_conjunto.rds"))
saveRDS(bic_curves, file.path(DIR_RDS, "bic_curves.rds"))
saveRDS(sil_curves, file.path(DIR_RDS, "sil_curves.rds"))

cat(sprintf("\nSalvo em '%s/':\n", DIR_RDS))
cat("  comparacao_variaveis.rds  · nb_fits_por_conjunto.rds\n")
cat("  km_results_por_conjunto.rds · bic_curves.rds · sil_curves.rds\n")


# =============================================================================
# ── 8. Paleta, tema e helpers de visualização ─────────────────────────────────
# =============================================================================

# Paleta numérica (1–4) — pastéis consistentes com EDA
# Ordem da simulação: 1=baixo_uso  2=ambulatorial_coordenado  3=agudo_hospitalar  4=atipico
COR <- setNames(
  c("#7BAFD4", "#6EBF8B", "#E8B96B", "#D98585"),
  as.character(1:4)
)

TEMA <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "#f0f0f0", colour = NA),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    plot.caption     = element_text(colour = "grey55", size = 8, hjust = 0)
  )

# Helper: data frame de matriz de confusão com percentuais
conf_df <- function(real, estimado) {
  est_par  <- parear_otimo(real, estimado)
  g_levels <- as.character(sort(unique(c(real, estimado))))
  ct <- as.data.frame(table(
    Real     = factor(as.character(real),     levels = g_levels),
    Estimado = factor(as.character(est_par),  levels = g_levels)
  ))
  ct$pct      <- ct$Freq / sum(ct$Freq) * 100
  ct$Real_rev <- factor(ct$Real, levels = rev(g_levels))
  ct
}

# Helper: plot de uma matriz de confusão
plot_conf <- function(ct, titulo = NULL, subtitulo = NULL) {
  ggplot(ct, aes(Estimado, Real_rev)) +
    geom_tile(aes(fill = Freq), colour = "white", linewidth = 0.7) +
    geom_text(
      aes(label = ifelse(Freq > 0, paste0(Freq, "\n(", round(pct, 1), "%)"), "")),
      size = 2.7, lineheight = 1.05, colour = "grey10"
    ) +
    scale_fill_gradient(low = "#f0f7ff", high = "#2980b9", name = "n") +
    scale_x_discrete(position = "top") +
    labs(x = "Grupo estimado", y = "Grupo real",
         title = titulo, subtitle = subtitulo) +
    TEMA +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 9))
}

# Identifica o melhor conjunto NB (maior ARI) para uso nos plots PCA e perfis
tab_nb   <- tab[tab$Método == "Mistura NB", ]
label_best <- tab_nb$Conjunto[which.max(tab_nb$ARI)]
BEST_V   <- names(LABELS_CONJUNTOS)[LABELS_CONJUNTOS == label_best]
cat(sprintf("\nMelhor conjunto (NB · ARI): %s — %s\n\n", BEST_V, label_best))


# =============================================================================
# ── 9. VIS 1 · Matrizes de confusão — todos os conjuntos × métodos ────────────
# =============================================================================

cat("Gerando figuras...\n")

conf_plots <- list()
for (nome_v in names(CONJUNTOS)) {
  if (!nome_v %in% names(fits_nb) || !nome_v %in% names(km_results)) next

  ct_nb <- conf_df(grupo_real_int, fits_nb[[nome_v]]$classification)
  ct_km <- conf_df(grupo_real_int, km_results[[nome_v]]$cluster)
  met_nb <- calc_metricas(grupo_real_int, fits_nb[[nome_v]]$classification)
  met_km <- calc_metricas(grupo_real_int, km_results[[nome_v]]$cluster)

  conf_plots[[nome_v]] <- list(
    nb = plot_conf(
      ct_nb,
      titulo    = paste0(nome_v, " · Mistura NB (g=", fits_nb[[nome_v]]$g, ")"),
      subtitulo = sprintf("ARI=%.3f  ·  Acurácia=%.3f", met_nb$ari, met_nb$acuracia)
    ),
    km = plot_conf(
      ct_km,
      titulo    = paste0(nome_v, " · K-means (k=", km_results[[nome_v]]$k, ")"),
      subtitulo = sprintf("ARI=%.3f  ·  Acurácia=%.3f", met_km$ari, met_km$acuracia)
    )
  )
}

# Grade 5 linhas × 2 colunas: [V1_NB, V1_KM] / [V2_NB, V2_KM] / ...
plots_conf_list <- unlist(
  lapply(names(conf_plots), function(v) list(conf_plots[[v]]$nb, conf_plots[[v]]$km)),
  recursive = FALSE
)

fig_conf <- wrap_plots(plots_conf_list, ncol = 2) +
  plot_annotation(
    title    = "Matrizes de Confusão — Mistura NB vs. K-means",
    subtitle = sprintf(
      "Rótulos pareados via algoritmo húngaro  ·  n = %d beneficiários  ·  %d grupos reais",
      nrow(df), G_REAL
    ),
    caption  = paste0(
      "Diagonal = classificações corretas  ·  Off-diagonal = trocas de grupo\n",
      LEGENDA_GRUPOS
    ),
    theme    = TEMA
  )

ggsave(file.path(DIR_FIGS, "01_matrizes_confusao.png"),
       fig_conf, width = 11, height = 4.2 * length(conf_plots), dpi = 150, bg = "white")
cat("  01_matrizes_confusao.png\n")


# =============================================================================
# ── 10. VIS 2 · PCA — verdadeiro | EM | K-means lado a lado ──────────────────
# =============================================================================

vars_best <- CONJUNTOS[[BEST_V]]
Y_best    <- as.matrix(df[, vars_best])

pca_res <- prcomp(log1p(Y_best), scale. = TRUE)
var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
lab_pc  <- c(paste0("PC1 (", var_exp[1], "% var.)"),
             paste0("PC2 (", var_exp[2], "% var.)"))

pca_df <- data.frame(
  PC1        = pca_res$x[, 1],
  PC2        = pca_res$x[, 2],
  Verdadeiro = as.character(grupo_real_int),
  EM         = as.character(parear_otimo(grupo_real_int, fits_nb[[BEST_V]]$classification)),
  KMeans     = as.character(parear_otimo(grupo_real_int, km_results[[BEST_V]]$cluster))
)

make_pca <- function(col_grupo, titulo, subtitulo = NULL) {
  ggplot(pca_df, aes(.data[["PC1"]], .data[["PC2"]],
                     colour = .data[[col_grupo]])) +
    geom_point(alpha = 0.50, size = 1.3) +
    scale_colour_manual(values = COR, name = "Grupo") +
    labs(x = lab_pc[1], y = lab_pc[2], title = titulo, subtitle = subtitulo) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    TEMA + theme(legend.position = "right")
}

met_nb_best <- calc_metricas(grupo_real_int, fits_nb[[BEST_V]]$classification)
met_km_best <- calc_metricas(grupo_real_int, km_results[[BEST_V]]$cluster)

fig_pca <- (
  make_pca("Verdadeiro", "Grupo verdadeiro", "referência") |
  make_pca("EM",
           paste0("Mistura NB  (g=", fits_nb[[BEST_V]]$g, ")"),
           sprintf("ARI=%.3f  ·  Acurácia=%.3f", met_nb_best$ari, met_nb_best$acuracia)) |
  make_pca("KMeans",
           paste0("K-means  (k=", km_results[[BEST_V]]$k, ")"),
           sprintf("ARI=%.3f  ·  Acurácia=%.3f", met_km_best$ari, met_km_best$acuracia))
) +
  plot_annotation(
    title    = paste0("Projeção PCA — ", BEST_V, ": ", LABELS_CONJUNTOS[BEST_V]),
    subtitle = "Espaço log1p padronizado  ·  Rótulos pareados via algoritmo húngaro",
    theme    = TEMA
  )

ggsave(file.path(DIR_FIGS, "02_pca_classificacao.png"),
       fig_pca, width = 15, height = 5, dpi = 150, bg = "white")
cat("  02_pca_classificacao.png\n")


# =============================================================================
# ── 10b. VIS 2b · PCA por grupo em destaque (4 linhas × 3 painéis) ────────────
# =============================================================================
# Cada linha isola um grupo real: pontos do grupo em foco = cor viva;
# demais pontos = cinza claro · permite avaliar a recuperação grupo a grupo

# Cores vivas por grupo (mais saturadas que a paleta EDA)
# Ordem: 1=baixo_uso  2=ambulatorial  3=agudo  4=atipico
COR_VIVID <- c("1" = "#1552A0", "2" = "#1A7A3E", "3" = "#C07000", "4" = "#B52415")

# Rótulos legíveis para os títulos de linha
NOME_GRUPO_LABEL <- setNames(
  c("G1 · baixo uso",
    "G2 · ambulatorial coordenado",
    "G3 · agudo / hospitalar",
    "G4 · padrão atípico"),
  as.character(1:4)
)

make_pca_focus <- function(col_grupo, titulo, subtitulo = NULL, fg) {
  fg_chr  <- as.character(fg)
  df_back <- pca_df[pca_df[[col_grupo]] != fg_chr, ]
  df_fore <- pca_df[pca_df[[col_grupo]] == fg_chr, ]
  ggplot() +
    geom_point(data = df_back,
               aes(PC1, PC2),
               colour = "#C8C8C8", alpha = 0.20, size = 0.9) +
    geom_point(data = df_fore,
               aes(PC1, PC2),
               colour = COR_VIVID[fg_chr], alpha = 0.80, size = 1.5) +
    labs(x = lab_pc[1], y = lab_pc[2], title = titulo, subtitle = subtitulo) +
    TEMA + theme(legend.position = "none")
}

focus_rows <- lapply(seq_len(G_REAL), function(fg) {
  lbl <- NOME_GRUPO_LABEL[as.character(fg)]
  (
    make_pca_focus("Verdadeiro",
                   paste0(lbl, " — verdadeiro"),
                   "referência", fg) |
    make_pca_focus("EM",
                   paste0(lbl, " — Mistura NB (g=", g_best, ")"),
                   sprintf("ARI=%.3f  ·  Acurácia=%.3f",
                           met_nb_best$ari, met_nb_best$acuracia), fg) |
    make_pca_focus("KMeans",
                   paste0(lbl, " — K-means (k=", km_results[[BEST_V]]$k, ")"),
                   sprintf("ARI=%.3f  ·  Acurácia=%.3f",
                           met_km_best$ari, met_km_best$acuracia), fg)
  )
})

fig_focus <- wrap_plots(focus_rows, ncol = 1) +
  plot_annotation(
    title    = paste0("Classificação por grupo — projeção PCA · ", BEST_V,
                      ": ", LABELS_CONJUNTOS[BEST_V]),
    subtitle = paste0(
      "Cada linha destaca um grupo: cor viva = grupo em foco  ·  cinza = demais\n",
      LEGENDA_GRUPOS
    ),
    theme    = TEMA
  )

ggsave(file.path(DIR_FIGS, "08_pca_por_grupo.png"),
       fig_focus, width = 15, height = 20, dpi = 150, bg = "white")
cat("  08_pca_por_grupo.png\n")


# =============================================================================
# ── 11. VIS 3 · Curvas de seleção: BIC (EM) e Silhouette (K-means) ────────────
# =============================================================================

# BIC —————————————————————————————————————————————————————————————————————————
bic_long <- do.call(rbind, lapply(names(bic_curves), function(v) {
  bic_v <- bic_curves[[v]]
  data.frame(
    ConjID      = v,
    Conjunto    = LABELS_CONJUNTOS[v],
    g           = seq_len(G_MAX),
    BIC         = bic_v,
    selecionado = seq_len(G_MAX) == which.min(bic_v)
  )
}))
bic_long$Conjunto <- factor(bic_long$Conjunto, levels = as.character(LABELS_CONJUNTOS))

p_bic <- ggplot(bic_long, aes(g, BIC, colour = Conjunto, group = Conjunto)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  geom_point(data = bic_long[bic_long$selecionado, ],
             shape = 21, size = 5, stroke = 1.6,
             colour = "grey20", fill = NA) +
  scale_x_continuous(breaks = seq_len(G_MAX)) +
  scale_colour_brewer(palette = "Set1", name = "") +
  labs(x = "Número de componentes (g)", y = "BIC",
       title = "Curvas BIC — Mistura NB",
       subtitle = "Círculo destacado = g* selecionado (mínimo BIC)") +
  TEMA

# Silhouette ——————————————————————————————————————————————————————————————————
sil_long <- do.call(rbind, lapply(names(sil_curves), function(v) {
  sil_v <- sil_curves[[v]]
  k_seq <- 2:G_MAX
  data.frame(
    ConjID      = v,
    Conjunto    = LABELS_CONJUNTOS[v],
    k           = k_seq,
    Sil         = sil_v,
    selecionado = k_seq == k_seq[which.max(sil_v)]
  )
}))
sil_long$Conjunto <- factor(sil_long$Conjunto, levels = as.character(LABELS_CONJUNTOS))

p_sil <- ggplot(sil_long, aes(k, Sil, colour = Conjunto, group = Conjunto)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  geom_point(data = sil_long[sil_long$selecionado, ],
             shape = 21, size = 5, stroke = 1.6,
             colour = "grey20", fill = NA) +
  scale_x_continuous(breaks = 2:G_MAX) +
  scale_colour_brewer(palette = "Set1", name = "") +
  labs(x = "Número de clusters (k)", y = "Silhouette médio",
       title = "Curvas Silhouette — K-means",
       subtitle = "Círculo destacado = k* selecionado (máximo silhouette)") +
  TEMA

fig_sel <- (p_bic | p_sil) +
  plot_annotation(
    title   = "Critérios de seleção do número de grupos",
    caption = "BIC penaliza complexidade via d(g)=2gp+(g−1)  ·  Silhouette mede coesão/separação intra-cluster",
    theme   = TEMA
  )

ggsave(file.path(DIR_FIGS, "03_curvas_selecao.png"),
       fig_sel, width = 13, height = 5, dpi = 150, bg = "white")
cat("  03_curvas_selecao.png\n")


# =============================================================================
# ── 12. VIS 4 · Perfis dos componentes (μ̂_hj e θ̂_hj) — melhor modelo EM ─────
# =============================================================================

mu_mat    <- fits_nb[[BEST_V]]$mu      # g × p
theta_mat <- fits_nb[[BEST_V]]$theta   # g × p
g_best    <- fits_nb[[BEST_V]]$g
vars_best_lbl <- CONJUNTOS[[BEST_V]]

comp_labels <- paste0("G", seq_len(g_best))
rownames(mu_mat)    <- comp_labels
rownames(theta_mat) <- comp_labels
colnames(mu_mat)    <- vars_best_lbl
colnames(theta_mat) <- vars_best_lbl

# Reshape μ → long
mu_long <- data.frame(
  Componente = rep(comp_labels, each = length(vars_best_lbl)),
  Variavel   = rep(vars_best_lbl, times = g_best),
  mu         = as.vector(t(mu_mat))
)
mu_long$Variavel   <- factor(mu_long$Variavel,   levels = vars_best_lbl)
mu_long$Componente <- factor(mu_long$Componente, levels = rev(comp_labels))

# z-score por variável (permite comparar magnitude relativa entre variáveis)
mu_long$mu_z <- ave(mu_long$mu, mu_long$Variavel,
                    FUN = function(x) (x - mean(x)) / (sd(x) + 1e-9))

p_mu_abs <- ggplot(mu_long, aes(Variavel, Componente)) +
  geom_tile(aes(fill = mu), colour = "white", linewidth = 0.7) +
  geom_text(aes(label = round(mu, 1)), size = 2.9, colour = "grey10") +
  scale_fill_gradient(low = "#f0f7ff", high = "#1a5276", name = "μ̂") +
  labs(title = "Médias estimadas μ̂_hj",
       subtitle = "Valores absolutos por componente × variável",
       x = NULL, y = NULL) +
  TEMA +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

p_mu_z <- ggplot(mu_long, aes(Variavel, Componente)) +
  geom_tile(aes(fill = mu_z), colour = "white", linewidth = 0.7) +
  geom_text(aes(label = round(mu_z, 2)), size = 2.9, colour = "grey10") +
  scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#2980b9",
                       midpoint = 0, name = "z") +
  labs(title = "Médias padronizadas (z-score por variável)",
       subtitle = "Destaca quais variáveis mais distinguem cada componente",
       x = NULL, y = NULL) +
  TEMA +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# Reshape θ → long
theta_long <- data.frame(
  Componente = rep(comp_labels, each = length(vars_best_lbl)),
  Variavel   = rep(vars_best_lbl, times = g_best),
  theta      = as.vector(t(theta_mat))
)
theta_long$Variavel   <- factor(theta_long$Variavel,   levels = vars_best_lbl)
theta_long$Componente <- factor(theta_long$Componente, levels = rev(comp_labels))
# Variância = μ + μ²/θ; θ pequeno = alta sobredispersão
theta_long$sobredisp <- ifelse(theta_long$theta < 1, "alta",
                        ifelse(theta_long$theta < 5, "moderada", "baixa"))

p_theta <- ggplot(theta_long, aes(Variavel, Componente)) +
  geom_tile(aes(fill = log10(theta)), colour = "white", linewidth = 0.7) +
  geom_text(aes(label = round(theta, 1)), size = 2.9, colour = "grey10") +
  scale_fill_gradient(low = "#fdf2e9", high = "#b7410e",
                      name = "log₁₀(θ̂)",
                      labels = function(x) round(10^x, 1)) +
  labs(title = "Dispersão estimada θ̂_hj",
       subtitle = "Cor na escala log  ·  θ→0 = alta sobredispersão  ·  θ→∞ ≈ Poisson",
       x = NULL, y = NULL) +
  TEMA +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

fig_prof <- (p_mu_abs | p_mu_z | p_theta) +
  plot_annotation(
    title    = paste0("Perfis dos componentes — ", BEST_V, ": ", LABELS_CONJUNTOS[BEST_V]),
    subtitle = paste0("Mistura NB  ·  g=", g_best, " componentes"),
    caption  = paste0("Var[Y] = μ + μ²/θ  ·  Estimativas obtidas pelo algoritmo EM  ·  n=",
                      nrow(df), " beneficiários"),
    theme    = TEMA
  )

ggsave(file.path(DIR_FIGS, "04_perfis_componentes.png"),
       fig_prof, width = 15, height = 3 + g_best * 0.6, dpi = 150, bg = "white")
cat("  04_perfis_componentes.png\n")


# =============================================================================
# ── 13. VIS 5 · Incerteza de classificação — distribuição de max(τ̂_ih) ────────
# =============================================================================

tau_df <- do.call(rbind, lapply(names(fits_nb), function(v) {
  tau_v   <- fits_nb[[v]]$tau
  g_v     <- fits_nb[[v]]$g
  max_tau <- apply(tau_v, 1, max)
  data.frame(
    ConjID   = v,
    Conjunto = factor(LABELS_CONJUNTOS[v], levels = as.character(LABELS_CONJUNTOS)),
    max_tau  = max_tau,
    limiar   = 1 / g_v,        # chance uniforme (incerteza máxima)
    g        = g_v
  )
}))

# Mediana por conjunto
med_tau <- aggregate(max_tau ~ Conjunto, tau_df, median)
# Proporção com alta certeza (max_tau > 0.9)
cert_tau <- aggregate(max_tau ~ Conjunto, tau_df,
                      FUN = function(x) round(100 * mean(x > 0.9), 1))
names(cert_tau)[2] <- "pct_cert"
tau_info <- merge(med_tau, cert_tau, by = "Conjunto")
tau_info$label <- sprintf("med=%.2f · %.0f%% c/ τ>0.9", tau_info$max_tau, tau_info$pct_cert)

# Limiares únicos por conjunto (para geom_vline)
limiares <- unique(tau_df[, c("Conjunto", "limiar")])

p_unc <- ggplot(tau_df, aes(max_tau)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "#2980b9", colour = "white", alpha = 0.70) +
  geom_density(colour = "#1a5276", linewidth = 0.85) +
  geom_vline(data = med_tau,
             aes(xintercept = max_tau),
             colour = "#c0392b", linetype = "dashed", linewidth = 0.75) +
  geom_vline(data = limiares,
             aes(xintercept = limiar),
             colour = "#e6a817", linetype = "dotted", linewidth = 0.75) +
  geom_text(data = tau_info,
            aes(x = 0.02, y = Inf, label = label),
            hjust = 0, vjust = 1.4, size = 2.8, colour = "grey30") +
  facet_wrap(~ Conjunto, ncol = 3, scales = "free_y") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x       = "max τ̂_ih  (certeza de pertencimento ao componente mais provável)",
    y       = "Densidade",
    title   = "Incerteza de classificação — Mistura NB",
    subtitle = "Vermelho tracejado = mediana  ·  Amarelo pontilhado = chance aleatória (1/g)",
    caption  = "Valores próximos de 1 = beneficiário bem separado  ·  Valores próximos de 1/g = sobreposição entre grupos"
  ) +
  TEMA

ggsave(file.path(DIR_FIGS, "05_incerteza_classificacao.png"),
       p_unc, width = 13, height = 7, dpi = 150, bg = "white")
cat("  05_incerteza_classificacao.png\n")


# =============================================================================
# ── 14. VIS 6 · Evolução das métricas V1 → V5 ────────────────────────────────
# =============================================================================

tab_evo <- tab
tab_evo$Conjunto <- factor(tab_evo$Conjunto,
                            levels = as.character(LABELS_CONJUNTOS))

# Long: ARI, Acurácia, NMI, Pureza
tab_long <- pivot_longer(
  tab_evo,
  cols      = c("ARI", "Acuracia", "NMI", "Pureza"),
  names_to  = "Metrica",
  values_to = "Valor"
)
tab_long$Metrica <- factor(tab_long$Metrica,
                            levels = c("ARI", "Acuracia", "NMI", "Pureza"),
                            labels = c("ARI (concordância de pares)",
                                       "Acurácia (após pareamento húngaro)",
                                       "NMI (informação mútua norm.)",
                                       "Pureza"))

p_evo <- ggplot(tab_long, aes(Conjunto, Valor,
                               colour = Método, group = Método)) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(shape = Método), size = 3.2) +
  geom_text(
    data = tab_long[tab_long$Conjunto == levels(tab_long$Conjunto)[nlevels(tab_long$Conjunto)], ],
    aes(label = sprintf("%.3f", Valor)),
    hjust = -0.15, size = 2.7, show.legend = FALSE
  ) +
  facet_wrap(~ Metrica, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c("Mistura NB" = "#4a7c59", "K-means" = "#c0392b"),
                      name = "") +
  scale_shape_manual(values  = c("Mistura NB" = 16,        "K-means" = 17),
                     name = "") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    x       = NULL,
    y       = "Valor",
    title   = "Evolução das métricas conforme variáveis são adicionadas (V1 → V3)",
    subtitle = "V1 = baseline 3 variáveis  ·  cada Vi adiciona uma variável ao conjunto anterior",
    caption  = "V4 e V5 desativados: não melhoram o modelo (confirmado empiricamente)"
  ) +
  TEMA +
  theme(axis.text.x = element_text(angle = 28, hjust = 1))

ggsave(file.path(DIR_FIGS, "06_evolucao_metricas.png"),
       p_evo, width = 13, height = 7, dpi = 150, bg = "white")
cat("  06_evolucao_metricas.png\n")


# =============================================================================
# ── 15. VIS 7 · Painel resumo — melhor modelo ────────────────────────────────
# =============================================================================
# Combina num único painel de alta síntese: PCA (3 col) + confusão (2 col) +
# incerteza posterior do melhor conjunto

ct_nb_best <- conf_df(grupo_real_int, fits_nb[[BEST_V]]$classification)
ct_km_best <- conf_df(grupo_real_int, km_results[[BEST_V]]$cluster)

p_conf_nb_best <- plot_conf(
  ct_nb_best,
  titulo    = paste0("Mistura NB (g=", g_best, ")"),
  subtitulo = sprintf("ARI=%.3f · Acurácia=%.3f", met_nb_best$ari, met_nb_best$acuracia)
)
p_conf_km_best <- plot_conf(
  ct_km_best,
  titulo    = paste0("K-means (k=", km_results[[BEST_V]]$k, ")"),
  subtitulo = sprintf("ARI=%.3f · Acurácia=%.3f", met_km_best$ari, met_km_best$acuracia)
)

tau_best_df <- data.frame(
  max_tau = apply(fits_nb[[BEST_V]]$tau, 1, max),
  Grupo   = as.character(fits_nb[[BEST_V]]$classification)
)
p_tau_best <- ggplot(tau_best_df, aes(max_tau, fill = Grupo, colour = Grupo)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values   = COR, name = "Componente") +
  scale_colour_manual(values = COR, name = "Componente") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  geom_vline(xintercept = 1 / g_best,
             colour = "#e6a817", linetype = "dotted", linewidth = 0.8) +
  labs(x = "max τ̂_ih", y = "Densidade",
       title = "Certeza por componente",
       subtitle = "Amarelo pontilhado = chance aleatória (1/g)") +
  TEMA + theme(legend.position = "right")

fig_resumo <- (
  (make_pca("Verdadeiro", "Verdadeiro") |
   make_pca("EM",         paste0("EM (g=", g_best, ")")) |
   make_pca("KMeans",     paste0("K-means (k=", km_results[[BEST_V]]$k, ")"))) /
  (p_conf_nb_best | p_conf_km_best | p_tau_best)
) +
  plot_annotation(
    title    = paste0("Painel resumo — melhor conjunto: ", BEST_V,
                      " (", LABELS_CONJUNTOS[BEST_V], ")"),
    subtitle = paste0("Linha superior: projeção PCA  ·  ",
                      "Linha inferior: matrizes de confusão e incerteza posterior"),
    caption  = LEGENDA_GRUPOS,
    theme    = TEMA
  )

ggsave(file.path(DIR_FIGS, "07_painel_resumo.png"),
       fig_resumo, width = 15, height = 11, dpi = 150, bg = "white")
cat("  07_painel_resumo.png\n")


# ── Resumo final ──────────────────────────────────────────────────────────────

cat(sprintf("\n%s\n", strrep("═", 60)))
cat(sprintf("Figuras salvas em '%s/':\n", DIR_FIGS))
cat("  01_matrizes_confusao.png    — grade 5×2 (conjuntos × métodos)\n")
cat("  02_pca_classificacao.png    — PCA: verdadeiro | EM | K-means\n")
cat("  03_curvas_selecao.png       — BIC (EM) e Silhouette (K-means)\n")
cat("  04_perfis_componentes.png   — μ̂ absoluto | μ̂ z-score | θ̂\n")
cat("  05_incerteza_classificacao.png — distribuição de max(τ̂_ih)\n")
cat("  06_evolucao_metricas.png    — ARI/Acurácia/NMI/Pureza: V1→V5\n")
cat("  07_painel_resumo.png        — painel de síntese (melhor conjunto)\n")
cat(strrep("═", 60), "\n")
cat("Próximo passo: usar rds/comparacao_variaveis.rds no analises.Rmd\n")
