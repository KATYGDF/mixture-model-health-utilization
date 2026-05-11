# =============================================================================
# eda_v3.R
# Análise Exploratória de Dados — V3: 5 variáveis clínicas
#   consultas · ps · exames · internacoes · terapias
#
# Saída: figuras PNG em docs/assets/img/eda/
#        tabela CSV  em docs/assets/img/eda/tabela_descritiva.csv
#
# USO: Source (Ctrl+Shift+S) no RStudio a partir de src/02_eda/
# =============================================================================

# ── 1. Parâmetros ─────────────────────────────────────────────────────────────

CAMINHO_CSV <- "../../dados/dados_simulados.csv"
DIR_FIGS    <- "../../docs/assets/img/eda"
SEED        <- 42L
VARS_V3     <- c("consultas", "ps", "exames", "internacoes", "terapias")

GRUPOS_LABELS <- c(
  agudo_hospitalar        = "Agudo/hospitalar",
  ambulatorial_coordenado = "Ambulatorial",
  atipico                 = "Atípico",
  baixo_uso               = "Baixo uso"
)

COR_G1 <- "#7BAFD4"   # baixo_uso
COR_G2 <- "#6EBF8B"   # ambulatorial_coordenado
COR_G3 <- "#E8B96B"   # agudo_hospitalar
COR_G4 <- "#D98585"   # atipico

COR_GRUPO <- c(
  agudo_hospitalar        = COR_G3,
  ambulatorial_coordenado = COR_G2,
  atipico                 = COR_G4,
  baixo_uso               = COR_G1
)

# ── 2. Pacotes ────────────────────────────────────────────────────────────────

pkgs  <- c("ggplot2", "patchwork", "scales", "dplyr", "tidyr", "corrplot")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Instalando: ", paste(novos, collapse = ", "))
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
  library(dplyr);   library(tidyr);    library(corrplot)
})

if (!dir.exists(DIR_FIGS)) dir.create(DIR_FIGS, recursive = TRUE)

# ── 3. Dados ──────────────────────────────────────────────────────────────────

df_raw <- read.csv(CAMINHO_CSV, stringsAsFactors = FALSE)
candidatos <- c("grupo_latente", "grupo", "grupo_real", "grupo_idx")
col_grupo  <- candidatos[candidatos %in% names(df_raw)][1]
df <- df_raw
names(df)[names(df) == col_grupo] <- "grupo"
df$grupo <- factor(df$grupo)

n  <- nrow(df)
ng <- nlevels(df$grupo)
cat(sprintf("Dados: %d obs · %d grupos · V3 = %d variáveis\n\n", n, ng, length(VARS_V3)))

# ── 4. Tema e utilitários ─────────────────────────────────────────────────────

TEMA <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 9),
    strip.background  = element_rect(fill = "#f0f0f0", colour = NA),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(colour = "grey40", size = 9),
    plot.caption      = element_text(colour = "grey55", size = 8, hjust = 0),
    legend.position   = "bottom"
  )

salvar <- function(p, nome, w = 10, h = 5.5) {
  path <- file.path(DIR_FIGS, nome)
  ggsave(path, p, width = w, height = h, dpi = 150, bg = "white")
  cat(sprintf("  Salvo: %s\n", nome))
}

# ── 5. Tabela descritiva ──────────────────────────────────────────────────────

cat("── Tabela descritiva ────────────────────────────────────────\n")

tab_desc <- do.call(rbind, lapply(VARS_V3, function(v) {
  x <- df[[v]]
  data.frame(
    variavel  = v,
    n         = length(x),
    media     = round(mean(x),   2),
    dp        = round(sd(x),     2),
    mediana   = round(median(x), 2),
    p25       = round(quantile(x, .25), 1),
    p75       = round(quantile(x, .75), 1),
    p99       = round(quantile(x, .99), 1),
    zeros_pct = round(100 * mean(x == 0), 1),
    razao_v_m = round(var(x) / mean(x), 1)
  )
}))

print(tab_desc, row.names = FALSE)
write.csv(tab_desc, file.path(DIR_FIGS, "tabela_descritiva.csv"), row.names = FALSE)


# ── 6. FIG 01 · Distribuição marginal das 5 variáveis ────────────────────────

cat("\n── Fig 01: distribuição marginal ────────────────────────────\n")

plots_marginal <- lapply(VARS_V3, function(v) {
  x    <- df[[v]]
  x_p99 <- quantile(x, 0.97)
  x_tr  <- x[x <= x_p99]

  ggplot(data.frame(x = x_tr), aes(x)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill = "#2980b9", colour = "white", alpha = 0.75) +
    geom_density(colour = "#1a5276", linewidth = 0.8) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x        = v,
      y        = "Densidade",
      title    = v,
      subtitle = sprintf("média=%.2f · DP=%.2f · zeros=%.1f%% · Var/Média=%.1f",
                         mean(x), sd(x), 100 * mean(x == 0), var(x) / mean(x))
    ) +
    TEMA + theme(legend.position = "none",
                 plot.title    = element_text(size = 10, face = "bold"),
                 plot.subtitle = element_text(size = 8))
})

fig01 <- wrap_plots(plots_marginal, ncol = 3) +
  plot_annotation(
    title    = "Distribuição marginal das variáveis clínicas (V3)",
    subtitle = sprintf("n = %d beneficiários  ·  eixo x truncado no percentil 97", n),
    caption  = "Histograma + densidade kernel  ·  dados simulados calibrados ao Painel ANS",
    theme    = TEMA
  )

salvar(fig01, "fig01_marginal.png", w = 14, h = 8)


# ── 6b. FIG 08 · Ridgeline por grupo (log1p) ─────────────────────────────────

cat("── Fig 08: ridgeline por grupo ────────────────────────────────\n")

if (!requireNamespace("ggridges", quietly = TRUE))
  install.packages("ggridges", repos = "https://cloud.r-project.org", quiet = TRUE)
library(ggridges)

VARS_LABELS <- setNames(
  paste0("log(", VARS_V3, " + 1)"),
  VARS_V3
)

df_ridge <- do.call(rbind, lapply(VARS_V3, function(v) {
  data.frame(
    variavel = VARS_LABELS[v],
    valor    = log1p(df[[v]]),
    grupo    = df$grupo
  )
}))
df_ridge$grupo <- factor(df_ridge$grupo,
  levels = c("baixo_uso", "ambulatorial_coordenado", "agudo_hospitalar", "atipico"))
df_ridge$variavel <- factor(df_ridge$variavel, levels = VARS_LABELS)

fig08 <- ggplot(df_ridge, aes(x = valor, y = grupo, fill = grupo)) +
  ggridges::geom_density_ridges(
    alpha          = 0.75,
    scale          = 1.8,
    colour         = "white",
    linewidth      = 0.35,
    rel_min_height = 0.01
  ) +
  scale_fill_manual(values = COR_GRUPO, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  facet_wrap(~ variavel, scales = "free_x", ncol = 2) +
  labs(
    x        = NULL,
    y        = NULL,
    title    = "Distribuição das variáveis por grupo latente",
    subtitle = "Ridgeline: separação entre grupos evidencia estrutura de mistura",
    caption  = sprintf("n = %d  ·  escala log(y + 1)  ·  grupos verdadeiros (gabarito simulado)", n)
  ) +
  TEMA +
  theme(
    axis.text.y   = element_text(size = 8.5),
    panel.spacing = unit(1, "lines")
  )

salvar(fig08, "fig08_densidade_grupo.png", w = 10, h = 12)


# ── 7. FIG 02 · Razão Var/Média (sobredispersão) ──────────────────────────────

cat("── Fig 02: sobredispersão ────────────────────────────────────\n")

df_sd <- tab_desc[, c("variavel", "razao_v_m")]
df_sd$cor <- ifelse(df_sd$razao_v_m > 15, "alta",
             ifelse(df_sd$razao_v_m > 8,  "moderada-alta", "moderada"))
df_sd$variavel <- factor(df_sd$variavel, levels = rev(VARS_V3))

fig02 <- ggplot(df_sd, aes(razao_v_m, variavel, fill = cor)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.1f×", razao_v_m)), hjust = -0.15, size = 3.5) +
  geom_vline(xintercept = 1, colour = "grey30", linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 1.6, y = 0.55, label = "Poisson\n(razão = 1)",
           hjust = 0, size = 3.1, colour = "grey35") +
  scale_fill_manual(
    values  = c(alta = "#c0392b", `moderada-alta` = "#e6a817", moderada = "#4a7c59"),
    name    = "Sobredispersão",
    guide   = guide_legend(reverse = TRUE)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x        = "Razão Variância / Média",
    y        = NULL,
    title    = "Sobredispersão das variáveis clínicas (V3)",
    subtitle = "Todas as variáveis violam a equidispersão da Poisson (razão = 1)",
    caption  = sprintf("n = %d  ·  grupos simulados = %d", n, ng)
  ) +
  TEMA + theme(legend.position = "right")

salvar(fig02, "fig02_sobredispersao.png", w = 8, h = 4.5)


# ── 8. FIG 03 · Zeros por variável e grupo ────────────────────────────────────

cat("── Fig 03: proporção de zeros ────────────────────────────────\n")

df_zeros <- do.call(rbind, lapply(VARS_V3, function(v) {
  do.call(rbind, lapply(levels(df$grupo), function(g) {
    sub <- df[df$grupo == g, v]
    data.frame(
      variavel = v,
      grupo    = g,
      pct_zero = 100 * mean(sub == 0)
    )
  }))
}))
df_zeros$variavel <- factor(df_zeros$variavel, levels = VARS_V3)
df_zeros$grupo    <- factor(df_zeros$grupo)

fig03 <- ggplot(df_zeros, aes(variavel, pct_zero, fill = grupo)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.4) +
  scale_fill_manual(values = COR_GRUPO, name = "Grupo") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    x        = NULL,
    y        = "% zeros",
    title    = "Proporção de zeros por variável e grupo latente",
    subtitle = "Internações e terapias têm excesso de zeros estrutural; padrão varia entre grupos",
    caption  = sprintf("n = %d  ·  grupos verdadeiros (gabarito simulado)", n)
  ) +
  TEMA

salvar(fig03, "fig03_zeros_grupo.png", w = 10, h = 5.5)


# ── 9. FIG 04 · Boxplots por grupo ────────────────────────────────────────────

cat("── Fig 04: boxplots por grupo ────────────────────────────────\n")

df_long <- df %>%
  dplyr::select(all_of(c("grupo", VARS_V3))) %>%
  pivot_longer(cols = all_of(VARS_V3), names_to = "variavel", values_to = "valor")
df_long$variavel <- factor(df_long$variavel, levels = VARS_V3)

# Truncar ao P97 por variável para visualização
df_long <- df_long %>%
  group_by(variavel) %>%
  mutate(p97 = quantile(valor, 0.97)) %>%
  filter(valor <= p97) %>%
  ungroup()

fig04 <- ggplot(df_long, aes(grupo, valor, fill = grupo)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.5) +
  scale_fill_manual(values = COR_GRUPO, guide = "none") +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
  facet_wrap(~ variavel, scales = "free_y", ncol = 3) +
  labs(
    x        = NULL,
    y        = "Contagem",
    title    = "Distribuição das variáveis clínicas por grupo latente (V3)",
    subtitle = "Eixo y livre por variável  ·  Outliers acima do P97 ocultados",
    caption  = sprintf("n = %d  ·  grupos verdadeiros (gabarito simulado)", n)
  ) +
  TEMA + theme(axis.text.x = element_text(size = 7.5))

salvar(fig04, "fig04_boxplots_grupo.png", w = 14, h = 9)


# ── 10. FIG 05 · Médias por grupo (perfil de utilização) ─────────────────────

cat("── Fig 05: perfil de utilização por grupo ────────────────────\n")

df_medias <- df %>%
  group_by(grupo) %>%
  summarise(across(all_of(VARS_V3), mean), .groups = "drop") %>%
  pivot_longer(cols = all_of(VARS_V3), names_to = "variavel", values_to = "media")

# Normalizar pela média global para comparação radar-like
df_global <- df %>% summarise(across(all_of(VARS_V3), mean))
df_medias <- df_medias %>%
  group_by(variavel) %>%
  mutate(media_rel = media / df_global[[unique(variavel)]]) %>%
  ungroup()
df_medias$variavel <- factor(df_medias$variavel, levels = VARS_V3)

fig05 <- ggplot(df_medias, aes(variavel, media_rel, colour = grupo, group = grupo)) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.7, linetype = "dashed") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_colour_manual(values = COR_GRUPO, name = "Grupo") +
  scale_y_continuous(labels = function(x) sprintf("%.1f×", x)) +
  labs(
    x        = NULL,
    y        = "Média relativa à média global",
    title    = "Perfil de utilização por grupo latente (médias relativas)",
    subtitle = "Linha tracejada = média global (1,0×)  ·  Valores acima = uso acima da média",
    caption  = sprintf("n = %d  ·  grupos verdadeiros (gabarito simulado)", n)
  ) +
  TEMA

salvar(fig05, "fig05_perfil_grupo.png", w = 10, h = 5.5)


# ── 11. FIG 06 · Correlação de Spearman ───────────────────────────────────────

cat("── Fig 06: correlação de Spearman ───────────────────────────\n")

mat_cor <- cor(df[, VARS_V3], method = "spearman")
colnames(mat_cor) <- rownames(mat_cor) <- VARS_V3

# Salvar como PNG via corrplot
png(file.path(DIR_FIGS, "fig06_correlacao.png"),
    width = 900, height = 800, res = 150, bg = "white")
corrplot::corrplot(
  mat_cor,
  method   = "color",
  type     = "upper",
  addCoef.col = "black",
  number.cex  = 0.85,
  tl.col   = "black",
  tl.srt   = 45,
  tl.cex   = 0.9,
  cl.cex   = 0.8,
  col      = colorRampPalette(c("#c0392b", "white", "#2980b9"))(200),
  title    = "Correlação de Spearman — variáveis clínicas V3",
  mar      = c(0, 0, 2, 0)
)
dev.off()
cat(sprintf("  Salvo: fig06_correlacao.png\n"))


# ── 12. FIG 07 · Distribuição log1p por grupo (violin + jitter) ──────────────

cat("── Fig 07: violin log1p por grupo ────────────────────────────\n")

df_log <- df %>%
  dplyr::select(all_of(c("grupo", VARS_V3))) %>%
  mutate(across(all_of(VARS_V3), log1p)) %>%
  pivot_longer(cols = all_of(VARS_V3), names_to = "variavel", values_to = "log1p_valor")
df_log$variavel <- factor(df_log$variavel, levels = VARS_V3)

fig07 <- ggplot(df_log, aes(grupo, log1p_valor, fill = grupo, colour = grupo)) +
  geom_violin(alpha = 0.35, linewidth = 0.4, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
               colour = "grey30", linewidth = 0.5) +
  scale_fill_manual(values   = COR_GRUPO, guide = "none") +
  scale_colour_manual(values = COR_GRUPO, guide = "none") +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
  facet_wrap(~ variavel, scales = "free_y", ncol = 3) +
  labs(
    x        = NULL,
    y        = "log(1 + contagem)",
    title    = "Distribuição de log(1 + Y) por grupo e variável",
    subtitle = "Violino = forma da distribuição  ·  Caixa interna = Q1–Q3  ·  Linha = mediana",
    caption  = sprintf("n = %d  ·  escala log1p para visualização", n)
  ) +
  TEMA + theme(axis.text.x = element_text(size = 7.5))

salvar(fig07, "fig07_violin_log1p.png", w = 14, h = 9)


# ── 13. Resumo no console ─────────────────────────────────────────────────────

cat(sprintf("\n%s\n", strrep("═", 60)))
cat("RESUMO EDA — V3\n")
cat(strrep("═", 60), "\n\n")
cat(sprintf("n = %d  ·  grupos = %d\n\n", n, ng))

cat("Tamanho por grupo:\n")
print(table(df$grupo))

cat("\nMédias por grupo:\n")
for (v in VARS_V3) {
  cat(sprintf("  %-14s", v))
  m <- round(tapply(df[[v]], df$grupo, mean), 2)
  cat(paste(sprintf("%s=%.2f", names(m), m), collapse = "  "), "\n")
}

cat("\nFiguras salvas em:", DIR_FIGS, "\n")
cat("  fig01_marginal.png\n  fig02_sobredispersao.png\n  fig03_zeros_grupo.png\n")
cat("  fig04_boxplots_grupo.png\n  fig05_perfil_grupo.png\n")
cat("  fig06_correlacao.png\n  fig07_violin_log1p.png\n")
cat("  fig08_densidade_grupo.png\n")
cat("  tabela_descritiva.csv\n")
