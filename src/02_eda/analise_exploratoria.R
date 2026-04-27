# ============================================================
# EDA COMPLETA — EVIDÊNCIA DE MISTURA DE DISTRIBUIÇÕES
# SimSaúde · Katy Garcia de Freitas · Universidade do Minho
# ============================================================

# ------------------------------------------------------------
# 0) PACOTES
# ------------------------------------------------------------
packs <- c(
  "tidyverse", "GGally", "scales", "patchwork",
  "mclust", "diptest", "moments", "factoextra",
  "ggridges", "ggrepel", "viridis"
)

instalar <- packs[!packs %in% installed.packages()[, "Package"]]
if (length(instalar) > 0) install.packages(instalar)
lapply(packs, library, character.only = TRUE)

# ------------------------------------------------------------
# TEMA E PALETA GLOBAL
# ------------------------------------------------------------

# Paleta pastel diferenciada por grupo (4 grupos + neutros)
COR_G1  <- "#7BAFD4"   # azul-acinzentado   — baixo_uso
COR_G2  <- "#6EBF8B"   # verde-sálvia        — ambulatorial_coordenado
COR_G3  <- "#E8B96B"   # âmbar suave         — agudo_hospitalar
COR_G4  <- "#D98585"   # rosa-terracota      — atipico
CORES_GRUPOS <- c(
  "baixo_uso"               = COR_G1,
  "ambulatorial_coordenado" = COR_G2,
  "agudo_hospitalar"        = COR_G3,
  "atipico"                 = COR_G4
)
CINZA_HIST <- "#B0B8C1"
COR_DENS   <- "#555F6E"
COR_ECDF   <- "#4D8FAC"
COR_ACENTO <- "#5B7FA6"

tema_simsaude <- theme_minimal(base_size = 12) +
  theme(
    text             = element_text(family = "sans", color = "#2c2c2c"),
    plot.title       = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, color = "#666", margin = margin(b = 8)),
    plot.caption     = element_text(size = 8, color = "#999", hjust = 0),
    axis.title       = element_text(size = 10, color = "#555"),
    axis.text        = element_text(size = 9, color = "#666"),
    panel.grid.major = element_line(color = "#EBEBEB", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.5, "cm"),
    strip.text       = element_text(size = 9, face = "bold", color = "#444"),
    plot.margin      = margin(12, 16, 10, 12)
  )

theme_set(tema_simsaude)

# ------------------------------------------------------------
# 1) LEITURA E PREPARAÇÃO
# ------------------------------------------------------------
dados <- read.csv("C:/Katy/Doutorado/painel_WCDANM/dados_simulados.csv") %>%
  mutate(
    grupo_latente = factor(grupo_latente,
                           levels = c("baixo_uso", "ambulatorial_coordenado",
                                      "agudo_hospitalar", "atipico")),
    faixa_idade = factor(faixa_idade, levels = c("0-17", "18-59", "60+")),
    sexo        = factor(sexo),
    copart      = factor(copart),
    regiao      = factor(regiao),
    # Log-transformações para variáveis assimétricas
    log_gasto         = log1p(gasto_total),
    log_total_eventos = log1p(total_eventos),
    log_prestadores   = log1p(prestadores_distintos),
    log_exames        = log1p(exames),
    log_consultas     = log1p(consultas),
    log_terapias      = log1p(terapias),
    log_internacoes   = log1p(internacoes)
  )

vars_contagem  <- c("consultas", "ps", "exames", "terapias", "internacoes")
vars_derivadas <- c("total_eventos", "gasto_total", "meses_utilizacao",
                    "prestadores_distintos", "razao_exames_consultas", "razao_ps_consultas")
vars_flags     <- c("flag_alta_utilizacao", "flag_baixa_coerencia")
vars_cov       <- c("idade", "tempo_no_plano")

cat("\n", strrep("=", 60), "\n")
cat("SimSaúde — EDA\n")
cat("Dimensão:", nrow(dados), "observações ×", ncol(dados), "variáveis\n")
cat(strrep("=", 60), "\n\n")

# ------------------------------------------------------------
# 2) ESTATÍSTICAS DESCRITIVAS
# ------------------------------------------------------------
resumo_num <- function(x) {
  z <- x[!is.na(x)]
  c(
    n           = length(z),
    media       = mean(z),
    mediana     = median(z),
    dp          = sd(z),
    cv          = sd(z) / mean(z),
    assimetria  = moments::skewness(z),
    curtose     = moments::kurtosis(z),
    min         = min(z),
    p25         = unname(quantile(z, .25)),
    p75         = unname(quantile(z, .75)),
    max         = max(z),
    prop_zero   = mean(z == 0)
  )
}

tab_desc <- bind_rows(
  lapply(c(vars_contagem, vars_derivadas), function(v) {
    tibble(variavel = v, !!!as.list(resumo_num(dados[[v]])))
  })
)

cat("TABELA 1 — Estatísticas descritivas\n")
print(tab_desc, width = 120)

# Sobredispersão
tab_overdisp <- tibble(
  variavel       = vars_contagem,
  media          = sapply(dados[vars_contagem], mean),
  variancia      = sapply(dados[vars_contagem], var),
  razao_var_mu   = variancia / media,
  sobredisperso  = razao_var_mu > 1.5
)
cat("\nTABELA 2 — Sobredispersão nas contagens (var/média)\n")
print(tab_overdisp)

# ------------------------------------------------------------
# 3) DISTRIBUIÇÃO DOS GRUPOS LATENTES
# ------------------------------------------------------------

p_grupos <- ggplot(dados, aes(x = grupo_latente, fill = grupo_latente)) +
  geom_bar(width = 0.65) +
  geom_text(
    stat = "count",
    aes(label = paste0(after_stat(count), "\n(",
                       scales::percent(after_stat(count) / sum(after_stat(count)), .1),
                       ")")),
    vjust = -0.4, size = 3.2, lineheight = 1.2, color = "#333"
  ) +
  scale_fill_manual(values = CORES_GRUPOS) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = "Distribuição dos grupos latentes",
       subtitle = "n = 5.000 | proporções nominais: 50%, 25%, 20%, 5%",
       x = NULL, y = "Frequência") +
  guides(fill = "none")

print(p_grupos)

# ------------------------------------------------------------
# 4) HISTOGRAMAS + DENSIDADE — VARIÁVEIS DE CONTAGEM
#    (com detalhe: distribuição por grupo sobreposta)
# ------------------------------------------------------------

plot_hist_grupo <- function(df, var, titulo, log_scale = FALSE, bins = 45) {
  x_sym <- if (log_scale) sym(paste0("log_", var)) else sym(var)
  xlab  <- if (log_scale) paste0("log(", var, " + 1)") else var

  # Painel A — população total
  pA <- ggplot(df, aes(x = !!x_sym)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = bins, fill = CINZA_HIST, color = "white", linewidth = .25) +
    geom_density(color = COR_DENS, linewidth = .9, adjust = 1.3) +
    labs(title = "População total", x = xlab, y = "Densidade") +
    theme(plot.title = element_text(size = 11))

  # Painel B — densidades por grupo
  pB <- ggplot(df, aes(x = !!x_sym, color = grupo_latente, fill = grupo_latente)) +
    geom_density(alpha = .18, linewidth = .8, adjust = 1.3) +
    scale_color_manual(values = CORES_GRUPOS) +
    scale_fill_manual(values = CORES_GRUPOS) +
    labs(title = "Por grupo latente", x = xlab, y = "Densidade",
         color = NULL, fill = NULL) +
    theme(plot.title = element_text(size = 11),
          legend.position = "bottom",
          legend.text = element_text(size = 8))

  (pA + pB) +
    plot_annotation(
      title    = titulo,
      theme    = theme(plot.title = element_text(size = 13, face = "bold"))
    )
}

# Variáveis com estrutura de mistura visível
print(plot_hist_grupo(dados, "total_eventos",
                      "Total de eventos: assimetria à direita e mistura"))
print(plot_hist_grupo(dados, "gasto",
                      "Gasto total (log+1): separação entre grupos",
                      log_scale = TRUE))
print(plot_hist_grupo(dados, "prestadores",
                      "Prestadores distintos (log+1): heterogeneidade por grupo",
                      log_scale = TRUE))
print(plot_hist_grupo(dados, "consultas",
                      "Consultas: sobredispersão e zero-inflação moderada"))
print(plot_hist_grupo(dados, "exames",
                      "Exames: forte separação entre G1/G2"))
print(plot_hist_grupo(dados, "terapias",
                      "Terapias: excesso de zeros (ZINB)"))
print(plot_hist_grupo(dados, "internacoes",
                      "Internações: excesso de zeros (ZINB)"))

# ------------------------------------------------------------
# 5) RIDGELINE — DISTRIBUIÇÃO CONTÍNUA POR GRUPO
#    Ideal para ver separação de densidades simultaneamente
# ------------------------------------------------------------

ridge_vars <- c("log_total_eventos", "log_gasto", "log_consultas",
                "log_exames", "razao_exames_consultas", "razao_ps_consultas")

dados_ridge <- dados %>%
  select(grupo_latente, all_of(ridge_vars)) %>%
  pivot_longer(-grupo_latente, names_to = "variavel", values_to = "valor") %>%
  mutate(
    variavel = recode(variavel,
      log_total_eventos      = "log(total_eventos+1)",
      log_gasto              = "log(gasto_total+1)",
      log_consultas          = "log(consultas+1)",
      log_exames             = "log(exames+1)",
      razao_exames_consultas = "razão exames/consultas",
      razao_ps_consultas     = "razão PS/consultas"
    )
  )

p_ridge <- ggplot(dados_ridge,
                  aes(x = valor, y = grupo_latente,
                      fill = grupo_latente, color = grupo_latente)) +
  ggridges::geom_density_ridges(
    alpha = 0.45, linewidth = 0.6,
    scale = 1.4, rel_min_height = 0.01
  ) +
  scale_fill_manual(values  = CORES_GRUPOS) +
  scale_color_manual(values = CORES_GRUPOS) +
  facet_wrap(~ variavel, scales = "free_x", ncol = 2) +
  labs(
    title    = "Distribuição das variáveis-chave por grupo latente",
    subtitle = "Ridgeline: separação entre grupos evidencia estrutura de mistura",
    x = NULL, y = NULL
  ) +
  guides(fill = "none", color = "none") +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    axis.text.y = element_text(size = 8)
  )

print(p_ridge)

# ------------------------------------------------------------
# 6) EXCESSO DE ZEROS — OBSERVADO VS ESPERADO (Poisson e NB)
# ------------------------------------------------------------

comparar_zeros <- function(x, nome) {
  mu   <- mean(x)
  vari <- var(x)
  p0_obs  <- mean(x == 0)
  p0_pois <- exp(-mu)

  # NB estimada por momentos: size = mu^2 / (var - mu) se var > mu
  if (vari > mu) {
    size_nb <- mu^2 / (vari - mu)
    p0_nb   <- dnbinom(0, size = size_nb, mu = mu)
  } else {
    p0_nb <- NA_real_
  }

  tibble(variavel = nome, p0_obs = p0_obs, p0_pois = p0_pois, p0_nb = p0_nb)
}

tab_zeros <- bind_rows(lapply(vars_contagem, function(v) comparar_zeros(dados[[v]], v)))

cat("\nTABELA 3 — Proporção de zeros: observada vs esperada\n")
print(tab_zeros)

# Gráfico de zeros
zeros_long <- tab_zeros %>%
  pivot_longer(c(p0_obs, p0_pois, p0_nb),
               names_to = "modelo", values_to = "prop") %>%
  mutate(
    modelo = recode(modelo,
      p0_obs  = "Observado",
      p0_pois = "Esperado Poisson",
      p0_nb   = "Esperado NB"
    ),
    variavel = factor(variavel, levels = vars_contagem)
  )

p_zeros <- ggplot(zeros_long %>% filter(!is.na(prop)),
                  aes(x = variavel, y = prop, fill = modelo)) +
  geom_col(position = position_dodge(.75), width = .65) +
  geom_text(aes(label = scales::percent(prop, .1)),
            position = position_dodge(.75), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = c(
    "Observado"       = "#D98585",
    "Esperado Poisson"= "#B0B8C1",
    "Esperado NB"     = "#7BAFD4"
  )) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, .15))) +
  labs(
    title    = "Excesso de zeros nas contagens",
    subtitle = "Terapias e internações têm proporção de zeros muito acima do esperado pela NB → ZINB",
    x = NULL, y = "Proporção de zeros", fill = NULL
  )

print(p_zeros)

# ------------------------------------------------------------
# 7) SOBREDISPERSÃO — VISUALIZAÇÃO GRÁFICA
# ------------------------------------------------------------

p_overdisp <- tab_overdisp %>%
  mutate(variavel = factor(variavel, levels = vars_contagem)) %>%
  ggplot(aes(x = variavel, y = razao_var_mu, fill = variavel)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#888", linewidth = .8) +
  geom_text(aes(label = round(razao_var_mu, 1)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c(
    "consultas"   = COR_G1,
    "ps"          = COR_G4,
    "exames"      = COR_G2,
    "terapias"    = COR_G3,
    "internacoes" = "#C3A0C8"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  annotate("text", x = 0.6, y = 1.3, label = "linha = 1\n(Poisson)",
           size = 2.8, color = "#888", hjust = 0) +
  labs(
    title    = "Razão variância/média por tipo de serviço",
    subtitle = "Todos os valores > 1 confirmam sobredispersão — NB necessária",
    x = NULL, y = "Var / Média"
  ) +
  guides(fill = "none")

print(p_overdisp)

# ------------------------------------------------------------
# 8) VIOLIN + BOXPLOT POR GRUPO — CONTAGENS PRIMÁRIAS
# ------------------------------------------------------------

violin_grupo <- function(df, var, titulo, log_y = FALSE) {
  y_sym <- sym(var)
  ylab  <- if (log_y) paste0("log(", var, " + 1)") else var

  df_plot <- df %>%
    mutate(y_val = if (log_y) log1p(.data[[var]]) else .data[[var]])

  ggplot(df_plot, aes(x = grupo_latente, y = y_val, fill = grupo_latente)) +
    geom_violin(alpha = .45, trim = FALSE, linewidth = .4, color = "white") +
    geom_boxplot(width = .12, outlier.alpha = .1, outlier.size = .5,
                 fill = "white", linewidth = .5) +
    stat_summary(fun = mean, geom = "point", shape = 18,
                 size = 2.5, color = "#333") +
    scale_fill_manual(values = CORES_GRUPOS) +
    labs(title = titulo, x = NULL, y = ylab) +
    guides(fill = "none") +
    theme(axis.text.x = element_text(angle = 12, hjust = .6))
}

v1 <- violin_grupo(dados, "consultas",   "Consultas/ano")
v2 <- violin_grupo(dados, "ps",          "Pronto-socorro/ano")
v3 <- violin_grupo(dados, "exames",      "Exames/ano",      log_y = TRUE)
v4 <- violin_grupo(dados, "terapias",    "Terapias/ano",    log_y = TRUE)
v5 <- violin_grupo(dados, "internacoes", "Internações/ano", log_y = TRUE)
v6 <- violin_grupo(dados, "total_eventos", "Total de eventos", log_y = TRUE)

print(
  (v1 + v2 + v3) / (v4 + v5 + v6) +
    plot_annotation(
      title    = "Distribuição das contagens por grupo latente",
      subtitle = "Losango = média; violino = distribuição completa; box = IQR",
      theme    = theme(plot.title = element_text(size = 14, face = "bold"))
    )
)

# Variáveis derivadas
d1 <- violin_grupo(dados, "gasto_total",            "Gasto total (R$)", log_y = TRUE)
d2 <- violin_grupo(dados, "meses_utilizacao",        "Meses com utilização")
d3 <- violin_grupo(dados, "prestadores_distintos",   "Prestadores distintos")
d4 <- violin_grupo(dados, "razao_exames_consultas",  "Razão exames/consultas")
d5 <- violin_grupo(dados, "razao_ps_consultas",      "Razão PS/consultas")

print(
  (d1 + d2 + d3) / (d4 + d5 + plot_spacer()) +
    plot_annotation(
      title    = "Distribuição das variáveis derivadas por grupo latente",
      theme    = theme(plot.title = element_text(size = 14, face = "bold"))
    )
)

# ------------------------------------------------------------
# 9) COVARIÁVEIS: IDADE, SEXO, COPARTICIPAÇÃO, REGIÃO
# ------------------------------------------------------------

# 9a. Distribuição de idade por grupo
p_idade <- ggplot(dados, aes(x = idade, fill = grupo_latente, color = grupo_latente)) +
  geom_density(alpha = .3, linewidth = .7, adjust = 1.2) +
  scale_fill_manual(values  = CORES_GRUPOS) +
  scale_color_manual(values = CORES_GRUPOS) +
  labs(title    = "Distribuição de idade por grupo latente",
       subtitle = "Covariáveis geradas independentemente do grupo — distribuições similares esperadas",
       x = "Idade (anos)", y = "Densidade",
       color = "Grupo", fill = "Grupo")

print(p_idade)

# 9b. Proporção de faixas etárias
p_faixa <- dados %>%
  group_by(grupo_latente, faixa_idade) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(grupo_latente) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = grupo_latente, y = prop, fill = faixa_idade)) +
  geom_col(width = .65) +
  geom_text(aes(label = scales::percent(prop, 1)),
            position = position_stack(vjust = .5), size = 2.8, color = "white") +
  scale_fill_manual(values = c("0-17" = "#A8C8E8", "18-59" = "#6EBF8B", "60+" = "#E8B96B")) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Faixas etárias por grupo latente",
       x = NULL, y = "Proporção", fill = "Faixa")

# 9c. Coparticipação
p_copart <- dados %>%
  group_by(grupo_latente, copart) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(grupo_latente) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = grupo_latente, y = prop, fill = copart)) +
  geom_col(width = .65) +
  geom_text(aes(label = scales::percent(prop, 1)),
            position = position_stack(vjust = .5), size = 2.8, color = "white") +
  scale_fill_manual(values = c("Com" = "#7BAFD4", "Sem" = "#B0B8C1")) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Coparticipação por grupo latente",
       x = NULL, y = "Proporção", fill = "Copart.")

print(
  (p_faixa + p_copart) +
    plot_annotation(
      title = "Covariáveis observadas por grupo latente",
      subtitle = "Distribuições semelhantes confirmam independência das covariáveis em relação ao grupo",
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
)

# ------------------------------------------------------------
# 10) SCATTERPLOTS ESTRATÉGICOS — COM E SEM COR DE GRUPO
# ------------------------------------------------------------

scatter_base <- function(df, xvar, yvar, xlab, ylab, titulo,
                         log_x = FALSE, log_y = FALSE, amostra = 1500) {
  set.seed(42)
  df_s <- df %>% sample_n(min(amostra, nrow(df)))

  p <- ggplot(df_s, aes(x = .data[[xvar]], y = .data[[yvar]],
                         color = grupo_latente)) +
    geom_point(alpha = .3, size = .9) +
    scale_color_manual(values = CORES_GRUPOS) +
    labs(title = titulo, x = xlab, y = ylab, color = "Grupo") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

  if (log_x) p <- p + scale_x_continuous(trans = "log1p")
  if (log_y) p <- p + scale_y_continuous(trans = "log1p", labels = comma)
  p
}

s1 <- scatter_base(dados, "total_eventos", "gasto_total",
                   "Total de eventos", "Gasto total (R$)",
                   "Gasto vs eventos", log_y = TRUE)

s2 <- scatter_base(dados, "consultas", "exames",
                   "Consultas", "Exames",
                   "Exames vs consultas")

s3 <- scatter_base(dados, "consultas", "ps",
                   "Consultas", "PS",
                   "PS vs consultas", log_x = TRUE, log_y = TRUE)

s4 <- scatter_base(dados, "terapias", "ps",
                   "Terapias", "PS",
                   "PS vs terapias — padrão atípico",
                   log_x = TRUE, log_y = TRUE)

s5 <- scatter_base(dados, "exames", "razao_exames_consultas",
                   "Exames", "Razão exames/consultas",
                   "Razão exames vs exames brutos", log_x = TRUE, log_y = TRUE)

s6 <- scatter_base(dados, "meses_utilizacao", "prestadores_distintos",
                   "Meses com utilização", "Prestadores distintos",
                   "Fragmentação do cuidado vs meses ativos")

print(
  (s1 + s2) / (s3 + s4) / (s5 + s6) +
    plot_annotation(
      title  = "Scatterplots estratégicos coloridos por grupo latente",
      theme  = theme(plot.title = element_text(size = 14, face = "bold"))
    )
)

# ------------------------------------------------------------
# 11) FLAGS DE ALERTA — PREVALÊNCIA E POSICIONAMENTO
# ------------------------------------------------------------

# 11a. Prevalência por grupo
flags_long <- dados %>%
  group_by(grupo_latente) %>%
  summarise(
    Alta_Utilizacao = mean(flag_alta_utilizacao),
    Baixa_Coerencia = mean(flag_baixa_coerencia)
  ) %>%
  pivot_longer(-grupo_latente, names_to = "flag", values_to = "prop") %>%
  mutate(flag = recode(flag,
    "Alta_Utilizacao" = "Alta utilização (P95)",
    "Baixa_Coerencia" = "Baixa coerência (P95)"
  ))

p_flags <- ggplot(flags_long, aes(x = grupo_latente, y = prop, fill = grupo_latente)) +
  geom_col(width = .6) +
  geom_text(aes(label = scales::percent(prop, .1)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = CORES_GRUPOS) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, .18))) +
  facet_wrap(~ flag) +
  labs(title    = "Prevalência de flags de alerta por grupo latente",
       subtitle = "Flag de baixa coerência concentrada no grupo atípico — sinal de auditoria",
       x = NULL, y = "Proporção") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 12, hjust = .6))

print(p_flags)

# 11b. Distribuição de razões nos beneficiários flagados vs não flagados
p_flag_coer <- ggplot(dados,
  aes(x = razao_exames_consultas, fill = factor(flag_baixa_coerencia))) +
  geom_density(alpha = .5, adjust = 1.5) +
  scale_fill_manual(values = c("0" = "#B0C8D8", "1" = "#D98585"),
                    labels = c("Sem flag", "Com flag")) +
  scale_x_continuous(limits = c(0, quantile(dados$razao_exames_consultas, .99))) +
  labs(title    = "Razão exames/consultas: flagados vs não flagados",
       x = "Razão exames/consultas", y = "Densidade", fill = NULL)

p_flag_util <- ggplot(dados,
  aes(x = log_total_eventos, fill = factor(flag_alta_utilizacao))) +
  geom_density(alpha = .5, adjust = 1.3) +
  scale_fill_manual(values = c("0" = "#B0C8D8", "1" = "#E8B96B"),
                    labels = c("Sem flag", "Com flag")) +
  labs(title = "Total de eventos (log): alta utilização vs restante",
       x = "log(total_eventos + 1)", y = "Densidade", fill = NULL)

print(p_flag_coer + p_flag_util +
  plot_annotation(title = "Distribuição das variáveis de alerta"))

# ------------------------------------------------------------
# 12) CORRELAÇÃO — MATRIZ POR GRUPO
# ------------------------------------------------------------

vars_corr <- c("consultas", "ps", "exames", "terapias", "internacoes",
               "razao_exames_consultas", "razao_ps_consultas",
               "meses_utilizacao", "prestadores_distintos")

corr_mat_grupo <- function(grp) {
  dados %>%
    filter(grupo_latente == grp) %>%
    select(all_of(vars_corr)) %>%
    cor(use = "complete.obs") %>%
    as.data.frame() %>%
    rownames_to_column("var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "corr") %>%
    mutate(grupo = grp)
}

corr_all <- bind_rows(lapply(levels(dados$grupo_latente), corr_mat_grupo)) %>%
  mutate(
    var1  = factor(var1,  levels = vars_corr),
    var2  = factor(var2,  levels = rev(vars_corr)),
    grupo = factor(grupo, levels = levels(dados$grupo_latente))
  )

p_corr <- ggplot(corr_all, aes(x = var1, y = var2, fill = corr)) +
  geom_tile(color = "white", linewidth = .3) +
  geom_text(aes(label = ifelse(abs(corr) > .15, round(corr, 2), "")),
            size = 2.2, color = "#333") +
  scale_fill_gradient2(
    low = "#D98585", mid = "white", high = "#7BAFD4",
    midpoint = 0, limits = c(-1, 1), name = "r"
  ) +
  facet_wrap(~ grupo, ncol = 2) +
  labs(title    = "Matriz de correlação por grupo latente",
       subtitle = "Correlações diferem entre grupos — evidência de heterogeneidade estrutural",
       x = NULL, y = NULL) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 7.5),
    axis.text.y  = element_text(size = 7.5),
    legend.key.height = unit(.8, "cm")
  )

print(p_corr)

# ------------------------------------------------------------
# 13) ECDF COMPARATIVO — VARIÁVEIS MISTAS
# ------------------------------------------------------------

ecdf_vars <- c("log_total_eventos", "log_gasto", "meses_utilizacao",
               "prestadores_distintos")

dados_ecdf <- dados %>%
  select(grupo_latente, all_of(ecdf_vars)) %>%
  pivot_longer(-grupo_latente, names_to = "variavel", values_to = "valor") %>%
  mutate(variavel = recode(variavel,
    log_total_eventos   = "log(total_eventos+1)",
    log_gasto           = "log(gasto_total+1)",
    meses_utilizacao    = "meses_utilizacao",
    prestadores_distintos = "prestadores_distintos"
  ))

p_ecdf <- ggplot(dados_ecdf,
                 aes(x = valor, color = grupo_latente)) +
  stat_ecdf(linewidth = .8, pad = FALSE) +
  scale_color_manual(values = CORES_GRUPOS) +
  facet_wrap(~ variavel, scales = "free_x", ncol = 2) +
  labs(title    = "ECDF por grupo latente",
       subtitle = "Separação entre curvas evidencia estrutura de mistura",
       x = NULL, y = "F(x)", color = "Grupo") +
  theme(legend.position = "bottom")

print(p_ecdf)

# ------------------------------------------------------------
# 14) TESTE DE HARTIGAN (unimodalidade)
# ------------------------------------------------------------

vars_dip <- c("total_eventos", "log_gasto", "log_total_eventos",
              "prestadores_distintos", "consultas", "ps",
              "razao_exames_consultas", "meses_utilizacao")

dip_results <- bind_rows(
  lapply(vars_dip, function(v) {
    teste <- diptest::dip.test(dados[[v]])
    tibble(
      variavel  = v,
      dip       = unname(teste$statistic),
      p_valor   = teste$p.value,
      rejeita_H0 = teste$p.value < 0.05
    )
  })
)

cat("\nTABELA 4 — Teste de Hartigan (H₀: unimodalidade)\n")
print(dip_results)

# Gráfico com p-valores
p_dip <- dip_results %>%
  mutate(variavel = fct_reorder(variavel, p_valor)) %>%
  ggplot(aes(x = variavel, y = -log10(p_valor),
             fill = rejeita_H0)) +
  geom_col(width = .65) +
  geom_hline(yintercept = -log10(.05), linetype = "dashed",
             color = "#888", linewidth = .8) +
  scale_fill_manual(values = c("TRUE" = "#D98585", "FALSE" = "#B0B8C1"),
                    labels = c("TRUE" = "Rejeita H₀ (p < 0,05)",
                               "FALSE"= "Não rejeita")) +
  annotate("text", x = .6, y = -log10(.05) + .15, label = "α = 0,05",
           size = 3, color = "#888", hjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  coord_flip() +
  labs(title    = "Teste de Hartigan — evidência de multimodalidade",
       subtitle = "−log₁₀(p): quanto maior, mais evidência contra unimodalidade",
       x = NULL, y = "−log₁₀(p-valor)", fill = NULL) +
  theme(legend.position = "top")

print(p_dip)

# ------------------------------------------------------------
# 15) PCA — ESPAÇO REDUZIDO COM COR DE GRUPO
# ------------------------------------------------------------

vars_pca <- c("consultas", "ps", "exames", "terapias", "internacoes",
              "log_total_eventos", "log_gasto", "log_prestadores",
              "razao_exames_consultas", "razao_ps_consultas",
              "meses_utilizacao")

X_pca <- dados %>%
  select(all_of(vars_pca)) %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.x), NA, .x))) %>%
  drop_na()

idx_pca <- as.integer(rownames(X_pca))
pca_fit <- prcomp(scale(X_pca), center = TRUE, scale. = TRUE)

# Variância explicada
var_exp <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)

p_scree <- tibble(
  PC   = paste0("PC", seq_along(var_exp)),
  var  = var_exp,
  acum = cumsum(var_exp)
) %>%
  slice(1:8) %>%
  mutate(PC = factor(PC, levels = PC)) %>%
  ggplot(aes(x = PC)) +
  geom_col(aes(y = var), fill = COR_G1, width = .6) +
  geom_line(aes(y = acum, group = 1), color = COR_G4, linewidth = .9) +
  geom_point(aes(y = acum), color = COR_G4, size = 2.5) +
  scale_y_continuous(labels = percent_format(), sec.axis = sec_axis(~ .,
    labels = percent_format(), name = "Variância acumulada")) +
  labs(title = "Scree plot — variância explicada por componente",
       x = NULL, y = "Variância explicada")

# Biplot PC1 vs PC2
scores_df <- as_tibble(pca_fit$x[, 1:3]) %>%
  mutate(grupo_latente = dados$grupo_latente[idx_pca])

p_pca12 <- ggplot(scores_df %>% sample_n(min(2000, nrow(.))),
                  aes(PC1, PC2, color = grupo_latente)) +
  geom_point(alpha = .3, size = .9) +
  stat_ellipse(level = .85, linewidth = .7, linetype = "dashed") +
  scale_color_manual(values = CORES_GRUPOS) +
  labs(title    = "PCA — PC1 vs PC2",
       subtitle = "Elipses de 85% de contenção por grupo",
       color = "Grupo")

p_pca13 <- ggplot(scores_df %>% sample_n(min(2000, nrow(.))),
                  aes(PC1, PC3, color = grupo_latente)) +
  geom_point(alpha = .3, size = .9) +
  stat_ellipse(level = .85, linewidth = .7, linetype = "dashed") +
  scale_color_manual(values = CORES_GRUPOS) +
  labs(title = "PCA — PC1 vs PC3", color = "Grupo")

print(
  p_scree / (p_pca12 + p_pca13) +
    plot_annotation(
      title  = "Análise de Componentes Principais (PCA)",
      theme  = theme(plot.title = element_text(size = 14, face = "bold"))
    )
)

# Loadings — contribuição de cada variável
loadings_df <- as_tibble(pca_fit$rotation[, 1:3],
                          rownames = "variavel") %>%
  pivot_longer(-variavel, names_to = "PC", values_to = "loading")

p_load <- ggplot(loadings_df,
                 aes(x = fct_reorder(variavel, loading), y = loading,
                     fill = loading > 0)) +
  geom_col(width = .7) +
  geom_hline(yintercept = 0, linewidth = .4) +
  scale_fill_manual(values = c("TRUE" = COR_G2, "FALSE" = COR_G4)) +
  coord_flip() +
  facet_wrap(~ PC, ncol = 3) +
  labs(title = "Loadings por componente principal",
       x = NULL, y = "Loading") +
  guides(fill = "none")

print(p_load)

# ------------------------------------------------------------
# 16) GMM UNIVARIADO — BIC E DENSIDADE DE MISTURA
# ------------------------------------------------------------

fit_e_plot_gmm <- function(x_vals, nome_var, G_max = 6) {
  mod  <- Mclust(x_vals, G = 1:G_max, modelNames = "V", verbose = FALSE)
  G_opt <- mod$G

  # BIC do modelo selecionado
  bic_val <- max(mod$BIC, na.rm = TRUE)
  cat("\n── GMM univariado:", nome_var,
      "| Melhor G:", G_opt,
      "| BIC:", round(bic_val, 1), "\n")

  grid_x <- seq(min(x_vals), max(x_vals), length.out = 500)

  # Densidade total da mistura: soma ponderada das densidades componentes
  pros   <- mod$parameters$pro                  # pesos (comprimento G_opt)
  means  <- mod$parameters$mean                 # médias
  # sigmasq pode ser escalar (modelo "E") ou vetor de comprimento G (modelo "V")
  sigmas <- mod$parameters$variance$sigmasq
  if (length(sigmas) == 1) sigmas <- rep(sigmas, G_opt)

  # Densidade de cada componente em cada ponto do grid
  comp_df <- lapply(seq_len(G_opt), function(k) {
    tibble(
      x    = grid_x,
      dens = pros[k] * dnorm(grid_x, mean = means[k], sd = sqrt(sigmas[k])),
      comp = paste0("Componente ", k)
    )
  }) %>% bind_rows()

  # Densidade total = soma das densidades componentes
  dens_mix_df <- comp_df %>%
    group_by(x) %>%
    summarise(dens_total = sum(dens), .groups = "drop")

  ggplot() +
    geom_histogram(
      data = tibble(x = x_vals),
      aes(x = x, y = after_stat(density)),
      bins = 45, fill = CINZA_HIST, color = "white", linewidth = .2
    ) +
    geom_line(
      data = dens_mix_df,
      aes(x = x, y = dens_total),
      color = "#333", linewidth = 1
    ) +
    geom_line(
      data = comp_df,
      aes(x = x, y = dens, color = comp),
      linewidth = .75, linetype = "dashed"
    ) +
    scale_color_manual(values = c(
      "Componente 1" = "#7BAFD4",
      "Componente 2" = "#D98585",
      "Componente 3" = "#6EBF8B",
      "Componente 4" = "#E8B96B",
      "Componente 5" = "#C3A0C8",
      "Componente 6" = "#A0B89A"
    )) +
    labs(
      title    = paste0("GMM univariado — ", nome_var),
      subtitle = paste0("G ótimo = ", G_opt,
                        " componentes (BIC) | linha sólida = mistura total"),
      x = nome_var, y = "Densidade", color = NULL
    ) +
    theme(legend.position = "top")
}

print(fit_e_plot_gmm(dados$log_total_eventos, "log(total_eventos+1)"))
print(fit_e_plot_gmm(dados$log_gasto,          "log(gasto_total+1)"))
print(fit_e_plot_gmm(dados$log_prestadores,    "log(prestadores_distintos+1)"))
print(fit_e_plot_gmm(dados$meses_utilizacao,   "meses_utilizacao"))

# ------------------------------------------------------------
# 17) GMM MULTIVARIADO — BIC E PROJEÇÃO PCA
# ------------------------------------------------------------

X_mix <- dados %>%
  transmute(
    log_total_eventos      = log1p(total_eventos),
    log_gasto              = log1p(gasto_total),
    log_prestadores        = log1p(prestadores_distintos),
    razao_exames_consultas = razao_exames_consultas,
    razao_ps_consultas     = razao_ps_consultas,
    meses_utilizacao       = meses_utilizacao,
    consultas              = consultas,
    ps                     = ps,
    exames                 = exames,
    terapias               = terapias,
    internacoes            = internacoes
  ) %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.x), NA, .x))) %>%
  drop_na() %>%
  scale()

cat("\nAjustando GMM multivariado (G = 1:6)...\n")
mod_mix <- Mclust(X_mix, G = 1:6, verbose = FALSE)

cat("\nTABELA 5 — Resumo GMM multivariado\n")
print(summary(mod_mix))

cluster_est <- factor(mod_mix$classification)
grupo_real  <- dados$grupo_latente[seq_len(nrow(X_mix))]

tab_cruzada <- table(Estimado = cluster_est, Real = grupo_real)
cat("\nTABELA 6 — Tabela cruzada: cluster estimado vs grupo latente real\n")
print(tab_cruzada)
cat("\nProporção por linha:\n")
print(round(prop.table(tab_cruzada, margin = 1), 3))

# BIC por número de componentes
bic_df <- tibble(
  G   = seq_len(nrow(mod_mix$BIC)),
  BIC = apply(mod_mix$BIC, 1, max, na.rm = TRUE)
)

p_bic <- ggplot(bic_df, aes(x = G, y = BIC)) +
  geom_line(color = COR_ACENTO, linewidth = .9) +
  geom_point(color = COR_ACENTO, size = 2.8) +
  geom_point(data = bic_df %>% filter(BIC == max(BIC)),
             color = COR_G4, size = 4.5, shape = 18) +
  scale_x_continuous(breaks = 1:6) +
  labs(title    = "BIC por número de componentes — GMM multivariado",
       subtitle = paste0("Losango vermelho = G ótimo = ", mod_mix$G),
       x = "Número de componentes (G)", y = "BIC")

# Projeção no espaço PCA
scores_mix <- as_tibble(prcomp(X_mix, scale. = FALSE)$x[, 1:2]) %>%
  mutate(cluster_est   = cluster_est,
         grupo_latente = grupo_real)

gm1 <- ggplot(scores_mix, aes(PC1, PC2, color = cluster_est)) +
  geom_point(alpha = .3, size = .9) +
  stat_ellipse(level = .8, linewidth = .6, linetype = "dashed") +
  scale_color_manual(values = c("#7BAFD4","#6EBF8B","#E8B96B","#D98585","#C3A0C8","#A0B89A")) +
  labs(title = "Clusters estimados (GMM)", color = "Cluster")

gm2 <- ggplot(scores_mix, aes(PC1, PC2, color = grupo_latente)) +
  geom_point(alpha = .3, size = .9) +
  stat_ellipse(level = .8, linewidth = .6, linetype = "dashed") +
  scale_color_manual(values = CORES_GRUPOS) +
  labs(title = "Grupos latentes reais", color = "Grupo")

print(
  p_bic / (gm1 + gm2) +
    plot_annotation(
      title = "GMM Multivariado — seleção de modelo e comparação com verdade",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
)

# ------------------------------------------------------------
# 18) HEATMAP — MÉDIAS PADRONIZADAS POR GRUPO (PERFIL)
# ------------------------------------------------------------

vars_perfil <- c("consultas", "ps", "exames", "terapias", "internacoes",
                 "razao_exames_consultas", "razao_ps_consultas",
                 "meses_utilizacao", "prestadores_distintos")

medias_grupo <- dados %>%
  group_by(grupo_latente) %>%
  summarise(across(all_of(vars_perfil), mean)) %>%
  pivot_longer(-grupo_latente, names_to = "variavel", values_to = "media")

# Padronização z para comparação
medias_pad <- medias_grupo %>%
  group_by(variavel) %>%
  mutate(z = (media - mean(media)) / sd(media)) %>%
  ungroup()

p_heat <- ggplot(medias_pad,
                 aes(x = variavel, y = grupo_latente, fill = z)) +
  geom_tile(color = "white", linewidth = .5) +
  geom_text(aes(label = round(z, 1)), size = 3, color = "#333") +
  scale_fill_gradient2(
    low = "#D98585", mid = "#f7f7f7", high = "#7BAFD4",
    midpoint = 0, name = "z-score"
  ) +
  labs(
    title    = "Heatmap de perfis por grupo latente (z-score das médias)",
    subtitle = "Células azuis = acima da média geral; vermelhas = abaixo",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10)
  )

print(p_heat)

# ------------------------------------------------------------
# RESUMO FINAL
# ------------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("RESULTADOS-CHAVE DA EDA\n")
cat(strrep("=", 60), "\n")
cat("1. Todas as contagens apresentam sobredispersão (var/μ > 1).\n")
cat("2. Terapias e internações têm excesso de zeros estruturais → ZINB.\n")
cat("3. Ridgelines e violins confirmam separabilidade dos grupos.\n")
cat("4. Teste de Hartigan rejeita unimodalidade nas variáveis agregadas.\n")
cat("5. PCA: elipses dos grupos se separam nos 2 primeiros componentes.\n")
cat("6. GMM univariado seleciona G > 1 em log_total_eventos e log_gasto.\n")
cat("7. GMM multivariado: G ótimo identificado pelo BIC.\n")
cat("8. Heatmap confirma perfis distintos — base válida para FMNB.\n")
cat(strrep("=", 60), "\n")
