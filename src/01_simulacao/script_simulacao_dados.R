# ============================================================
# SIMULACAO DE DADOS DE UTILIZACAO EM SAUDE SUPLEMENTAR
# Katy - estudo preliminar para poster
# ============================================================
# VERSAO ATUALIZADA:
#  1. Geracao vetorizada de idades
#  2. Coeficientes grupo-especificos
#  3. meses_utilizacao com Beta-Binomial de fato
#  4. prestadores_distintos com cap inteiro
#  5. Diagnostico do gasto_total com log-transformacao
#  6. Validacao com medias simuladas e medias teoricas (mu)
#  7. Grafico adicional: gasto x eventos por grupo
#  8. Tabela de prevalencia das flags por grupo
#  9. Intercepto de exames do grupo atipico suavizado
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)

set.seed(1234)

# ------------------------------------------------------------
# 1) TAMANHO DA AMOSTRA
# ------------------------------------------------------------
n <- 5000

# ------------------------------------------------------------
# 2) PROPORCOES DOS GRUPOS LATENTES
# ------------------------------------------------------------
# 1 = baixo uso
# 2 = uso ambulatorial coordenado
# 3 = uso agudo/hospitalar
# 4 = padrao atipico

pi_grupos <- c(0.50, 0.25, 0.20, 0.05)

grupo_latente <- sample(
  x       = 1:4,
  size    = n,
  replace = TRUE,
  prob    = pi_grupos
)

# ------------------------------------------------------------
# 3) COVARIAVEIS OBSERVADAS
# ------------------------------------------------------------

# --- Idade: geracao vetorizada ------------------------------
faixa_idx <- sample(1:3, size = n, replace = TRUE, prob = c(0.20, 0.55, 0.25))

pool_0_17  <- 0:17
pool_18_59 <- 18:59
pool_60_85 <- 60:85

idade <- integer(n)
idade[faixa_idx == 1] <- sample(pool_0_17,  sum(faixa_idx == 1), replace = TRUE)
idade[faixa_idx == 2] <- sample(pool_18_59, sum(faixa_idx == 2), replace = TRUE)
idade[faixa_idx == 3] <- sample(pool_60_85, sum(faixa_idx == 3), replace = TRUE)

faixa_idade <- c("0-17", "18-59", "60+")[faixa_idx]
idade10     <- idade / 10

# --- Sexo: 1 = feminino, 0 = masculino ---------------------
sexo <- rbinom(n, size = 1, prob = 0.55)

# --- Tempo no plano em meses -------------------------------
tempo_no_plano <- sample(1:120, size = n, replace = TRUE)
tempo_anos     <- tempo_no_plano / 12

# --- Coparticipacao: 1 = com, 0 = sem ----------------------
copart <- rbinom(n, size = 1, prob = 0.40)

# --- Regiao ------------------------------------------------
regiao  <- sample(x = 1:3, size = n, replace = TRUE, prob = c(0.40, 0.35, 0.25))
regiao2 <- ifelse(regiao == 2, 1, 0)
regiao3 <- ifelse(regiao == 3, 1, 0)

# ------------------------------------------------------------
# 4) PARAMETROS DOS GRUPOS (interceptos no log da media)
# ------------------------------------------------------------
# Recalibrados apos comparar medias simuladas com alvos ANS

alpha_consultas   <- log(c(2.5,  8.0,  4.0,  1.5))  - 0.588
alpha_ps          <- log(c(0.2,  0.5,  4.0,  7.0))  - 0.321
alpha_exames      <- log(c(8.0, 50.0, 20.0, 52.0))  - 0.607
alpha_terapias    <- log(c(0.1,  3.0,  0.7,  6.5))  + 0.282
alpha_internacoes <- log(c(0.04, 0.06, 0.65, 0.35)) + 0.250



# ------------------------------------------------------------
# 5) COEFICIENTES GRUPO-ESPECIFICOS
# ------------------------------------------------------------

betas <- list(
  
  consultas = matrix(
    c(
      # idade10  sexo   tempo  copart  regiao2  regiao3
      0.08,    0.12,  0.03, -0.18,   0.05,    0.00,   # G1
      0.10,    0.14,  0.04, -0.10,   0.05,    0.00,   # G2
      0.07,    0.10,  0.02, -0.05,   0.05,    0.00,   # G3
      -0.05,    0.08,  0.01, -0.02,   0.00,    0.05    # G4
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("idade10","sexo","tempo","copart","regiao2","regiao3"))
  ),
  
  ps = matrix(
    c(
      0.03,   -0.03,  0.00, -0.10,   0.00,    0.08,
      0.04,   -0.02,  0.00, -0.06,   0.00,    0.08,
      0.06,    0.00,  0.01, -0.03,   0.00,    0.08,
      0.02,    0.05,  0.01, -0.01,   0.05,    0.10
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("idade10","sexo","tempo","copart","regiao2","regiao3"))
  ),
  
  exames = matrix(
    c(
      0.10,    0.08,  0.02, -0.15,   0.05,    0.00,
      0.12,    0.09,  0.03, -0.08,   0.05,    0.00,
      0.08,    0.06,  0.02, -0.04,   0.05,    0.00,
      -0.03,    0.05,  0.00, -0.01,   0.00,    0.03
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("idade10","sexo","tempo","copart","regiao2","regiao3"))
  ),
  
  terapias = matrix(
    c(
      0.05,    0.10,  0.01, -0.12,   0.00,    0.00,
      0.06,    0.12,  0.02, -0.07,   0.00,    0.00,
      0.04,    0.08,  0.01, -0.04,   0.00,    0.00,
      0.02,    0.15,  0.01, -0.02,   0.00,    0.05
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("idade10","sexo","tempo","copart","regiao2","regiao3"))
  ),
  
  internacoes = matrix(
    c(
      0.12,    0.00,  0.01, -0.03,   0.00,    0.08,
      0.13,    0.02,  0.01, -0.02,   0.00,    0.08,
      0.15,    0.03,  0.02, -0.01,   0.00,    0.08,
      0.05,    0.02,  0.00, -0.01,   0.05,    0.10
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(NULL, c("idade10","sexo","tempo","copart","regiao2","regiao3"))
  )
)

# ------------------------------------------------------------
# 6) PARAMETROS DE DISPERSAO (BINOMIAL NEGATIVA)
# ------------------------------------------------------------
theta_consultas   <- 2.5
theta_ps          <- 1.5
theta_exames      <- 2.0
theta_terapias    <- 1.2
theta_internacoes <- 0.8

# ------------------------------------------------------------
# 7) PROBABILIDADES DE ZERO ESTRUTURAL (ZINB)
# ------------------------------------------------------------
pzero_terapias    <- c(0.85, 0.45, 0.70, 0.35)
pzero_internacoes <- c(0.97, 0.94, 0.55, 0.75)

# ------------------------------------------------------------
# 8) FUNCOES AUXILIARES
# ------------------------------------------------------------

calc_mu <- function(alpha_vec, beta_mat, grupo, idade10, sexo,
                    tempo_anos, copart, regiao2, regiao3) {
  
  X <- cbind(idade10, sexo, tempo_anos, copart, regiao2, regiao3)
  colnames(X) <- c("idade10","sexo","tempo","copart","regiao2","regiao3")
  
  beta_ind <- beta_mat[grupo, ]
  eta <- alpha_vec[grupo] + rowSums(X * beta_ind)
  
  exp(eta)
}

rnbinom_mu <- function(mu, theta) {
  rnbinom(n = length(mu), size = theta, mu = mu)
}

rzinb_mu <- function(mu, theta, pzero) {
  zero_flag <- rbinom(length(mu), size = 1, prob = pzero)
  y_nb      <- rnbinom(n = length(mu), size = theta, mu = mu)
  ifelse(zero_flag == 1, 0, y_nb)
}

# ------------------------------------------------------------
# 9) MEDIAS CONDICIONAIS
# ------------------------------------------------------------

mu_consultas <- calc_mu(
  alpha_consultas, betas$consultas,
  grupo_latente, idade10, sexo, tempo_anos, copart, regiao2, regiao3
)

mu_ps <- calc_mu(
  alpha_ps, betas$ps,
  grupo_latente, idade10, sexo, tempo_anos, copart, regiao2, regiao3
)

mu_exames <- calc_mu(
  alpha_exames, betas$exames,
  grupo_latente, idade10, sexo, tempo_anos, copart, regiao2, regiao3
)

mu_terapias <- calc_mu(
  alpha_terapias, betas$terapias,
  grupo_latente, idade10, sexo, tempo_anos, copart, regiao2, regiao3
)

mu_internacoes <- calc_mu(
  alpha_internacoes, betas$internacoes,
  grupo_latente, idade10, sexo, tempo_anos, copart, regiao2, regiao3
)

# ------------------------------------------------------------
# 10) CONTAGENS
# ------------------------------------------------------------

consultas <- rnbinom_mu(mu_consultas, theta_consultas)
ps        <- rnbinom_mu(mu_ps, theta_ps)
exames    <- rnbinom_mu(mu_exames, theta_exames)

terapias <- rzinb_mu(
  mu    = mu_terapias,
  theta = theta_terapias,
  pzero = pzero_terapias[grupo_latente]
)

internacoes <- rzinb_mu(
  mu    = mu_internacoes,
  theta = theta_internacoes,
  pzero = pzero_internacoes[grupo_latente]
)

# ------------------------------------------------------------
# 11) VARIAVEIS DERIVADAS
# ------------------------------------------------------------

total_eventos <- consultas + ps + exames + terapias + internacoes

razao_exames_consultas <- exames / (consultas + 1)
razao_ps_consultas     <- ps / (consultas + 1)

# --- Gasto total simulado -----------------------------------
# Custos medios por evento (R$):
# consulta = 120 | ps = 250 | exame = 80 | terapia = 150 | internacao = 3500
# Ruido log-normal com CV aproximado de 20%

gasto_base <- consultas   * 120 +
  ps          * 250 +
  exames      * 80  +
  terapias    * 150 +
  internacoes * 3500

ruido_gasto <- rlnorm(n, meanlog = 0, sdlog = 0.20)
gasto_total <- round(gasto_base * ruido_gasto, 2)

# --- meses_utilizacao: Beta-Binomial de fato ----------------
shape_a <- pmax(0.5, total_eventos / 6)
shape_b <- pmax(0.5, (36 - total_eventos) / 6)

p_mes <- rbeta(n, shape1 = shape_a, shape2 = shape_b)
meses_utilizacao <- rbinom(n, size = 12, prob = p_mes)
meses_utilizacao <- pmax(1L, meses_utilizacao)

# --- prestadores_distintos com cap inteiro ------------------
prest_raw <- pmax(1L, rpois(n, lambda = pmax(1, total_eventos / 5)))
cap_prest <- as.integer(quantile(prest_raw, 0.99))
prestadores_distintos <- pmin(prest_raw, cap_prest)

# --- Flags --------------------------------------------------
flag_alta_utilizacao <- ifelse(
  total_eventos >= quantile(total_eventos, 0.95),
  1L, 0L
)

flag_baixa_coerencia <- ifelse(
  razao_exames_consultas > quantile(razao_exames_consultas, 0.95) |
    razao_ps_consultas   > quantile(razao_ps_consultas, 0.95),
  1L, 0L
)

# ------------------------------------------------------------
# 12) BASE FINAL
# ------------------------------------------------------------

dados_sim <- data.frame(
  id                     = 1:n,
  grupo_latente          = factor(
    grupo_latente,
    levels = 1:4,
    labels = c(
      "baixo_uso",
      "ambulatorial_coordenado",
      "agudo_hospitalar",
      "atipico"
    )
  ),
  idade                  = idade,
  faixa_idade            = factor(faixa_idade, levels = c("0-17","18-59","60+")),
  sexo                   = factor(sexo, levels = c(0,1), labels = c("Masculino","Feminino")),
  tempo_no_plano         = tempo_no_plano,
  copart                 = factor(copart, levels = c(0,1), labels = c("Sem","Com")),
  regiao                 = factor(regiao),
  consultas              = consultas,
  ps                     = ps,
  exames                 = exames,
  terapias               = terapias,
  internacoes            = internacoes,
  total_eventos          = total_eventos,
  razao_exames_consultas = razao_exames_consultas,
  razao_ps_consultas     = razao_ps_consultas,
  gasto_total            = gasto_total,
  meses_utilizacao       = meses_utilizacao,
  prestadores_distintos  = prestadores_distintos,
  flag_alta_utilizacao   = flag_alta_utilizacao,
  flag_baixa_coerencia   = flag_baixa_coerencia
)


# ------------------------------------------------------------
# CHECAGEM DOS AGREGADOS VS ALVOS ANS
# ------------------------------------------------------------

alvos_ans <- c(
  consultas   = 4.15,
  ps          = 1.34,
  exames      = 22.9,
  terapias    = 1.30,
  internacoes = 0.182
)

medias_obs <- c(
  consultas   = mean(consultas),
  ps          = mean(ps),
  exames      = mean(exames),
  terapias    = mean(terapias),
  internacoes = mean(internacoes)
)

comparacao_ans <- data.frame(
  variavel = names(alvos_ans),
  alvo_ans = as.numeric(alvos_ans),
  media_observada = as.numeric(medias_obs),
  diferenca = as.numeric(medias_obs - alvos_ans),
  diferenca_percentual = 100 * as.numeric((medias_obs - alvos_ans) / alvos_ans)
)

cat("\nComparacao com alvos ANS:\n")
print(
  comparacao_ans %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
)

# ------------------------------------------------------------
# 13) RESUMOS
# ------------------------------------------------------------

cat("\nDimensao da base:\n")
print(dim(dados_sim))

cat("\nPrimeiras linhas:\n")
print(head(dados_sim))

cat("\nDistribuicao dos grupos latentes:\n")
print(prop.table(table(dados_sim$grupo_latente)))

resumo_grupo <- dados_sim %>%
  group_by(grupo_latente) %>%
  summarise(
    n                          = n(),
    idade_media                = mean(idade),
    consultas_media            = mean(consultas),
    ps_media                   = mean(ps),
    exames_media               = mean(exames),
    terapias_media             = mean(terapias),
    internacoes_media          = mean(internacoes),
    total_eventos_medio        = mean(total_eventos),
    gasto_total_medio          = mean(gasto_total),
    razao_exames_consultas_med = mean(razao_exames_consultas),
    razao_ps_consultas_med     = mean(razao_ps_consultas),
    meses_utilizacao_medio     = mean(meses_utilizacao),
    .groups = "drop"
  )

cat("\nMedias observadas por grupo:\n")
print(resumo_grupo)

# ------------------------------------------------------------
# TABELA: medias das cinco contagens por grupo
# ------------------------------------------------------------

tabela_contagens <- dados_sim %>%
  group_by(grupo_latente) %>%
  summarise(
    n = n(),
    media_consultas = mean(consultas),
    media_ps = mean(ps),
    media_exames = mean(exames),
    media_terapias = mean(terapias),
    media_internacoes = mean(internacoes),
    .groups = "drop"
  )

cat("\nTabela com medias das cinco contagens por grupo:\n")
print(
  tabela_contagens %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
)

# ------------------------------------------------------------
# 14) VALIDACAO
# ------------------------------------------------------------

cat("\n--- VALIDACAO: medias-base dos interceptos vs. medias observadas ---\n")
cat("Observacao: as medias observadas incorporam efeitos das covariaveis,\n")
cat("portanto nao precisam coincidir exatamente com os valores-base.\n\n")

alvos <- data.frame(
  grupo_latente = c("baixo_uso","ambulatorial_coordenado","agudo_hospitalar","atipico"),
  consultas_alvo = c(3.0, 8.0, 5.0, 4.0),
  ps_alvo = c(0.4, 0.8, 3.0, 4.0),
  exames_alvo = c(2.0, 10.0, 6.0, 14.0),
  terapias_alvo = c(0.2, 1.5, 0.6, 3.0),
  internacoes_alvo = c(0.05, 0.10, 0.8, 0.5)
)

comparacao_obs <- resumo_grupo %>%
  left_join(alvos, by = "grupo_latente")

print(comparacao_obs)

cat("\n--- VALIDACAO: medias teoricas condicionais (mu) por grupo ---\n")

valid_mu <- data.frame(
  grupo_latente = factor(
    c("baixo_uso","ambulatorial_coordenado","agudo_hospitalar","atipico"),
    levels = c("baixo_uso","ambulatorial_coordenado","agudo_hospitalar","atipico")
  ),
  mu_consultas   = as.numeric(tapply(mu_consultas,   dados_sim$grupo_latente, mean)),
  mu_ps          = as.numeric(tapply(mu_ps,          dados_sim$grupo_latente, mean)),
  mu_exames      = as.numeric(tapply(mu_exames,      dados_sim$grupo_latente, mean)),
  mu_terapias    = as.numeric(tapply(mu_terapias,    dados_sim$grupo_latente, mean)),
  mu_internacoes = as.numeric(tapply(mu_internacoes, dados_sim$grupo_latente, mean))
)

print(valid_mu)

cat("\nProporcao de zeros por variavel:\n")
prop_zeros <- sapply(
  dados_sim[, c("consultas","ps","exames","terapias","internacoes")],
  function(x) mean(x == 0)
)
print(round(prop_zeros, 3))

# --- Diagnostico do gasto -----------------------------------
cat("\nResumo do gasto_total:\n")
print(summary(dados_sim$gasto_total))

cat("\nResumo do log(gasto_total + 1):\n")
print(summary(log(dados_sim$gasto_total + 1)))

# --- Flags por grupo ----------------------------------------
flags_grupo <- dados_sim %>%
  group_by(grupo_latente) %>%
  summarise(
    prop_alta_utilizacao = mean(flag_alta_utilizacao),
    prop_baixa_coerencia = mean(flag_baixa_coerencia),
    .groups = "drop"
  )

cat("\nPrevalencia das flags por grupo:\n")
print(flags_grupo)

# ------------------------------------------------------------
# 15) GRAFICOS
# ------------------------------------------------------------

# 1. Histograma do total de eventos
p1 <- ggplot(dados_sim, aes(x = total_eventos)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(
    title = "Distribuicao do total de eventos",
    x = "Total de eventos",
    y = "Frequencia"
  ) +
  theme_minimal()

print(p1)

# 2. Boxplot total de eventos por grupo
p2 <- ggplot(
  dados_sim,
  aes(x = grupo_latente, y = total_eventos, fill = grupo_latente, group = grupo_latente)
) +
  geom_boxplot() +
  labs(
    title = "Total de eventos por grupo latente",
    x = "Grupo",
    y = "Total de eventos"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

print(p2)

# 3. Boxplot razao exames/consultas por grupo
p3 <- ggplot(
  dados_sim,
  aes(x = grupo_latente, y = razao_exames_consultas, fill = grupo_latente, group = grupo_latente)
) +
  geom_boxplot() +
  labs(
    title = "Razao exames/consultas por grupo latente",
    x = "Grupo",
    y = "Razao exames/consultas"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

print(p3)

# 4. Gasto total x total de eventos por grupo
p4 <- ggplot(
  dados_sim,
  aes(x = total_eventos, y = gasto_total, color = grupo_latente)
) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.9) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Gasto total vs. total de eventos por grupo latente",
    x = "Total de eventos",
    y = "Gasto total (R$, escala log)",
    color = "Grupo"
  ) +
  theme_minimal()

print(p4)

# 5. Distribuicao de meses_utilizacao por grupo
p5 <- ggplot(dados_sim, aes(x = meses_utilizacao, fill = grupo_latente)) +
  geom_histogram(binwidth = 1, color = "white", position = "dodge") +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Distribuicao de meses com utilizacao por grupo",
    x = "Meses com utilizacao",
    y = "Frequencia",
    fill = "Grupo"
  ) +
  theme_minimal()

print(p5)

# 6. Prevalencia das flags por grupo
flags_long <- flags_grupo %>%
  tidyr::pivot_longer(
    cols = c(prop_alta_utilizacao, prop_baixa_coerencia),
    names_to = "tipo_flag",
    values_to = "proporcao"
  )

p6 <- ggplot(flags_long, aes(x = grupo_latente, y = proporcao, fill = tipo_flag)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Prevalencia de sinais de alerta por grupo",
    x = "Grupo latente",
    y = "Proporcao",
    fill = "Indicador"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

print(p6)

# 7. Boxplot razao PS/consultas por grupo
p7 <- ggplot(
  dados_sim,
  aes(x = grupo_latente, y = razao_ps_consultas, fill = grupo_latente, group = grupo_latente)
) +
  geom_boxplot() +
  labs(
    title = "Razao PS/consultas por grupo latente",
    x = "Grupo",
    y = "Razao PS/consultas"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

print(p7)


# ------------------------------------------------------------
# 16) EXPORTAR
# ------------------------------------------------------------

write.csv(dados_sim, "C:/Katy/Doutorado/painel_WCDANM/dados_simulados_saude_suplementar.csv", row.names = FALSE)

cat("\nArquivo salvo: dados_simulados_saude_suplementar.csv\n")