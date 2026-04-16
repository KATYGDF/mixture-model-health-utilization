# SimSaúde: Detecção de Padrões Atípicos em Saúde Suplementar via Modelos de Mistura

## Objetivo Científico

**Demonstrar que modelos de mistura finita (finite mixture models) podem identificar e caracterizar padrões de utilização de serviços de saúde que divergem significativamente dos perfis esperados, usando dados sintéticos calibrados a indicadores reais.**

Este projeto é uma **investigação metodológica em ambiente controlado** sobre como estruturas latentes em dados de utilização — particularmente grupos com comportamentos heterogêneos — podem ser recuperadas, caracterizadas e usadas como base para sistemas de detecção de anomalias.

## Contexto e Motivação

Na saúde suplementar brasileira, beneficiários exibem padrões de utilização extremamente heterogêneos:
- Alguns raramente utilizam o plano de saúde
- Outros concentram internações de alto custo
- Alguns apresentam proporções anormais entre exames e consultas

**Este projeto investiga**: 
- Quão bem conseguimos recuperar a estrutura latente (grupos heterogêneos verdadeiros) usando diferentes métodos estatísticos? 
- Qual é a "assinatura" de cada grupo em termos de padrões de utilização e custos? 
- Como algoritmos de detecção de anomalias se comportam quando sabemos a verdade sobre grupos latentes?

## Estrutura do Projeto:

### Fase 1: Modelagem Teórica e Simulação 
**Objetivo**: Criar um gerador de dados sintéticos fundamentado em teoria estatística, validado contra indicadores reais.

**Componentes**:
- **Modelo de mistura finita** com 4 grupos latentes ($G_1, G_2, G_3, G_4$)
- **Distribuições de contagem**: Binomial Negativa (NB₂) para eventos frequentes; Zero-Inflated Negative Binomial (ZINB) para eventos raros (terapias, internações)
- **Regressão log-linear** com coeficientes heterogêneos por grupo
- **Calibração ANS**: Ajuste de interceptos para reproduzir médias populacionais observadas (desvio < 3%)

**Dados simulados**:
- 5.000 beneficiários
- 5 desfechos de contagem: consultas, pronto-socorro, exames, terapias, internações
- Covariáveis: idade, sexo, tempo no plano, coparticipação, região
- Gasto simulado com componente de ruído log-normal

**Validação**:
- ✅ Médias por desfecho vs. Painel de Indicadores ANS
- ✅ Médias teóricas $\bar{\mu}_{gk}$ condicionais e não-enviesadas
- ✅ Estrutura de dispersão realista (θ estimados)
- ✅ Proporções de zeros estruturais condizentes com prática clínica

**Saída**: Base de dados sintética de 5.000 

---

### Fase 2: Recuperação de Estrutura Latente

**Objetivo**: Aplicar **3 algoritmos principais** de modelos de mistura aos dados simulados para demonstrar recuperação de grupos latentes verdadeiros.

**Métodos a implementar** :

#### A. Finite Mixture of Negative Binomial Regressions (FMNB)
- **Pacote**: `flexmix` ou `mixtools` (R)
- **Abordagem**: Estimar modelo de mistura diretamente na escala de contagem com log-link, componentes específicas por grupo
- **Hipótese**: Método mais apropriado para dados de contagem; deve recuperar $\alpha_{gk}$ e $\boldsymbol{\beta}_{gk}$
- **Métrica**: BIC/AIC para seleção de $K$; pureza de classificação posterior; ARI

#### B. Latent Class Analysis (LCA) com Variáveis Derivadas
- **Pacote**: `poLCA` (R)
- **Abordagem**: Discretizar contagens em categorias (baixo/médio/alto uso) e aplicar LCA
- **Hipótese**: LCA deve identificar ~4 classes correspondentes aos grupos verdadeiros
- **Métrica**: Pureza, entropia, índice de Rand ajustado (ARI)

#### C. Model-Based Clustering (Gaussian Mixture Models) em Espaço Transformado
- **Pacote**: `mclust` (R)
- **Abordagem**: Log-transformar contagens + PCA; aplicar GMM com seleção automática de $K$ via BIC
- **Hipótese**: Padrão de G4 (atípico) deve ser claramente separável mesmo em dimensão reduzida
- **Métrica**: Silhueta, pureza, ARI

#### D. Comparação de métodos



---

## Metodologia Estatística Resumida

**Referência Principal**: Ng, Xiang & Yau (2010), *Mixture Modelling for Medical and Health Sciences*, Chapman and Hall/CRC. 

### Modelo Generativo

Para cada beneficiário $i$:

1. **Atribuir grupo latente**: $g_i \sim \text{Categórica}(\pi_1 = 0.5, \pi_2 = 0.25, \pi_3 = 0.20, \pi_4 = 0.05)$
   - *Ref*: Cap. 1 — Fundamental Concepts of Finite Mixture Models

2. **Gerar covariáveis observadas**:
   - $\text{idade}_i \sim \text{Mistura de 3 faixas}$
   - $\text{sexo}_i \sim \text{Bernoulli}(0.55)$
   - $\text{tempo}_{ni}, \text{copart}_i, \text{regiao}_i$ independentes

3. **Gerar contagens de serviços** (Cap. 4 — Mixture of Generalized Linear Models for Count Data):

   **Log-link com intercepto e covariáveis grupo-específicas**: log(μ<sub>igk</sub>) = α<sub>gk</sub> + x<sub>i</sub><sup>T</sup> β<sub>gk</sub>

   **Distribuições por tipo de desfecho**:
   - Para consultas, PS e exames: Y<sub>igk</sub> ~ Poisson(μ<sub>igk</sub>) ou Y<sub>igk</sub> ~ NB₂(μ<sub>igk</sub>, θ<sub>k</sub>)
   - Para terapias e internações: Y<sub>igk</sub> ~ ZINB(μ<sub>igk</sub>, θ<sub>k</sub>, ψ<sub>gk</sub>)

   **Referências**:
   - Cap. 4.2 — *Poisson Mixture Regression Model*
   - Cap. 4.4 — *Zero-Inflated Negative Binomial Regression Models*

### Objetivos de Análise 

Reponder as seguintes perguntas científicas:

- Conseguimos encontrar os 4 grupos verdadeiros?
- Qual método é mais eficiente?
- Qual é a assinatura clínicoeconomica de cada grupo?
- G4 (atípico, 5%) é detectável?" 

---

## Resultados Esperados

### Fase 1 ✅
- **Calibração**: Desvios < 3% vs. ANS em todos os desfechos
- **Estrutura**: 4 grupos com assinaturas clínicas distintas
- **Reprodutibilidade**: Resultados idênticos com mesma semente aleatória

### Fase 2 
- **FMNB**: Recuperação de $\alpha_{gk}$ e $\boldsymbol{\beta}_{gk}$ com pureza > 85%
- **LCA**: Identificação de ~4 classes com entropia baixa (< 0.5)
- **GMM**: Separação clara em espaço transformado, com 4 componentes selecionados via BIC
- **Comparação**: Tabela mostrando qual método tem melhor desempenho por métrica
- **Visualizações**: 
  - Matrix plot: classificação real vs. estimada para cada método
  - Boxplot: distribuição de contagens por grupo inferido
  - Tabela: médias teóricas $\bar{\mu}_{gk}$ vs. estimadas por método

---

## Publicações e Apresentações

### Pôster Científico (WCDANM 2026)
**Título**: Finite Mixture Models for Heterogeneous Healthcare Utilization: Recovering Latent Groups in Synthetic  Health Care Data

**Seções**:
1. Introdução: Heterogeneidade em saúde suplementar
2. Métodos: Simulação + FMNB/LCA/GMM
3. Resultados: Pureza/ARI por método + assinaturas clínicas
4. Discussão: Implicações para detecção de padrões atípicos

---

## Escopo e Limitações

### Escopo Deliberado (Sintético)

Este é um **estudo metodológico em ambiente controlado**. O projeto investiga:
- ✅ Se modelos de mistura conseguem recuperar grupos latentes verdadeiros
- ✅ Qual algoritmo tem melhor desempenho em diferentes cenários
- ✅ Como padrões atípicos (G4) aparecem em métodos de detecção de anomalias
- ✅ Assinaturas clínicoeconomicas de cada grupo latente

Não se propõe a:
- ❌ Validar em dados reais de operadoras
- ❌ Demonstrar impacto operacional em auditoria
- ❌ Quantificar taxa de detecção de fraude em produção

### Limitações

1. **Dados sintéticos**: Calibrados a indicadores agregados; não capturam desvios idiossincráticos reais
2. **Coeficientes não-estimados**: $\boldsymbol{\beta}_{gk}$ baseados em julgamento especializado, não em ML
3. **Ausência de temporalidade**: Modelo transversal; não há dinâmica de transição entre grupos
4. **Covariáveis limitadas**: Apenas idade, sexo, tempo, copart, região; diagnósticos não inclusos
5. **Validação por simulação**: Grupos latentes são **conhecidos por construção** (não cegos)



---

## Autor

**Katy Garcia de Freitas**  
Doutoranda em Matemática (Universidade do Minho)  
Data Professional & Healthcare Analytics Consultant  


**Contato**: katygdf@gmail.com  
**GitHub**: [@KATYGDF](https://github.com/KATYGDF)

---

## Referências Principais

- **[LIVRO PRINCIPAL]** Ng, S. K. (2019). Mixture modelling for medical and health sciences. Chapman and Hall/CRC.
  
