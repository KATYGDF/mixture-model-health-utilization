# SimSaúde · Documentação do Projeto de Simulação

Site de documentação técnica para o estudo de simulação de dados de utilização em saúde suplementar.

## Estrutura

```
mixture-model-health-utilization/
├── index.html                  # Página inicial
├── css/
│   └── main.css                # Estilos globais
└── pages/
    ├── metodologia.html        # Fundamentos matemáticos e fórmulas
    ├── documentacao.html       # Passo a passo do script R
    ├── evolucao.html           # Histórico de versões v1→v3
    └── codigo.html             # Código-fonte com syntax highlighting
```

## Como publicar no GitHub Pages

1. Crie um repositório no GitHub (ex: `simsaude-docs`)
2. Faça upload de todos os arquivos mantendo a estrutura de pastas
3. Vá em **Settings → Pages**
4. Em **Source**, selecione `Deploy from a branch`
5. Escolha a branch `main` e pasta `/ (root)`
6. Clique em **Save**

O site ficará disponível em `https://SEU_USUARIO.github.io/simsaude-docs/`

## Tecnologias

- HTML5 + CSS3 puro (sem framework)
- [MathJax 3](https://www.mathjax.org/) — renderização de fórmulas LaTeX
- [highlight.js](https://highlightjs.org/) — syntax highlighting do código R
- Google Fonts: Playfair Display, Source Serif 4, IBM Plex Mono

## Referências

- Painel de Indicadores da ANS — indicadores de utilização
- Cameron & Trivedi (2013) — *Regression Analysis of Count Data*
- Lambert (1992) — Zero-Inflated Poisson Regression
