#!/usr/bin/perl
# _inject_sidebar.pl — re-injectável, UTF-8 seguro, suporta CRLF
use strict;
use warnings;
use utf8;          # literais neste script são UTF-8
binmode STDOUT, ':utf8';

# ── CSS override ─────────────────────────────────────────────────────────────
my $CSS_BLOCK = <<'END';
<!-- ss-style-start -->
<style>
/* SimSaúde: empurra conteúdo Bootstrap para dar espaço ao sidebar fixo */
.main-container {
  margin-left: calc(var(--sidebar-w) + 20px) !important;
  max-width: none !important;
}
</style>
<!-- ss-style-end -->
END

# ── Nav HTML ─────────────────────────────────────────────────────────────────
sub nav_html {
    my ($active) = @_;
    my @links = (
        [ "../index.html",    "In\x{ED}cio"      ],
        [ "estudo.html",      "O Estudo"          ],
        [ "evolucao.html",    "Evolu\x{E7}\x{E3}o" ],
        [ "simulacao.html",   "Simula\x{E7}\x{E3}o" ],
        [ "eda_simsaude.html","EDA"                ],
        [ "analises.html",   "An\x{E1}lises"      ],
    );
    my $items = "";
    for my $l (@links) {
        my ($href, $label) = @$l;
        my $cls = ($href eq $active) ? ' class="active"' : "";
        $items .= "    <li><a href=\"$href\"$cls>$label</a></li>\n";
    }
    return <<"END_NAV";
<!-- ss-nav-start -->
<nav class="sidebar">
  <a href="../index.html" class="sidebar-brand">SimSa\x{FA}de</a>
  <p class="sidebar-subtitle">Modelagem de Heterogeneidade<br>em Sa\x{FA}de Suplementar</p>
  <ul class="sidebar-nav">
$items  </ul>
  <div class="sidebar-footer">Katy Garcia de Freitas \x{B7} 2026</div>
</nav>
<!-- ss-nav-end -->
END_NAV
}

# ── Injeção ───────────────────────────────────────────────────────────────────
sub inject {
    my ($filepath, $active) = @_;

    # Lê como bytes brutos para não tocar no encoding existente
    open my $fh, '<:raw', $filepath
        or die "Nao consigo ler $filepath: $!";
    local $/;
    my $raw = <$fh>;
    close $fh;

    # Decodifica para string Unicode
    require Encode;
    my $html = Encode::decode('UTF-8', $raw, Encode::FB_CROAK());

    # ── Remove injeções anteriores (qualquer versão) ──────────────────────────
    # Bloco CSS entre marcadores
    $html =~ s/<!-- ss-style-start -->.*?<!-- ss-style-end -->\R?//s;
    # Nav entre marcadores
    $html =~ s/<!-- ss-nav-start -->.*?<!-- ss-nav-end -->\R?//s;
    # Nav antiga (ss-sidebar, sem marcadores)
    $html =~ s/<nav class="ss-sidebar">.*?<\/nav>//s;
    # Link/style antigo (injeção sem marcadores antes de </head>)
    $html =~ s/<link rel="preconnect" href="https:\/\/fonts\.googleapis\.com">\R.*?<\/style>\R(?=<\/head>)//s;

    # ── Injeta CSS antes de </head> ───────────────────────────────────────────
    $html =~ s{(</head>)}{$CSS_BLOCK$1};

    # ── Injeta nav logo após <body> ───────────────────────────────────────────
    my $nav = nav_html($active);
    $html =~ s{(<body[^>]*>)}{$1\n$nav}s;

    # Escreve de volta como UTF-8 puro
    my $out_bytes = Encode::encode('UTF-8', $html);
    open my $out, '>:raw', $filepath
        or die "Nao consigo escrever $filepath: $!";
    print $out $out_bytes;
    close $out;

    print "  $filepath: OK\n";
}

# ── Executa ───────────────────────────────────────────────────────────────────
my $pages = "docs/pages";
inject("$pages/simulacao.html",    "simulacao.html");
inject("$pages/eda_simsaude.html", "eda_simsaude.html");
print "Pronto.\n";
