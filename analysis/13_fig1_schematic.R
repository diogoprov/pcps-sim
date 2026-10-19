# Figure 1: schematic of the simulation design (tree shape -> traits ->
# assembly -> analysis -> inference). Hand-built SVG; the trees of panel 1 come
# from the tree-shape grid (analysis/06_tree_grid.R) and the matrices of panel 4
# from one simulated metacommunity. Rendered to PDF and PNG with rsvg-convert.
# Usage: Rscript analysis/13_fig1_schematic.R

for (f in list.files("R", full.names = TRUE)) source(f)
out_dir <- file.path("manuscript", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Style -------------------------------------------------------------------------
pal_env <- c("#26456e", "#1c7c9c", "#3fae9a", "#a6cf6b", "#f2c14e", "#e8833a")
col_phylo <- "#6c4f9e"
col_h1 <- "#d1495b"
col_h0 <- "#8d99ae"
col_text <- "#22303c"
col_card <- "#f5f6f8"
col_line <- "#3d4a5c"
font <- "Helvetica Neue, Helvetica, Arial, sans-serif"

env_col <- function(v) {
  m <- grDevices::colorRamp(pal_env)(pmin(pmax(v, 0), 1))
  grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
}
mix <- function(col, with = "#ffffff", w = 0.5) {
  a <- grDevices::col2rgb(col)
  b <- grDevices::col2rgb(with)
  m <- (1 - w) * a + w * b
  grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}
ramp <- function(v, col) {
  m <- grDevices::colorRamp(c("#ffffff", col))(pmin(pmax(v, 0), 1))
  grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
}

# SVG primitives ----------------------------------------------------------------
svg <- character(0)
add <- function(...) svg <<- c(svg, paste0(...))
num <- function(x) formatC(x, format = "f", digits = 1)
rect <- function(x, y, w, h, fill = "none", stroke = "none", sw = 1, rx = 0,
                 opacity = 1) {
  add(sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s" stroke="%s" stroke-width="%s" opacity="%s"/>',
              num(x), num(y), num(w), num(h), rx, fill, stroke, sw, opacity))
}
circle <- function(cx, cy, r, fill = "none", stroke = "none", sw = 1, opacity = 1) {
  add(sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s" stroke="%s" stroke-width="%s" opacity="%s"/>',
              num(cx), num(cy), num(r), fill, stroke, sw, opacity))
}
line <- function(x1, y1, x2, y2, stroke = col_line, sw = 1.5, dash = NULL,
                 marker = NULL, both = FALSE, cap = "round") {
  add(sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="%s" stroke-linecap="%s"%s%s%s/>',
              num(x1), num(y1), num(x2), num(y2), stroke, sw, cap,
              if (is.null(dash)) "" else sprintf(' stroke-dasharray="%s"', dash),
              if (is.null(marker)) "" else sprintf(' marker-end="url(#%s)"', marker),
              if (both) sprintf(' marker-start="url(#%s)"', marker) else ""))
}
path <- function(d, fill = "none", stroke = col_line, sw = 1.5, opacity = 1,
                 marker = NULL, join = "round") {
  add(sprintf('<path d="%s" fill="%s" stroke="%s" stroke-width="%s" opacity="%s" stroke-linejoin="%s" stroke-linecap="round"%s/>',
              d, fill, stroke, sw, opacity, join,
              if (is.null(marker)) "" else sprintf(' marker-end="url(#%s)"', marker)))
}
text <- function(x, y, label, size = 14, weight = "normal", anchor = "middle",
                 fill = col_text, style = "normal") {
  add(sprintf('<text x="%s" y="%s" font-family="%s" font-size="%s" font-weight="%s" font-style="%s" text-anchor="%s" fill="%s">%s</text>',
              num(x), num(y), font, size, weight, style, anchor, fill, label))
}
poly <- function(x, y) paste0("M", paste(num(x), num(y), collapse = " L"))

# Tree drawn upwards in the box (x, y, w, h); returns tip coordinates
tree_xy <- function(tree) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  ape::plot.phylo(tree, plot = FALSE)
  lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  list(depth = lp$xx / max(lp$xx), pos = (lp$yy - 1) / (ape::Ntip(tree) - 1))
}
draw_tree <- function(tree, x, y, w, h, col = col_line, sw = 1.2, stem = 0.06) {
  xy <- tree_xy(tree)
  X <- x + xy$pos * w
  Y <- y + h * (1 - stem) * (1 - xy$depth)
  e <- tree$edge
  d <- c(sprintf("M%s %s L%s %s", num(X[e[, 2]]), num(Y[e[, 1]]),
                 num(X[e[, 2]]), num(Y[e[, 2]])))
  for (p in unique(e[, 1])) {
    k <- e[e[, 1] == p, 2]
    d <- c(d, sprintf("M%s %s L%s %s", num(min(X[k])), num(Y[p]), num(max(X[k])), num(Y[p])))
  }
  root <- ape::Ntip(tree) + 1
  d <- c(d, sprintf("M%s %s L%s %s", num(X[root]), num(Y[root]), num(X[root]), num(y + h)))
  path(paste(d, collapse = " "), stroke = col, sw = sw)
  n <- ape::Ntip(tree)
  invisible(list(x = X[seq_len(n)], y = Y[seq_len(n)]))
}
heat <- function(mat, x, y, w, h, col_fun) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  cw <- w / nc
  ch <- h / nr
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    if (is.na(mat[i, j])) next
    rect(x + (j - 1) * cw, y + (i - 1) * ch, cw + 0.4, ch + 0.4, fill = col_fun(mat[i, j]))
  }
  rect(x, y, w, h, stroke = "#c5ccd6", sw = 1)
}
bell <- function(x0, base, width, mu, sd, height) {
  t <- seq(0, 1, length.out = 80)
  poly(x0 + t * width, base - height * exp(-(t - mu)^2 / (2 * sd^2)))
}
bell_icon <- function(x0, base, width, mu, sd, height, col, fill_op = 0.25) {
  d <- bell(x0, base, width, mu, sd, height)
  path(sprintf("%s L%s %s L%s %s Z", d, num(x0 + width), num(base), num(x0), num(base)),
       fill = col, stroke = "none", opacity = fill_op)
  path(d, stroke = col, sw = 2)
}
badge <- function(x, y, label, fill = col_text, r = 16, size = 18) {
  circle(x, y, r, fill = fill)
  text(x, y + size * 0.36, label, size = size, weight = "bold", fill = "#ffffff")
}
card <- function(x, y, w, h, number, title) {
  rect(x, y, w, h, fill = col_card, stroke = "#e3e7ec", rx = 14)
  badge(x + 28, y + 30, number)
  text(x + 54, y + 38, title, size = 22, weight = "bold", anchor = "start")
}
# mini tree glyph (a "Y") used for phylogenetic signal
glyph <- function(cx, cy, s = 6, col = "#ffffff") {
  path(sprintf("M%s %s L%s %s M%s %s L%s %s M%s %s L%s %s",
               num(cx), num(cy + s), num(cx), num(cy),
               num(cx), num(cy), num(cx - s * 0.8), num(cy - s * 0.8),
               num(cx), num(cy), num(cx + s * 0.8), num(cy - s * 0.8)),
       stroke = col, sw = 1.6)
}
signal_dot <- function(cx, cy, level, r = 11) {
  if (level == 0) {
    circle(cx, cy, r, fill = "#ffffff", stroke = "#b8c0cc", sw = 1.5)
  } else if (level < 1) {
    circle(cx, cy, r, fill = "#ffffff", stroke = col_phylo, sw = 1.5)
    path(sprintf("M%s %s A%s %s 0 0 1 %s %s Z", num(cx), num(cy - r), r, r, num(cx), num(cy + r)),
         fill = col_phylo, stroke = "none")
  } else {
    circle(cx, cy, r, fill = col_phylo)
    glyph(cx, cy + 1, s = 5)
  }
}
chevron <- function(x, y) path(sprintf("M%s %s L%s %s L%s %s", num(x - 5), num(y - 12),
                                       num(x + 6), num(y), num(x - 5), num(y + 12)),
                               stroke = "#9aa5b4", sw = 3.5)
arrow_curve <- function(x1, y1, x2, y2, bend = 0.25, sw = 1.6, col = col_line,
                        marker = "arr") {
  mx <- (x1 + x2) / 2 - bend * (y2 - y1)
  my <- (y1 + y2) / 2 + bend * (x2 - x1)
  path(sprintf("M%s %s Q%s %s %s %s", num(x1), num(y1), num(mx), num(my), num(x2), num(y2)),
       stroke = col, sw = sw, marker = marker)
}

# Data ------------------------------------------------------------------------------
grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds"))
grid50 <- grid[grid$n == 50, ]
pick_tree <- function(i) {
  s <- grid50$shape[[i]]
  sc <- (s$gamma - grid50$gamma_target[i])^2 / stats::var(s$gamma) +
    (s$colless - grid50$colless_target[i])^2 / stats::var(s$colless)
  grid50$trees[[i]][[which.min(sc)]]
}
grid_trees <- stats::setNames(lapply(seq_len(nrow(grid50)), pick_tree), grid50$cell)

set.seed(20261018)
mc <- sim_metacommunity_A(sim_tree_shape(16), lambda = 1, sigma_signal = FALSE,
                          h_signal = FALSE, env_original = TRUE, incidence = FALSE,
                          n_sites = 20, mu_tol = 14)
tr <- mc$tree
tip_order <- order(tree_xy(tr)$pos[seq_len(ape::Ntip(tr))])  # left-to-right tips
u <- stats::setNames(mc$traits$u, tr$tip.label)
site_order <- order(mc$env)
comm <- mc$comm[site_order, tr$tip.label[tip_order]]
env <- mc$env[site_order]
q <- phylo_q(tr)[tr$tip.label[tip_order], tr$tip.label[tip_order]]
p <- matrix_p(comm, q)
dp <- d_p(p)
scores <- pcps_axes(comm, q, 1:2)
env01 <- (env - min(env)) / diff(range(env))

# Canvas ------------------------------------------------------------------------
W <- 1800
H <- 820
add(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">', W, H, W, H))
add('<defs>')
add('<linearGradient id="envh" x1="0" y1="0" x2="1" y2="0">',
    paste(sprintf('<stop offset="%s" stop-color="%s"/>', seq(0, 1, length.out = 6), pal_env), collapse = ""),
    '</linearGradient>')
add('<linearGradient id="envv" x1="0" y1="0" x2="0" y2="1">',
    paste(sprintf('<stop offset="%s" stop-color="%s"/>', seq(0, 1, length.out = 6), pal_env), collapse = ""),
    '</linearGradient>')
for (m in list(c("arr", col_line), c("arrg", "#9aa5b4"), c("arrp", col_phylo))) {
  add(sprintf('<marker id="%s" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="%s"/></marker>',
              m[1], m[2]))
}
add('</defs>')
rect(0, 0, W, H, fill = "#ffffff")

cards <- data.frame(x = c(20, 355, 740, 1135, 1570), w = c(310, 360, 370, 410, 210))
titles <- c("Tree shape", "Traits", "Assembly", "Analysis", "Inference")
for (i in seq_len(nrow(cards))) card(cards$x[i], 20, cards$w[i], 780, i, titles[i])
for (i in 1:4) chevron((cards$x[i] + cards$w[i] + cards$x[i + 1]) / 2, 410)

# 1 Tree shape ------------------------------------------------------------------
cw <- 78
chh <- 150
gx <- 74
gy <- 84
rows <- c("yule", "mid", "low")    # top = Yule gamma (0), bottom = low gamma
cols <- c("yule", "mid", "high")   # left = Yule Colless, right = high
for (i in 1:3) for (j in 1:3) {
  x0 <- gx + (j - 1) * (cw + 9)
  y0 <- gy + (i - 1) * (chh + 10)
  is_yule <- i == 1 && j == 1
  rect(x0, y0, cw, chh, fill = "#ffffff", stroke = if (is_yule) col_phylo else "#e1e5ea",
       sw = if (is_yule) 2 else 1, rx = 8)
  draw_tree(grid_trees[[paste0("g_", rows[i], "__c_", cols[j])]], x0 + 6, y0 + 8,
            cw - 12, chh - 14, sw = 0.6)
}
gx2 <- gx + 3 * cw + 18
gy2 <- gy + 3 * chh + 20
line(48, gy2, 48, gy + 4, sw = 2.2, marker = "arr")
text(36, (gy + gy2) / 2 + 8, "γ", size = 26, style = "italic")
line(gx, gy2 + 22, gx2, gy2 + 22, sw = 2.2, marker = "arr")
text((gx + gx2) / 2, gy2 + 52, 'I<tspan font-size="15" dy="4">C</tspan>', size = 24,
     style = "italic")
text(gx + cw / 2, gy - 8, "Yule", size = 13, fill = col_phylo, weight = "bold")
# pool size
ny <- 715
text(48, ny + 8, "N", size = 24, style = "italic", weight = "bold")
for (k in 1:3) {
  cx <- 110 + (k - 1) * 82
  nt <- c(5, 8, 12)[k]
  set.seed(k)
  draw_tree(ape::rcoal(nt), cx - 25, ny - 40, 50, 40,
            sw = 1.4)
  for (t in seq_len(nt)) circle(cx - 25 + (t - 1) * 50 / (nt - 1), ny - 40, 2.6,
                                fill = col_phylo)
  text(cx, ny + 24, c("50", "100", "200")[k], size = 16, weight = "bold")
}

# 2 Traits ----------------------------------------------------------------------
tx <- 395
tw <- 290
tips <- draw_tree(tr, tx, 92, tw, 96, sw = 1.4)
ux <- (u - min(u)) / diff(range(u))
for (t in seq_along(tips$x)) circle(tips$x[t], tips$y[t] - 0.5, 6.5, fill = env_col(ux[t]),
                                    stroke = "#ffffff", sw = 1.5)
icons <- list(
  list(lab = "u", word = "optimum",
       curves = list(c(0.22, 0.10, 46, 0.1), c(0.5, 0.10, 46, 0.5), c(0.78, 0.10, 46, 0.9))),
  list(lab = "σ", word = "breadth",
       curves = list(c(0.5, 0.07, 46, 0.25), c(0.5, 0.2, 46, 0.75))),
  list(lab = "h", word = "abundance",
       curves = list(c(0.5, 0.11, 50, 0.85), c(0.5, 0.11, 22, 0.3)))
)
for (k in seq_along(icons)) {
  x0 <- 372 + (k - 1) * 112
  rect(x0, 222, 100, 76, fill = "#ffffff", stroke = "#e1e5ea", rx = 8)
  for (cv in icons[[k]]$curves) bell_icon(x0 + 8, 288, 84, cv[1], cv[2], cv[3], env_col(cv[4]))
  line(x0 + 8, 288, x0 + 92, 288, stroke = "#b8c0cc", sw = 1.2)
  text(x0 + 30, 324, icons[[k]]$lab, size = 20, style = "italic", weight = "bold")
  text(x0 + 44, 323, icons[[k]]$word, size = 12, fill = "#5c6878", anchor = "start")
}
# scenario matrix
sc <- scenarios_layerA()
lab <- c(S0 = "S0", S1 = "S1", S2 = "S2", S2sigma = "S2σ", S2h = "S2h",
         S3a = "S3a", S3b = "S3b", S3sigma = "S3σ")
colx <- c(470, 528, 586, 646)
pillx <- 690
hy <- 362
for (k in 1:3) text(colx[k], hy, c("u", "σ", "h")[k], size = 18, style = "italic", weight = "bold")
text(colx[4], hy, "E", size = 18, style = "italic", weight = "bold")
line(372, hy + 12, 700, hy + 12, stroke = "#d5dbe3", sw = 1)
for (i in seq_len(nrow(sc))) {
  y <- 396 + (i - 1) * 41
  if (i %in% c(3, 6)) line(372, y - 22, 700, y - 22, stroke = "#d5dbe3", sw = 1, dash = "3 4")
  text(372, y + 5, lab[[sc$scenario[i]]], size = 15, weight = "bold", anchor = "start")
  signal_dot(colx[1], y, sc$lambda[i])
  signal_dot(colx[2], y, as.numeric(sc$sigma_signal[i]))
  signal_dot(colx[3], y, as.numeric(sc$h_signal[i]))
  if (sc$env_original[i]) {
    rect(colx[4] - 18, y - 7, 36, 14, fill = "url(#envh)", rx = 4)
  } else {
    rect(colx[4] - 18, y - 7, 36, 14, fill = "#e4e8ee", stroke = "#b8c0cc", rx = 4)
    path(sprintf("M%s %s L%s %s M%s %s L%s %s", num(colx[4] - 8), num(y - 4), num(colx[4] + 8),
                 num(y + 4), num(colx[4] - 8), num(y + 4), num(colx[4] + 8), num(y - 4)),
         stroke = "#8d99ae", sw = 1.5)
  }
  h1 <- sc$estimates[i] == "power"
  rect(pillx - 15, y - 11, 30, 22, fill = if (h1) col_h1 else col_h0, rx = 11)
  text(pillx, y + 5, if (h1) "H₁" else "H₀", size = 13, weight = "bold", fill = "#ffffff")
}
# legend of the dots
ly <- 752
signal_dot(388, ly, 1, r = 8)
text(400, ly + 4, "phylogenetic signal", size = 12, anchor = "start", fill = "#5c6878")
signal_dot(526, ly, 0.5, r = 8)
text(538, ly + 4, "λ = 0.5", size = 12, anchor = "start", fill = "#5c6878")
signal_dot(606, ly, 0, r = 8)
text(618, ly + 4, "none", size = 12, anchor = "start", fill = "#5c6878")
rect(380, ly + 20, 26, 12, fill = "url(#envh)", rx = 3)
text(412, ly + 30, "generating E", size = 12, anchor = "start", fill = "#5c6878")
rect(518, ly + 20, 26, 12, fill = "#e4e8ee", stroke = "#b8c0cc", rx = 3)
text(550, ly + 30, "independent E", size = 12, anchor = "start", fill = "#5c6878")

# 3 Assembly --------------------------------------------------------------------
ax <- 760
badge(ax + 14, 88, "A", fill = col_line, r = 13, size = 15)
text(ax + 34, 93, "static", size = 13, anchor = "start", fill = "#5c6878")
bx <- 790
bw <- 292
base <- 270
sp <- list(c(0.12, 0.07, 70), c(0.33, 0.11, 110), c(0.55, 0.06, 92), c(0.74, 0.13, 60),
           c(0.9, 0.08, 120))
for (s in sp) bell_icon(bx, base, bw, s[1], s[2], s[3], env_col(s[1]), fill_op = 0.22)
rect(bx, base + 10, bw, 16, fill = "url(#envh)", rx = 8)
text(bx - 14, base + 24, "E", size = 18, style = "italic", weight = "bold")
for (k in 0:7) circle(bx + 18 + k * (bw - 36) / 7, base + 18, 4.5, fill = "#ffffff",
                      stroke = col_text, sw = 1.4)
# Poisson die
dx <- 1040
dy <- 92
rect(dx, dy, 36, 36, fill = "#ffffff", stroke = col_line, sw = 1.8, rx = 7)
for (pp in list(c(10, 10), c(26, 10), c(18, 18), c(10, 26), c(26, 26)))
  circle(dx + pp[1], dy + pp[2], 2.8, fill = col_line)
line(760, 330, 1090, 330, stroke = "#d5dbe3", sw = 1)

badge(ax + 14, 360, "B", fill = col_line, r = 13, size = 15)
text(ax + 34, 365, "dynamic", size = 13, anchor = "start", fill = "#5c6878")
patch <- data.frame(x = c(790, 908, 790, 908), y = c(392, 392, 510, 510),
                    e = c(0.1, 0.45, 0.65, 0.95))
set.seed(3)
for (k in seq_len(nrow(patch))) {
  rect(patch$x[k], patch$y[k], 80, 80, fill = mix(env_col(patch$e[k]), w = 0.72),
       stroke = env_col(patch$e[k]), sw = 1.5, rx = 12)
  sps <- c(rep(patch$e[k], 4), stats::runif(2))
  for (d in seq_along(sps)) circle(patch$x[k] + stats::runif(1, 16, 64),
                                   patch$y[k] + stats::runif(1, 16, 64), 5.5,
                                   fill = env_col(sps[d]), stroke = "#ffffff", sw = 1.2)
}
line(876, 432, 902, 432, sw = 2.2, marker = "arr", both = TRUE)
line(876, 550, 902, 550, sw = 2.2, marker = "arr", both = TRUE)
line(830, 478, 830, 504, sw = 2.2, marker = "arr", both = TRUE)
line(948, 478, 948, 504, sw = 2.2, marker = "arr", both = TRUE)
line(874, 476, 904, 506, sw = 2.2, marker = "arr", both = TRUE)
# loop: 1000 steps
lx <- 1046
ly2 <- 432
path(sprintf("M%s %s A26 26 0 1 1 %s %s", num(lx + 18), num(ly2 + 19), num(lx + 26), num(ly2 - 2)),
     stroke = col_line, sw = 2.4, marker = "arr")
text(lx, ly2 + 6, "1000", size = 14, weight = "bold")
text(lx, ly2 + 48, "time steps", size = 12, fill = "#5c6878")
# competition
cx3 <- 1046
cy3 <- 548
path(sprintf("M%s %s L%s %s L%s %s M%s %s L%s %s", num(cx3 - 18), num(cy3 - 14), num(cx3 - 18),
             num(cy3 - 30), num(cx3 + 18), num(cy3 - 30), num(cx3 + 18), num(cy3 - 30),
             num(cx3 + 18), num(cy3 - 14)), stroke = col_phylo, sw = 2)
line(cx3, cy3 - 30, cx3, cy3 - 40, stroke = col_phylo, sw = 2)
circle(cx3 - 18, cy3 - 4, 8, fill = env_col(0.3), stroke = "#ffffff", sw = 1.2)
circle(cx3 + 18, cy3 - 4, 8, fill = env_col(0.4), stroke = "#ffffff", sw = 1.2)
line(cx3 - 6, cy3 - 4, cx3 + 6, cy3 - 4, stroke = col_h1, sw = 2, marker = "arr", both = TRUE)
text(cx3, cy3 + 28, "α", size = 20, style = "italic", weight = "bold")
# archetypes
yy <- 640
line(760, yy - 22, 1090, yy - 22, stroke = "#d5dbe3", sw = 1, dash = "3 4")
text(775, yy + 26, "a", size = 18, style = "italic", weight = "bold")
line(808, yy + 20, 858, yy + 20, sw = 1.4, marker = "arr")
text(833, yy + 50, "sorting", size = 12, fill = "#5c6878")
for (o in c(-9, 0, 9)) line(908, yy + 20 + o, 958, yy + 20 + o, sw = 3.2, marker = "arr")
text(933, yy + 50, "mass effects", size = 12, fill = "#5c6878")
text(775, yy + 108, "α", size = 18, style = "italic", weight = "bold")
circle(818, yy + 100, 7, fill = env_col(0.3), stroke = "#ffffff")
circle(848, yy + 100, 7, fill = env_col(0.7), stroke = "#ffffff")
text(833, yy + 132, "none", size = 12, fill = "#5c6878")
path(sprintf("M%s %s L%s %s L%s %s L%s %s", 918, yy + 92, 918, yy + 80, 948, yy + 80, 948, yy + 92),
     stroke = col_phylo, sw = 2)
circle(918, yy + 100, 7, fill = env_col(0.3), stroke = "#ffffff")
circle(948, yy + 100, 7, fill = env_col(0.35), stroke = "#ffffff")
text(933, yy + 132, "phylogenetic", size = 12, fill = "#5c6878")
text(1012, yy + 66, "×", size = 30, fill = "#9aa5b4")
text(1050, yy + 60, "4", size = 26, weight = "bold")
text(1050, yy + 80, "archetypes", size = 11, fill = "#5c6878")

# 4 Analysis --------------------------------------------------------------------
x4 <- 1160
hy4 <- 118
hw <- 92
hh <- 118
rect(x4 - 10, hy4, 6, hh, fill = "url(#envv)", rx = 2)
heat(log1p(comm) / max(log1p(comm)), x4, hy4, hw, hh, \(v) ramp(v, "#1c7c9c"))
text(x4 + hw / 2, hy4 + hh + 28, "W", size = 20, style = "italic", weight = "bold")
text(x4 + hw + 20, hy4 + hh / 2 + 9, "×", size = 26, fill = "#5c6878")
qx <- x4 + hw + 40
draw_tree(tr, qx, hy4 - 2, hw, 26, col = col_phylo, sw = 1)
heat(q / max(q), qx, hy4 + 30, hw, hw, \(v) ramp(v^0.7, col_phylo))
text(qx + hw / 2, hy4 + hh + 28, "Q", size = 20, style = "italic", weight = "bold")
text(qx + hw + 20, hy4 + hh / 2 + 9, "=", size = 26, fill = "#5c6878")
px <- qx + hw + 40
heat(p / max(p), px, hy4, hw, hh, \(v) ramp(v^0.8, "#26456e"))
text(px + hw / 2, hy4 + hh + 28, "P", size = 20, style = "italic", weight = "bold")
# D_P ~ E
line(px + hw / 2, hy4 + hh + 40, px + hw / 2, 316, sw = 2, marker = "arr")
dy4 <- 322
ds <- 92
dlow <- dp
dlow[upper.tri(dlow, diag = TRUE)] <- NA
heat(dlow / max(dp), px, dy4, ds, ds, \(v) ramp(v, "#e8833a"))
text(px - 14, dy4 + ds + 4, 'D<tspan font-size="13" dy="4">P</tspan>', size = 20,
     style = "italic", weight = "bold", anchor = "end")
text(px - 40, dy4 + 56, "~", size = 30, fill = "#5c6878")
rect(px - 80, dy4, 12, ds, fill = "url(#envv)", rx = 3)
text(px - 74, dy4 + ds + 22, "E", size = 18, style = "italic", weight = "bold")
# tests
ty <- 470
tsz <- 62
tests <- c("PERMANOVA", "Mantel", "PCPS–GLM", "RDA", "PCC")
tcx <- x4 - 2 + (0:4) * 81 + tsz / 2
line(1150, ty - 18, 1530, ty - 18, stroke = "#d5dbe3", sw = 1)
for (k in 1:5) {
  rect(tcx[k] - tsz / 2, ty, tsz, tsz, fill = "#ffffff", stroke = "#e1e5ea", rx = 8)
  text(tcx[k], ty + tsz + 20, tests[k], size = 11.5, weight = "bold")
}
# PERMANOVA: three groups
set.seed(5)
for (g in 1:3) {
  gc <- list(c(-14, 10), c(12, 12), c(0, -13))[[g]]
  circle(tcx[1] + gc[1], ty + 31 + gc[2], 12, fill = mix(env_col(c(0.1, 0.55, 0.95)[g]), w = 0.7))
  for (k in 1:4) circle(tcx[1] + gc[1] + stats::rnorm(1, 0, 4), ty + 31 + gc[2] + stats::rnorm(1, 0, 4),
                        2.4, fill = env_col(c(0.1, 0.55, 0.95)[g]))
}
# Mantel: distance vs distance
mx <- tcx[2] - 24
my <- ty + 52
line(mx, my, mx + 48, my, stroke = "#9aa5b4", sw = 1.2)
line(mx, my, mx, my - 42, stroke = "#9aa5b4", sw = 1.2)
dd <- as.vector(stats::dist(env))[1:40]
dpv <- as.matrix(dp)[lower.tri(dp)][1:40]
for (k in seq_along(dd)) circle(mx + 3 + 42 * dd[k] / max(dd), my - 3 - 36 * dpv[k] / max(dpv), 1.8,
                                fill = "#e8833a", opacity = 0.8)
line(mx + 3, my - 6, mx + 46, my - 38, stroke = col_text, sw = 1.6)
# PCPS-GLM: PCPS scores coloured by E
sx <- scores[, 1]
sy <- scores[, 2]
for (k in seq_along(sx)) circle(tcx[3] - 24 + 48 * (sx[k] - min(sx)) / diff(range(sx)),
                                ty + 54 - 46 * (sy[k] - min(sy)) / diff(range(sy)), 3.2,
                                fill = env_col(env01[k]), stroke = "#ffffff", sw = 0.6)
# RDA: triplot (sites, species, E arrow)
line(tcx[4] - 26, ty + 31, tcx[4] + 26, ty + 31, stroke = "#e1e5ea", sw = 0.8)
line(tcx[4], ty + 6, tcx[4], ty + 56, stroke = "#e1e5ea", sw = 0.8)
for (k in seq_along(sx)) circle(tcx[4] - 24 + 48 * (sx[k] - min(sx)) / diff(range(sx)),
                                ty + 54 - 46 * (sy[k] - min(sy)) / diff(range(sy)), 2,
                                fill = "#b8c0cc")
# species: crosses, as in vegan's triplot
spx <- c(-18, -10, -3, 6, 12, 19, -14, 9)
spy <- c(8, -12, 14, -6, 10, -14, -2, 3)
for (k in seq_along(spx)) {
  line(tcx[4] + spx[k] - 3, ty + 31 + spy[k], tcx[4] + spx[k] + 3, ty + 31 + spy[k],
       stroke = col_phylo, sw = 1.3)
  line(tcx[4] + spx[k], ty + 28 + spy[k], tcx[4] + spx[k], ty + 34 + spy[k],
       stroke = col_phylo, sw = 1.3)
}
line(tcx[4], ty + 31, tcx[4] + 22, ty + 13, stroke = col_h1, sw = 2.4, marker = "arr")
text(tcx[4] + 22, ty + 12, "E", size = 11, style = "italic", weight = "bold", anchor = "start",
     fill = col_h1)
# PCC: tree + eigenvector bars
draw_tree({set.seed(2); ape::rcoal(6)}, tcx[5] - 24, ty + 8, 48, 22,
          col = col_phylo, sw = 1.2)
bars <- c(0.8, 0.6, 0.2, -0.3, -0.6, -0.8)
for (k in 1:6) {
  bxk <- tcx[5] - 24 + (k - 1) * 48 / 5
  rect(bxk - 3, ty + 44 - pmax(bars[k], 0) * 14, 6, abs(bars[k]) * 14,
       fill = if (bars[k] > 0) "#1c7c9c" else "#e8833a")
}
line(tcx[5] - 28, ty + 44, tcx[5] + 28, ty + 44, stroke = "#9aa5b4", sw = 1)
# null models
ny4 <- 610
line(1150, ny4 - 14, 1530, ny4 - 14, stroke = "#d5dbe3", sw = 1)
nmx <- c(1200, 1335, 1470)
for (k in 1:3) rect(nmx[k] - 50, ny4, 100, 92, fill = "#ffffff", stroke = "#e1e5ea", rx = 8)
# site shuffle: rows permuted
for (r in 1:5) rect(nmx[1] - 14, ny4 + 12 + (r - 1) * 14, 44, 10,
                    fill = env_col(c(0.1, 0.7, 0.3, 0.9, 0.5)[r]), rx = 2)
path(sprintf("M%s %s C%s %s %s %s %s %s", nmx[1] - 22, ny4 + 16, nmx[1] - 42, ny4 + 30,
             nmx[1] - 42, ny4 + 60, nmx[1] - 22, ny4 + 72), stroke = col_line, sw = 1.8,
     marker = "arr")
text(nmx[1], ny4 + 118, "site shuffle", size = 12, weight = "bold")
# taxa shuffle: tips swapped
set.seed(4)
st <- ape::rcoal(5)
tp <- draw_tree(st, nmx[2] - 36, ny4 + 10, 72, 44, col = col_phylo, sw = 1.4)
tipc <- c(0.1, 0.3, 0.55, 0.75, 0.95)
for (t in 1:5) circle(tp$x[t], tp$y[t] + 2, 5, fill = env_col(tipc[t]), stroke = "#ffffff")
arrow_curve(tp$x[1], ny4 + 66, tp$x[4], ny4 + 66, bend = 0.25, sw = 1.6)
arrow_curve(tp$x[4], ny4 + 66, tp$x[1], ny4 + 66, bend = -0.25, sw = 1.6)
text(nmx[2], ny4 + 118, "taxa shuffle", size = 12, weight = "bold")
# max(p)
rect(nmx[3] - 26, ny4 + 54, 20, 18, fill = mix(col_h0, w = 0.4), rx = 2)
rect(nmx[3] + 6, ny4 + 30, 20, 42, fill = col_h0, rx = 2)
path(sprintf("M%s %s L%s %s", nmx[3] - 34, ny4 + 72, nmx[3] + 34, ny4 + 72), stroke = "#9aa5b4", sw = 1)
line(nmx[3] + 16, ny4 + 8, nmx[3] + 16, ny4 + 24, sw = 2, marker = "arr")
text(nmx[3] - 16, ny4 + 86, 'p<tspan font-size="9" dy="3">site</tspan>', size = 12, style = "italic")
text(nmx[3] + 16, ny4 + 86, 'p<tspan font-size="9" dy="3">taxa</tspan>', size = 12, style = "italic")
text(nmx[3], ny4 + 118, "max(p)", size = 12, weight = "bold", style = "italic")
text(1340, 792, "999 permutations", size = 12, fill = "#5c6878")

# 5 Inference -------------------------------------------------------------------
x5 <- 1595
# type I
rect(x5, 86, 160, 160, fill = "#ffffff", stroke = "#e1e5ea", rx = 8)
line(x5 + 18, 222, x5 + 146, 222, stroke = "#9aa5b4", sw = 1.2)
bh <- c(46, 50, 44, 98, 49)
for (k in 1:5) rect(x5 + 26 + (k - 1) * 24, 222 - bh[k], 16, bh[k],
                    fill = if (bh[k] > 70) col_h1 else col_h0, rx = 2)
line(x5 + 18, 172, x5 + 146, 172, stroke = col_text, sw = 1.6, dash = "5 4")
text(x5 + 150, 168, "α", size = 15, style = "italic", weight = "bold", anchor = "end")
text(x5 + 80, 274, "type I error", size = 15, weight = "bold")
text(x5 + 80, 292, 'under H<tspan font-size="10" dy="3">0</tspan>', size = 12, fill = "#5c6878")
# power
rect(x5, 316, 160, 160, fill = "#ffffff", stroke = "#e1e5ea", rx = 8)
line(x5 + 18, 452, x5 + 146, 452, stroke = "#9aa5b4", sw = 1.2, marker = "arrg")
line(x5 + 18, 452, x5 + 18, 330, stroke = "#9aa5b4", sw = 1.2, marker = "arrg")
t <- seq(0, 1, length.out = 60)
for (k in 1:2) {
  yv <- 1 / (1 + exp(-(t - c(0.4, 0.62)[k]) * 10))
  path(poly(x5 + 22 + t * 118, 446 - yv * 100), stroke = c(col_h1, "#f0a0aa")[k], sw = 3)
}
text(x5 + 140, 468, "signal", size = 11, fill = "#5c6878", anchor = "end")
text(x5 + 80, 504, "power", size = 15, weight = "bold")
text(x5 + 80, 522, 'under H<tspan font-size="10" dy="3">1</tspan>', size = 12, fill = "#5c6878")
# tree-shape surface
rect(x5, 546, 160, 160, fill = "#ffffff", stroke = "#e1e5ea", rx = 8)
vals <- matrix(c(0.9, 0.6, 0.35, 0.75, 0.5, 0.25, 0.55, 0.3, 0.1), 3, byrow = TRUE)
for (i in 1:3) for (j in 1:3) rect(x5 + 40 + (j - 1) * 34, 562 + (i - 1) * 34, 32, 32,
                                    fill = ramp(vals[i, j], col_h1), rx = 4)
line(x5 + 26, 662, x5 + 26, 566, sw = 1.6, marker = "arr")
text(x5 + 18, 620, "γ", size = 16, style = "italic", anchor = "middle")
line(x5 + 40, 680, x5 + 140, 680, sw = 1.6, marker = "arr")
text(x5 + 90, 700, 'I<tspan font-size="10" dy="3">C</tspan>', size = 15, style = "italic")
# replicates
for (k in 2:0) rect(x5 + 38 + k * 8, 724 + k * 6, 70, 46, fill = "#ffffff", stroke = "#b8c0cc",
                    sw = 1.2, rx = 6)
text(x5 + 73, 753, "×1000", size = 16, weight = "bold")
text(x5 + 136, 746, "replicates", size = 11, fill = "#5c6878", anchor = "start")
text(x5 + 136, 760, "per cell", size = 11, fill = "#5c6878", anchor = "start")

add("</svg>")

svg_file <- file.path(out_dir, "fig1_schematic.svg")
writeLines(svg, svg_file, useBytes = TRUE)
rsvg <- Sys.which("rsvg-convert")
if (!nzchar(rsvg)) rsvg <- "/opt/homebrew/bin/rsvg-convert"
system2(rsvg, c("-f", "pdf", "-o", shQuote(file.path(out_dir, "fig1_schematic.pdf")),
                shQuote(svg_file)))
system2(rsvg, c("-w", "2400", "-o", shQuote(file.path(out_dir, "fig1_schematic.png")),
                shQuote(svg_file)))
message("Written: ", svg_file, " (+ .pdf, .png)")
