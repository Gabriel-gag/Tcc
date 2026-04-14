
# =========================================================
# BIBLIOTECAS
# =========================================================

library(ggplot2)
library(gridExtra)
library(caret)
library(dplyr)

# =========================================================
# FUNÇÕES AUXILIARES
# =========================================================

glr_statistic <- function(x) {
  N <- length(x)
  best <- 0
  
  for (k in 15:(N - 15)) {
    x1 <- x[1:k]
    x2 <- x[(k + 1):N]
    
    mu1 <- mean(x1)
    mu2 <- mean(x2)
    s2  <- var(x)
    
    if (s2 > 0) {
      stat <- k * (N - k) / N * (mu1 - mu2)^2 / s2
      best <- max(best, stat)
    }
  }
  best
}

# =========================================================
# DADOS
# =========================================================

# Primeiro tenta carregar o arquivo local (completo, baixado manualmente)
# Depois o data() como fallback da versão do pacote
local_data <- "oil_3w_Type_1.RData"
if (file.exists(local_data)) {
  load(local_data)
  cat("Carregado do arquivo local:", local_data, "\n")
} else {
  ("nao carregou")
}

# Escolha do Poço
df <- oil_3w_Type_1[["WELL-00001_20140124213136"]]
#df <- oil_3w_Type_1[["WELL-00002_20140126200050"]]

# Lista de sensores disponíveis
sensor_list <- c("p_pdg", "p_tpt", "t_tpt", "p_mon_ckp", "t_jus_ckp", "p_jus_ckgl", "qgl")

# ---------------------------------------------------------
# SELECIONE O SENSOR AQUI (1 a 7)
# 1: p_pdg      (Pressão Fundo)
# 2: p_tpt      (Pressão Cabeça)
# 3: t_tpt      (Temp. Cabeça)
# 4: p_mon_ckp  (Pressão Montante Choke)
# 5: t_jus_ckp  (Temp. Jusante Choke)
# 6: p_jus_ckgl (Pressão Jusante GL)
# 7: qgl        (Vazão Gas Lift)
# ---------------------------------------------------------
i_sensor <- 2

series_name <- sensor_list[i_sensor]
series <- df[[series_name]]

labels <- as.integer(df$event)
n <- length(series)

# =========================================================
# HIPERPARÂMETROS FIXOS
# =========================================================

WARMUP        <- 500
GLR_WINDOW    <- 60
GLR_THRESHOLD <- 8
MIN_SEG       <- 15

# =========================================================
# MONITORAMENTO
# =========================================================

results_glr   <- rep(NA, n)
results_drift <- rep(0, n)
drift_indices <- c()

for (t in (WARMUP + GLR_WINDOW):n) {
  
  win <- series[(t - GLR_WINDOW + 1):t]
  glr_val <- glr_statistic(win)
  
  results_glr[t] <- glr_val
  
  if (!is.na(glr_val) && glr_val > GLR_THRESHOLD) {
    results_drift[t] <- 1
    drift_indices <- c(drift_indices, t)
  }
}

# =========================================================
# AVALIAÇÃO (DRIFT)
# =========================================================

# Expande rótulos para janela de tolerância
label_drift <- rep(0, n)
event_idx <- which(labels == 1)

for (i in event_idx) {
  label_drift[max(1, i-25):min(n, i+25)] <- 1
}

valid_idx <- which(!is.na(results_glr))

pred_vec <- factor(results_drift[valid_idx], levels = c(0,1))
ref_vec  <- factor(label_drift[valid_idx], levels = c(0,1))

cm <- confusionMatrix(pred_vec, ref_vec, positive = "1")

# -------------------------
# MATRIZ DE CONFUSÃO (LEGÍVEL)
# -------------------------
TN <- cm$table["0","0"]
FP <- cm$table["1","0"]
FN <- cm$table["0","1"]
TP <- cm$table["1","1"]

cat("\n================ MATRIZ DE CONFUSÃO ================\n")
cat("Verdadeiro Negativo (TN):", TN, "\n")
cat("Falso Positivo      (FP):", FP, "\n")
cat("Falso Negativo      (FN):", FN, "\n")
cat("Verdadeiro Positivo (TP):", TP, "\n")

cat("\n================ MÉTRICAS ==========================\n")
cat("Precisão :", round(cm$byClass["Precision"], 4), "\n")
cat("Recall   :", round(cm$byClass["Sensitivity"], 4), "\n")
cat("F1-score :", round(cm$byClass["F1"], 4), "\n")

cat("\nPrimeiro drift detectado em t =",
    ifelse(length(drift_indices) > 0, drift_indices[1], "Nenhum"), "\n")

# =========================================================
# PREPARAÇÃO DOS DADOS PARA GRÁFICOS
# =========================================================

sensor_labels <- c(
  p_pdg      = "Pressão de Fundo (PDG)",
  p_tpt      = "Pressão de Cabeça (TPT)",
  t_tpt      = "Temperatura de Cabeça (TPT)",
  p_mon_ckp  = "Pressão Montante do Choke",
  t_jus_ckp  = "Temperatura Jusante do Choke",
  p_jus_ckgl = "Pressão Jusante GL",
  qgl        = "Vazão de Gás Lift"
)
y_label <- sensor_labels[[series_name]] %||% series_name

df_plot <- data.frame(
  Index  = 1:n,
  Series = series,
  GLR    = results_glr,
  Alarm  = results_drift,
  Real   = labels
)

df_alarm <- df_plot[df_plot$Alarm == 1, ]
df_anom  <- df_plot[df_plot$Real == 1, ]

delay_glr <- sapply(event_idx, function(e) {
  det <- drift_indices[drift_indices >= e][1]
  ifelse(is.na(det), NA, det - e)
})

detected <- sum(!is.na(delay_glr))
total    <- length(event_idx)
f1_val   <- round(cm$byClass["F1"], 4)
prec_val <- round(cm$byClass["Precision"], 4)
rec_val  <- round(cm$byClass["Sensitivity"], 4)
subtitle_txt <- paste0(
  "Detectados: ", detected, "/", total,
  " | F1=", f1_val, " | Prec=", prec_val, " | Recall=", rec_val
)

# =========================================================
# GRÁFICO 1 — SÉRIE TEMPORAL + GLR
# =========================================================

g1 <- ggplot(df_plot, aes(x = Index)) +

  # Sombra de anomalia (mais suave)
  geom_ribbon(
    data = df_plot[df_plot$Real == 1, ],
    aes(ymin = min(series, na.rm = TRUE),
        ymax = max(series, na.rm = TRUE)),
    fill = "#ff6b6b",
    alpha = 0.05
  ) +

  # Série principal (mais destaque)
  geom_line(aes(y = Series),
            color = "#1b4332",
            size = 0.8) +

  # Alarmes GLR (mais discretos)
  geom_point(
    data = df_alarm,
    aes(y = Series),
    color = "#d00000",
    shape = 17,
    size = 2.5
  ) +

  # Anomalias reais (mais importantes)
  geom_point(
    data = df_anom,
    aes(y = Series),
    color = "#1a75ff",
    shape = 16,
    size = 3,
    alpha = 0.8
  ) +

  labs(
    title    = "Série Temporal — Detecção de Mudanças (GLR)",
    subtitle = subtitle_txt,
    x        = "Índice Temporal",
    y        = y_label
  ) +

  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 11),
    panel.grid.minor = element_blank(),   # remove ruído
    panel.grid.major = element_line(color = "gray85"),
    legend.position = "none"
  )

# =========================================================
# GRÁFICO 2 — ESTATÍSTICA GLR
# =========================================================

g2 <- ggplot() +
  geom_line(data = df_plot, aes(x = Index, y = GLR),
            color = "#2c6e49", size = 0.6) +

  geom_hline(
    yintercept = GLR_THRESHOLD,
    linetype = "dashed", color = "red", size = 0.8
  ) +

  labs(
    title = "Estatística GLR",
    x     = "Índice Temporal",
    y     = "Estatística GLR"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80"),
    legend.position  = "none"
  )

#grid.arrange(g1, g2, ncol = 1, heights = c(1.3, 1.1))
grid.arrange(g1)
# =========================================================
# 7. AVALIAÇÃO DE DELAY (ATRASO)
# =========================================================

stats_delay <- summary(delay_glr)

cat("\n================ DELAY (ATRASO) ====================\n")
if(all(is.na(delay_glr))) {
  cat("Nenhum evento detectado para calcular delay.\n")
} else {
  get_stat <- function(s, name) ifelse(name %in% names(s), as.numeric(s[name]), NA)
  cat("Mínimo   :", get_stat(stats_delay, "Min."), "\n")
  cat("Mediana  :", get_stat(stats_delay, "Median"), "\n")
  cat("Média    :", round(get_stat(stats_delay, "Mean"), 2), "\n")
  cat("Máximo   :", get_stat(stats_delay, "Max."), "\n")
}

cat("\n----------------------------------------------------\n")
cat("Eventos reais totais   :", length(event_idx), "\n")
cat("Eventos detectados     :", sum(!is.na(delay_glr)), "\n")
cat("Eventos não detectados :", sum(is.na(delay_glr)), "\n")


# =========================================================
# 8. AVALIAÇÃO COM TOLERÂNCIA (+/- 250)
# =========================================================

TOL <- 250

pred_tol_glr <- rep(0, n)
# Preenche "1" ao redor das detecções do GLR com a margem de tolerância
for (d in drift_indices) {
  pred_tol_glr[max(1, d - TOL):min(n, d + TOL)] <- 1
}

ref_tol_glr <- label_drift # Usa a referência padrão definida acima

cm_tol_glr <- confusionMatrix(
  factor(pred_tol_glr, levels = c(0,1)),
  factor(ref_tol_glr, levels = c(0,1)),
  positive = "1"
)

TN_tol <- cm_tol_glr$table["0","0"]
FP_tol <- cm_tol_glr$table["1","0"]
FN_tol <- cm_tol_glr$table["0","1"]
TP_tol <- cm_tol_glr$table["1","1"]

cat("\n================ MATRIZ (COM TOLERÂNCIA +/- 250) ===\n")
cat("Verdadeiro Negativo (TN):", TN_tol, "\n")
cat("Falso Positivo      (FP):", FP_tol, "\n")
cat("Falso Negativo      (FN):", FN_tol, "\n")
cat("Verdadeiro Positivo (TP):", TP_tol, "\n")

cat("\n================ MÉTRICAS (COM TOLERÂNCIA) =========\n")
cat("Precisão :", round(cm_tol_glr$byClass["Precision"], 4), "\n")
cat("Recall   :", round(cm_tol_glr$byClass["Sensitivity"], 4), "\n")
cat("F1-score :", round(cm_tol_glr$byClass["F1"], 4), "\n")
cat("====================================================\n")


#Deteccao

event_detected <- sapply(event_idx, function(e) {
  any(drift_indices >= e & drift_indices <= e + TOL)
})

event_recall <- mean(event_detected)

#Atraso

delay <- sapply(event_idx, function(e) {
  det <- drift_indices[drift_indices >= e][1]
  ifelse(is.na(det), NA, det - e)
})
summary(delay)

#Falso alarme
FAR_per_1000 <- sum(results_drift == 1) / n * 1000
print(FAR_per_1000)
print(event_recall)
