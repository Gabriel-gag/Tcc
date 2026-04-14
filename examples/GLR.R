
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

data(oil_3w_Type_1)

# Escolha do Poço
#df <- oil_3w_Type_1[["WELL-00001_20140124213136"]]
df <- oil_3w_Type_1[["WELL-00002_20140126200050"]]

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

df_plot <- data.frame(
  Index  = 1:n,
  Series = series,
  GLR    = results_glr,
  Alarm  = results_drift,
  Real   = labels
)

# Alarmes confirmados pelo GLR
df_alarm <- df_plot[df_plot$Alarm == 1, ]

# Eventos reais (ground truth)
df_anom  <- df_plot[df_plot$Real == 1, ]

# =========================================================
# GRÁFICO 1 — SÉRIE TEMPORAL + GLR
# =========================================================

g1 <- ggplot(df_plot, aes(x = Index)) +
  geom_line(aes(y = Series), color = "gray40") +
  
  geom_point(
    data = df_alarm,
    aes(y = Series, color = "Alarme GLR"),
    shape = 17,
    size = 3
  ) +
  
  geom_point(
    data = df_anom,
    aes(y = Series, color = "Anomalia Real"),
    shape = 4,
    size = 3
  ) +
  
  scale_color_manual(
    name = "Legenda",
    values = c(
      "Alarme GLR"    = "darkgreen",
      "Anomalia Real" = "black"
    )
  ) +
  
  labs(
    title = paste("Série:", series_name, "- Alarmes GLR"),
    y = "Valor",
    x = "Índice"
  ) +
  
  theme_minimal() +
  theme(legend.position = "bottom")

# =========================================================
# GRÁFICO 2 — ESTATÍSTICA GLR
# =========================================================

g2 <- ggplot(df_plot, aes(x = Index, y = GLR)) +
  geom_line(color = "darkgreen") +
  
  geom_hline(
    yintercept = GLR_THRESHOLD,
    linetype = "dashed",
    color = "red"
  ) +
  
  labs(
    title = "Estatística GLR",
    y = "GLR Score",
    x = "Índice"
  ) +
  
  theme_minimal()

grid.arrange(g1, g2, ncol = 1, heights = c(1.2, 1))

# =========================================================
# 7. AVALIAÇÃO DE DELAY (ATRASO)
# =========================================================

delay_glr <- sapply(event_idx, function(e) {
  det <- drift_indices[drift_indices >= e][1]
  ifelse(is.na(det), NA, det - e)
})

stats_delay <- summary(delay_glr)

cat("\n================ DELAY (ATRASO) ====================\n")
if(all(is.na(delay_glr))) {
  cat("Nenhum evento detectado para calcular delay.\n")
} else {
  # Verifica se existe o valor no summary antes de imprimir para evitar erro
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