
# =========================================================
# BIBLIOTECAS
# =========================================================

library(daltoolbox)
library(daltoolboxdp)
library(harbinger)
#library(united)
library(ggplot2)
library(gridExtra)
library(caret)
library(dplyr)


# =========================================================
# AJUSTE DE TUNINGS
# =========================================================

# -------------------------
# Tuning 1 (Conservador)
# -------------------------
tuning_1 <- list(
  WINDOW_SIZE = 10,
  LATENT_SIZE = 3,
  EPOCHS_INIT = 30,
  EPOCHS_RET  = 15,
  K_FACTOR    = 0.25,
  H_LOW       = 3,
  H_HIGH      = 5,
  PROB_TAU    = 0.995
)

# -------------------------
# Tuning 2 (Mais Sensível)
# -------------------------
tuning_2 <- list(
  WINDOW_SIZE = 8,
  LATENT_SIZE = 2,
  EPOCHS_INIT = 25,
  EPOCHS_RET  = 10,
  K_FACTOR    = 0.15,
  H_LOW       = 2.5,
  H_HIGH      = 4,
  PROB_TAU    = 0.99
)

# -------------------------
# Tuning 3 (Mais Estável)
# -------------------------
tuning_3 <- list(
  WINDOW_SIZE = 12,
  LATENT_SIZE = 4,
  EPOCHS_INIT = 40,
  EPOCHS_RET  = 20,
  K_FACTOR    = 0.35,
  H_LOW       = 3.5,
  H_HIGH      = 6,
  PROB_TAU    = 0.997
)

# Escolha do tuning
TuningEscolhido <- 3
TUNINGS <- list(tuning_1, tuning_2, tuning_3)
params <- TUNINGS[[TuningEscolhido]]

# =========================================================
# FUNÇÕES AUXILIARES
# =========================================================

ts_data_window <- function(data, sw) {
  n <- length(data)
  if (n < sw) return(NULL)
  mat <- matrix(NA, nrow = n - sw + 1, ncol = sw)
  for (i in 1:nrow(mat)) {
    mat[i, ] <- data[i:(i + sw - 1)]
  }
  as.data.frame(mat)
}

ts_norm_gminmax <- function() {
  structure(list(min = NULL, max = NULL), class = "ts_norm_gminmax")
}

fit.ts_norm_gminmax <- function(obj, data) {
  m <- as.matrix(data)
  obj$min <- min(m, na.rm = TRUE)
  obj$max <- max(m, na.rm = TRUE)
  obj
}

transform.ts_norm_gminmax <- function(obj, data) {
  m <- as.matrix(data)
  as.data.frame((m - obj$min) / (obj$max - obj$min + 1e-8))
}

cusum_cross_tc_step <- function(x, k, h_low, h_high, pos, neg) {
  pos <- max(0, pos + x - k)
  neg <- min(0, neg + x + k)

  warning <- (pos > h_low) | (neg < -h_low)
  drift   <- (pos > h_high) | (neg < -h_high)

  list(pos = pos, neg = neg, warning = warning, drift = drift)
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
# 2: p_tpt      (Pressão Cabeça)#
# 3: t_tpt      (Temp. Cabeça)###
# 4: p_mon_ckp  (Pressão Montante Choke)###
# 5: t_jus_ckp  (Temp. Jusante Choke)
# 6: p_jus_ckgl (Pressão Jusante GL)###
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

WARMUP       <- 500
RETRAIN_SIZE <- 300

# =========================================================
# TREINAMENTO INICIAL
# =========================================================
library(zoo)

# 1. interpolação linear
series <- zoo::na.approx(series, na.rm = FALSE)

# 2. forward fill (preenche para frente)
series <- zoo::na.locf(series, na.rm = FALSE)

# 3. backward fill (preenche para trás)
series <- zoo::na.locf(series, fromLast = TRUE)

# 🔍 verificação
cat("NA após tratamento:", sum(is.na(series)), "\n")

train_data <- series[1:WARMUP]
ts_train   <- ts_data_window(train_data, params$WINDOW_SIZE)
sum(is.na(ts_train))

norm_model <- fit(ts_norm_gminmax(), ts_train)
ts_train_n <- transform(norm_model, ts_train)

ae_model <- autoenc_ed(params$WINDOW_SIZE, params$LATENT_SIZE)
ae_model <- fit(ae_model, ts_train_n, epochs = params$EPOCHS_INIT)

rec_train <- transform(ae_model, ts_train_n)
mse_train <- rowMeans((as.matrix(ts_train_n) - as.matrix(rec_train))^2)
sum(is.na(mse_train))
sum(is.nan(mse_train))
summary(mse_train)
sum(is.na(ts_train))
sum(is.na(ts_train_n))
summary(ts_train_n)
tau <- quantile(mse_train, probs = params$PROB_TAU)

# =========================================================
# MONITORAMENTO
# =========================================================

results_mse   <- rep(NA, n)
results_tau   <- rep(NA, n)
results_warning <- rep(0, n)
results_drift <- rep(0, n)
drift_indices <- c()

pos <- 0
neg <- 0
last_retrain <- WARMUP

for (t in (WARMUP + 1):n) {

  win <- series[(t - params$WINDOW_SIZE + 1):t]
  win_df <- as.data.frame(t(win))
  win_n  <- transform(norm_model, win_df)

  rec <- transform(ae_model, win_n)
  mse_t <- mean((as.matrix(win_n) - as.matrix(rec))^2)

  results_mse[t] <- mse_t
  results_tau[t] <- tau

  if (t > last_retrain + 50) {

    buffer <- na.omit(results_mse[(t-100):(t-1)])

    if (length(buffer) > 20 && sd(buffer) > 0) {

      z <- (mse_t - mean(buffer)) / (sd(buffer) + 1e-5)

      cs <- cusum_cross_tc_step(
        z,
        params$K_FACTOR,
        params$H_LOW,
        params$H_HIGH,
        pos,
        neg
      )

      pos <- cs$pos
      neg <- cs$neg
      if (cs$warning) {
        results_warning[t] <- 1
      }
      if (cs$drift && (t - last_retrain > RETRAIN_SIZE)) {

        results_drift[t] <- 1
        drift_indices <- c(drift_indices, t)

        new_data <- series[(t - RETRAIN_SIZE + 1):t]
        ts_new   <- ts_data_window(new_data, params$WINDOW_SIZE)

        norm_model <- fit(ts_norm_gminmax(), ts_new)
        ts_new_n   <- transform(norm_model, ts_new)

        ae_model <- fit(ae_model, ts_new_n, epochs = params$EPOCHS_RET)

        rec_new <- transform(ae_model, ts_new_n)
        mse_new <- rowMeans((as.matrix(ts_new_n) - as.matrix(rec_new))^2)

        tau <- quantile(mse_new, probs = params$PROB_TAU)

        pos <- 0
        neg <- 0
        last_retrain <- t
      }
    }
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

valid_idx <- which(!is.na(results_mse))

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

# Erro normalizado (Z-score usado no CUSUM)
error_norm <- rep(NA, n)

for (t in 1:n) {
  if (t > 150 && !is.na(results_mse[t])) {
    buf <- na.omit(results_mse[(t-100):(t-1)])
    if (length(buf) > 20 && sd(buf) > 0) {
      error_norm[t] <- (results_mse[t] - mean(buf)) / (sd(buf) + 1e-5)
    }
  }
}

# Sensor label dinâmico
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
  Index     = 1:n,
  Series    = series,
  ErrorNorm = error_norm,
  Alarm     = results_drift,
  Real      = labels
)

# Alarmes confirmados pelo CUSUM (h_high)
df_alarm <- df_plot[df_plot$Alarm == 1, ]
df_alarm_low <- df_plot[df_plot$Alarm == 0 & results_warning == 1, ]
# Eventos reais (ground truth)
df_anom <- df_plot[df_plot$Real == 1, ]

# Métricas para exibir no gráfico
delay_glr <- sapply(event_idx, function(e) {
  det <- drift_indices[drift_indices >= e][1]
  ifelse(is.na(det), NA, det - e)
})
detected   <- sum(!is.na(delay_glr))
total      <- length(event_idx)
f1_val     <- round(cm$byClass["F1"], 4)
prec_val   <- round(cm$byClass["Precision"], 4)
rec_val    <- round(cm$byClass["Sensitivity"], 4)
subtitle_txt <- paste0(
  "Tuning ", TuningEscolhido, " - ",
  "Detecção de drift com CUSUM + Autoencoder"
)

# =========================================================
# GRÁFICO 1 — SÉRIE TEMPORAL + CUSUM
# =========================================================

g1 <- ggplot() +
  geom_line(data = df_plot, aes(x = Index, y = Series),
            color = "steelblue4", size = 1, alpha = 0.8) +
  
  # Regiões de anomalia real (sombreado)
  geom_ribbon(
    data = df_plot[df_plot$AlarmLow == 1, ],
    aes(x = Index,
        ymin = min(series),
        ymax = max(series)),
    fill = "orange",
    alpha = 0.05
  ) +
  # 2º: Alarme do CUSUM (Low TC)
  geom_point(
    data = df_alarm_low, # Atenção: certifique-se de que essa variável existe no seu código
    aes(x = Index, y = Series, color = "Alarme Low TC"),
    shape = 16, size = 2, alpha = 0.4 # 18 = Losango cheio
  ) +
  
  
  # 1º: Alarme do CUSUM (High TC)
  geom_point(
    data = df_alarm, # Atenção: certifique-se de que essa variável existe no seu código
    aes(x = Index, y = Series, color = "Alarme High TC"),
    shape = 17, size = 4 # 17 = Triângulo cheio
  ) +
  
  
  # 3º: Pontos de anomalia real (Por último para ficar por cima)
  geom_point(
    data = df_anom,
    aes(x = Index, y = Series, color = "Anomalia Real"),
    shape = 16, size = 4, alpha = 1 # alpha = 1 para não vazar a cor de baixo
  ) +
  
  # Atualização da Legenda com as 3 categorias
  scale_color_manual(
    name = "Eventos",
    values = c(
      "Alarme High TC" = "red",
      "Alarme Low TC"  = "orange", # Escolha a cor que preferir para o Low TC
      "Anomalia Real"  = "#1a75ff"
    ),
    guide = guide_legend(override.aes = list(
      shape = c(17, 18, 16), # Ordem alfabética: Alarme High, Alarme Low, Anomalia
      size  = c(4, 4, 3),
      alpha = c(1, 1, 1)
    ))
  ) +
  
  labs(
    title    = "Série Temporal — Detecção de Mudanças (CUSUM + Autoencoder)",
    subtitle = subtitle_txt,
    x        = "Índice Temporal",
    y        = y_label
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 10, color = "gray30"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80")
  )

# =========================================================
# GRÁFICO 2 — ERRO NORMALIZADO + THRESHOLDS CUSUM
# =========================================================

g2 <- ggplot() +
  geom_line(data = df_plot, aes(x = Index, y = ErrorNorm),
            color = "#d63031", size = 0.6, alpha = 0.85) +

  # Threshold H_HIGH
  geom_hline(
    yintercept = c(params$H_HIGH, -params$H_HIGH),
    linetype = "dashed", color = "red", size = 0.9, alpha = 0.8
  ) +

  # Threshold H_LOW
  geom_hline(
    yintercept = c(params$H_LOW, -params$H_LOW),
    linetype = "dotted", color = "orange", size = 0.7, alpha = 0.7
  ) +

  # Zona de aviso (entre h_low e h_high)
  annotate("rect",
    xmin = -Inf, xmax = Inf,
    ymin = params$H_LOW, ymax = params$H_HIGH,
    fill = "orange", alpha = 0.06
  ) +
  annotate("rect",
    xmin = -Inf, xmax = Inf,
    ymin = -params$H_HIGH, ymax = -params$H_LOW,
    fill = "orange", alpha = 0.06
  ) +

  # Alarmes detectados
  geom_point(
    data = df_alarm,
    aes(x = Index, y = ErrorNorm),
    color = "red", shape = 17, size = 3.5
  ) +

  scale_y_continuous(limits = c(-5, 5), expand = c(0, 0)) +

  labs(
    title    = "Erro de Reconstrução Normalizado (Z-score) — CUSUM",
    subtitle = paste0(
      "--- Zona de Alerta (h_low = ", params$H_LOW, ")  ---  --- Zona de Drift (h_high = ", params$H_HIGH, ") ---"
    ),
    x = "Índice Temporal",
    y = "Z-score"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "gray40", face = "italic"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80"),
    legend.position  = "none"
  )

grid.arrange(g1, g2, ncol = 1, heights = c(1.3, 1.1))

# =========================================================
# 7. AVALIAÇÃO DE DELAY (ATRASO)
# =========================================================

stats_delay <- summary(delay_glr)
detected   <- sum(!is.na(delay_glr))

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

pred_tol <- rep(0, n)
# Preenche "1" ao redor das detecções com a margem de tolerância
for (d in drift_indices) {
  pred_tol[max(1, d - TOL):min(n, d + TOL)] <- 1
}

ref_tol <- label_drift

cm_tol <- confusionMatrix(
  factor(pred_tol, levels = c(0,1)),
  factor(ref_tol, levels = c(0,1)),
  positive = "1"
)

TN_tol <- cm_tol$table["0","0"]
FP_tol <- cm_tol$table["1","0"]
FN_tol <- cm_tol$table["0","1"]
TP_tol <- cm_tol$table["1","1"]

cat("\n================ MATRIZ (COM TOLERÂNCIA +/- 250) ===\n")
cat("Verdadeiro Negativo (TN):", TN_tol, "\n")
cat("Falso Positivo      (FP):", FP_tol, "\n")
cat("Falso Negativo      (FN):", FN_tol, "\n")
cat("Verdadeiro Positivo (TP):", TP_tol, "\n")

cat("\n================ MÉTRICAS (COM TOLERÂNCIA) =========\n")
cat("Precisão :", round(cm_tol$byClass["Precision"], 4), "\n")
cat("Recall   :", round(cm_tol$byClass["Sensitivity"], 4), "\n")
cat("F1-score :", round(cm_tol$byClass["F1"], 4), "\n")
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
