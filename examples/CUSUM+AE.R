
# =========================================================
# BIBLIOTECAS
# =========================================================

library(daltoolbox)
library(daltoolboxdp)
library(harbinger)
library(united)
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

data(oil_3w_Type_1)

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

train_data <- series[1:WARMUP]
ts_train   <- ts_data_window(train_data, params$WINDOW_SIZE)

norm_model <- fit(ts_norm_gminmax(), ts_train)
ts_train_n <- transform(norm_model, ts_train)

ae_model <- autoenc_ed(params$WINDOW_SIZE, params$LATENT_SIZE)
ae_model <- fit(ae_model, ts_train_n, epochs = params$EPOCHS_INIT)

rec_train <- transform(ae_model, ts_train_n)
mse_train <- rowMeans((as.matrix(ts_train_n) - as.matrix(rec_train))^2)

tau <- quantile(mse_train, probs = params$PROB_TAU)

# =========================================================
# MONITORAMENTO
# =========================================================

results_mse   <- rep(NA, n)
results_tau   <- rep(NA, n)
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

df_plot <- data.frame(
  Index     = 1:n,
  Series    = series,
  ErrorNorm = error_norm,
  Alarm     = results_drift,
  Real      = labels
)

# Alarmes confirmados pelo CUSUM (h_high)
df_alarm <- df_plot[df_plot$Alarm == 1, ]

# Eventos reais (ground truth)
df_anom <- df_plot[df_plot$Real == 1, ]

# =========================================================
# GRÁFICO 1 — SÉRIE TEMPORAL + CUSUM
# =========================================================

g1 <- ggplot(df_plot, aes(x = Index)) +
  geom_line(aes(y = Series), color = "gray40") +
  
  geom_point(
    data = df_alarm,
    aes(y = Series, color = "Alarme CUSUM"),
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
      "Alarme CUSUM"  = "blue",
      "Anomalia Real" = "black"
    )
  ) +
  
  labs(
    title = "Série temporal com alarmes do CUSUM sobre o erro do Autoencoder",
    y = "Pressão",
    x = "Índice"
  ) +
  
  theme_minimal() +
  theme(legend.position = "bottom")

# =========================================================
# GRÁFICO 2 — ERRO NORMALIZADO (CUSUM)
# =========================================================

g2 <- ggplot(df_plot, aes(x = Index, y = ErrorNorm)) +
  geom_line(color = "red") +
  
  geom_hline(
    yintercept = c(3, -3),
    linetype = "dashed",
    color = "gray40"
  ) +
  
  geom_point(
    data = df_alarm,
    aes(y = ErrorNorm),
    color = "blue",
    shape = 17,
    size = 3
  ) +
  
  coord_cartesian(ylim = c(-5, 5)) +
  
  labs(
    title = "Erro de reconstrução normalizado utilizado pelo CUSUM",
    y = "Erro normalizado",
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