library(daltoolbox)
library(tspredit)
library(daltoolboxdp)
library(ggplot2)
library(united)


data(oil_3w_Type_1)
data <- oil_3w_Type_1[[1]]

series <- data$p_tpt
plot (as.ts(series))

model <- hcp_binseg()
model <- fit(model, series)

detection <- detect(model, series)

print(detection |> dplyr::filter(event == TRUE))

grf <- har_plot(model, series, detection, data$event)
plot(grf)

ev <- evaluate(model, detection$event, data$event)
print(ev$confMatrix)

ev$accuracy
ev$F1

ev_soft <- evaluate(har_eval_soft(sw = 90), detection$event, data$event)
print(ev_soft$confMatrix)

ev_soft$accuracy
ev_soft$F1


