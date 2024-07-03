library(forecast)

datacor <-read.csv("data/timeseries.csv", header=TRUE)


###Sum of variance gradual VS temperature
# fit an ARIMA model
x_model <- auto.arima(datacor$Temperature)
# keep the residuals ("white noise")
x_residuals <- x_model$residuals
# fit the same ARIMA model to the Y-series
# by using the "Arima" function in forecast
y_model <- Arima(datacor$sumvar_g, model=x_model)
# keep the residuals
y_filtered <- residuals(y_model)
# apply the ccf to the residuals
ccf(x_residuals, y_filtered)


###Sum of variance gradual VS CO2
# fit an ARIMA model
x_model <- auto.arima(datacor$CO2)
# keep the residuals ("white noise")
x_residuals <- x_model$residuals
# fit the same ARIMA model to the Y-series
# by using the "Arima" function in forecast
y_model <- Arima(datacor$sumvar_g, model=x_model)
# keep the residuals
y_filtered <- residuals(y_model)
# apply the ccf to the residuals
ccf(x_residuals, y_filtered)


###Centroid distance gradual VS temperature
# fit an ARIMA model
x_model <- auto.arima(datacor$Temperature)
# keep the residuals ("white noise")
x_residuals <- x_model$residuals
# fit the same ARIMA model to the Y-series
# by using the "Arima" function in forecast
y_model <- Arima(datacor$centr_g, model=x_model)
# keep the residuals
y_filtered <- residuals(y_model)
# apply the ccf to the residuals
ccf(x_residuals, y_filtered)


###Centroid distance gradual VS CO2
# fit an ARIMA model
x_model <- auto.arima(datacor$CO2)
# keep the residuals ("white noise")
x_residuals <- x_model$residuals
# fit the same ARIMA model to the Y-series
# by using the "Arima" function in forecast
y_model <- Arima(datacor$centr_g, model=x_model)
# keep the residuals
y_filtered <- residuals(y_model)
# apply the ccf to the residuals
ccf(x_residuals, y_filtered)



###Sum of variance punctuated VS temperature
# fit an ARIMA model
x_model <- auto.arima(datacor$Temperature)
# keep the residuals ("white noise")
x_residuals <- x_model$residuals
# fit the same ARIMA model to the Y-series
# by using the "Arima" function in forecast
y_model <- Arima(datacor$sumvar_p, model=x_model)
# keep the residuals
y_filtered <- residuals(y_model)
# apply the ccf to the residuals
ccf(x_residuals, y_filtered)



###Sum of variance gradual VS CO2
# fit an ARIMA model
x_model <- auto.arima(datacor$CO2)
# keep the residuals ("white noise")
x_residuals <- x_model$residuals
# fit the same ARIMA model to the Y-series
# by using the "Arima" function in forecast
y_model <- Arima(datacor$sumvar_p, model=x_model)
# keep the residuals
y_filtered <- residuals(y_model)
# apply the ccf to the residuals
ccf(x_residuals, y_filtered)


###Centroid distance gradual VS temperature
# fit an ARIMA model
x_model <- auto.arima(datacor$Temperature)
# keep the residuals ("white noise")
x_residuals <- x_model$residuals
# fit the same ARIMA model to the Y-series
# by using the "Arima" function in forecast
y_model <- Arima(datacor$centr_p, model=x_model)
# keep the residuals
y_filtered <- residuals(y_model)
# apply the ccf to the residuals
ccf(x_residuals, y_filtered)


###Centroid distance gradual VS CO2
# fit an ARIMA model
x_model <- auto.arima(datacor$CO2)
# keep the residuals ("white noise")
x_residuals <- x_model$residuals
# fit the same ARIMA model to the Y-series
# by using the "Arima" function in forecast
y_model <- Arima(datacor$centr_p, model=x_model)
# keep the residuals
y_filtered <- residuals(y_model)
# apply the ccf to the residuals
ccf(x_residuals, y_filtered)

