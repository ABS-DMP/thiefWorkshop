# ============================================================
# THieF aplicado al mercado eléctrico español
# Datos ESIOS / MIBGAS
# ============================================================


# ============================================================
# 1. LIBRERÍAS
# ============================================================

library(MAPA)
library(tsutils)
library(forecast)


# ============================================================
# 2. CARGA DE DATOS DE PRECIO
# ============================================================

data_dir <- "data/raw"

precio_raw <- read.csv(
  file.path(data_dir, "precio_2020_2026.csv"),
  check.names = FALSE
)


# ============================================================
# 3. SELECCIÓN DE ESPAÑA
# ============================================================

precio_es <- subset(precio_raw, geo_name == "España")

# Convertimos el timestamp UTC de texto a fecha-hora
precio_es$datetime_utc <- as.POSIXct(
  precio_es$datetime_utc,
  format = "%Y-%m-%dT%H:%M:%SZ",
  tz = "UTC"
)

# Ordenamos cronológicamente
precio_es <- precio_es[order(precio_es$datetime_utc), ]


# ============================================================
# 4. COMPROBACIONES DE CALIDAD
# ============================================================

cat("Observaciones:", nrow(precio_es), "\n")
cat("Duplicados:", sum(duplicated(precio_es$datetime_utc)), "\n")
cat("NA precios:", sum(is.na(precio_es$value)), "\n")

# Comprobar continuidad horaria
saltos <- diff(precio_es$datetime_utc)
print(table(saltos))

# Rango temporal
print(range(precio_es$datetime_utc))


# ============================================================
# 5. CREAR SERIE TEMPORAL HORARIA
# ============================================================

# frequency = 24:
# consideramos un ciclo diario compuesto por 24 observaciones horarias
precio_ts <- ts(
  precio_es$value,
  frequency = 24
)

cat("Frecuencia:", frequency(precio_ts), "\n")
cat("Longitud:", length(precio_ts), "\n")


# ============================================================
# 6. VISUALIZACIÓN INICIAL
# ============================================================

# Mostramos las últimas dos semanas
plot(
  tail(precio_ts, 24 * 14),
  type = "l",
  main = "Precio SPOT España - últimas 2 semanas",
  xlab = "Hora",
  ylab = "Precio"
)


# ============================================================
# 7. NIVELES DE AGREGACIÓN TEMPORAL
# ============================================================

# Mismo procedimiento utilizado en el ejemplo de Kourentzes
f <- frequency(precio_ts)

# Calculamos los divisores válidos de la frecuencia
k <- rev(f / (1:f))
k <- k[k %% 1 == 0]

# Número de niveles de agregación
p <- length(k)

cat("Frecuencia original:", f, "\n")
cat("Niveles de agregación:", k, "\n")
cat("Número de niveles:", p, "\n")


# ============================================================
# 8. AGREGACIÓN TEMPORAL
# ============================================================

# IMPORTANTE:
# fmean = TRUE -> agregamos mediante la media
# Para precios tiene más sentido económico que sumar precios.

Y <- tsaggr(
  precio_ts,
  k,
  TRUE,
  FALSE
)

if (exists("Y")) {
  cat("Agregación temporal completada.\n")
}
# Inspeccionar el resultado
# print(Y)

# ============================================================
# 9. HORIZONTE DE PREDICCIÓN
# ============================================================

# Queremos predecir las siguientes 24 horas
h <- 24

# Horizonte equivalente para cada nivel de agregación
H <- h / k

cat("Horizonte original:", h, "horas\n")

for (j in 1:p) {
  cat(
    "k =", k[j],
    "| horizonte =", H[j],
    "observaciones\n"
  )
}

# ============================================================
# 10. BASE FORECASTS: ARIMA EN CADA NIVEL
# ============================================================

forecast_file <- "data/processed/frc_arima_base.rds"

if (file.exists(forecast_file)) {

  cat("Cargando forecasts ARIMA guardados...\n")
  frc <- readRDS(forecast_file)

} else {

  cat("Calculando forecasts ARIMA...\n")

  frc <- list()

  for (j in 1:p) {

    temp <- Y[[1]][[j]]
    n <- length(temp)

    tempTrn <- head(temp, n - H[j])

    fit <- auto.arima(tempTrn)

    frc[[j]] <- forecast(
      fit,
      h = H[j]
    )$mean
  }

  dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
  saveRDS(frc, forecast_file)
}

# ============================================================
# 11. MATRIZ S PARA PRECIOS MEDIOS
# ============================================================

# 11. MATRIZ S para THieF
# La S debe reflejar cómo cada nivel agregado se relaciona con la serie original
S <- tsutils::Sthief(Y[[1]][[1]])
# Si en tu objeto aparece nombrado, también vale:
# S <- tsutils::Sthief(Y[[1]]$AL1)

# comprobación
dim(S)

# estimación de W y G
W <- diag(1 / rowSums(S))
G <- solve(t(S) %*% W %*% S) %*% t(S) %*% W

# base forecasts apilados
fbase <- cbind(unlist(rev(frc)))

# reconciliación THieF
freco <- S %*% G %*% fbase