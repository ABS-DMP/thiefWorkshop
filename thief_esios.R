# ============================================================
# THieF aplicado al mercado eléctrico español
# Datos ESIOS / MIBGAS
# ============================================================


# ============================================================
# 1. LIBRERÍAS
# ============================================================

library(MAPA)
library(tsutils)


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

# Creamos una serie para cada nivel de agregación:
# k = 1  -> 1 hora
# k = 2  -> bloques de 2 horas
# k = 3  -> bloques de 3 horas
# ...
# k = 24 -> bloques de 24 horas

Y <- tsaggr(
  precio_ts,
  k,
  FALSE,
  FALSE
)

# Inspeccionar el resultado
print(Y)