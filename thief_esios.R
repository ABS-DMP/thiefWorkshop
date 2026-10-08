# ============================================================
# THieF aplicado al precio electrico de Espana (ESIOS)
# Referencia: Forecasting with Temporal Hierarchies, Kourentzes,
# ISF 2025. Los numeros de diapositiva se indican en cada bloque.
#
# Uso: source("thief_esios.R")
# Opciones, ANTES de source():
# options(thief.config = list(
#   mode = "holdout",          # "future": siguientes 24 horas reales
#   last_hours = NULL,         # NULL: todos los datos; 24*60: prueba corta
#   calibration_days = 0L,     # >=20: backtest diario + intervalos empiricos
#   include_mapa = FALSE,      # MAPA original por estados (coste adicional)
#   arima_args = list()        # argumentos adicionales de auto.arima
# ))
# Cada configuracion necesita su propia cache si cambia el entrenamiento.
# Una cache existente se carga o se rechaza; nunca se reentrena en silencio.
# Para pruebas/importar funciones: options(thief.autorun = FALSE).
# ============================================================

thief_script_dir <- local({
  files <- vapply(sys.frames(), function(x)
    if (is.null(x$ofile)) "" else as.character(x$ofile)[1], "")
  files <- files[nzchar(files)]
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  path <- if (length(files)) tail(files, 1) else
    if (length(cli)) sub("^--file=", "", cli[1]) else "thief_esios.R"
  dirname(normalizePath(path, mustWork = TRUE))
})

# 1. FUNCIONES AUXILIARES
# ------------------------------------------------------------
thief_save_cache <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("thief-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("No se pudo guardar: ", path)
}

thief_aggregate <- function(y, k) {
  MAPA::tsaggr(y = y, fout = k, fmean = TRUE, outplot = FALSE)[[1]]
}

thief_valid_frc <- function(x, H) {
  is.list(x) && length(x) == length(H) &&
    identical(as.integer(lengths(x)), as.integer(H)) &&
    all(vapply(x, function(z)
      is.numeric(z) && all(is.finite(z)), logical(1)))
}

thief_fit <- function(y_train, k, H, arima_args) {
  # Agregar el entrenamiento por separado evita usar valores del test.
  train <- thief_aggregate(y_train, k)
  frc <- residuals <- vector("list", length(k))
  model_rows <- vector("list", length(k))
  for (j in seq_along(k)) {
    message("ARIMA: k=", k[j], " horas; n=", length(train[[j]]))
    fit <- do.call(forecast::auto.arima, c(list(y = train[[j]]), arima_args))
    frc[[j]] <- forecast::forecast(fit, h = H[j])$mean
    residuals[[j]] <- as.numeric(stats::residuals(fit))
    model_rows[[j]] <- data.frame(
      k = k[j], n_train = length(train[[j]]),
      model = paste(names(forecast::arimaorder(fit)),
                    forecast::arimaorder(fit), collapse = "; "),
      AICc = if (is.null(fit$aicc)) NA_real_ else fit$aicc,
      sigma2 = fit$sigma2
    )
  }
  attr(frc, "residuals") <- residuals
  attr(frc, "models") <- do.call(rbind, model_rows)
  frc
}

thief_load_or_fit <- function(path, context, y_train, k, H, arima_args) {
  if (file.exists(path)) {
    message("Cargando forecasts guardados: ", path)
    frc <- readRDS(path)
    if (!thief_valid_frc(frc, H)) stop("Cache incompatible: ", path)
    saved <- attr(frc, "cache_context")
    if (is.null(saved)) {
      warning("Cache antigua sin metadatos: no se puede verificar su origen.")
    } else if (!identical(saved, context)) {
      stop("La cache corresponde a otros datos/configuracion. Use otra ruta ",
           "forecast_file o archive manualmente la cache: ", path)
    }
    return(frc)
  }
  frc <- thief_fit(y_train, k, H, arima_args)
  stopifnot(thief_valid_frc(frc, H))
  attr(frc, "cache_context") <- context
  thief_save_cache(frc, path)
  frc
}

thief_residual_matrix <- function(residuals, k, f) {
  if (is.null(residuals)) return(NULL)
  # Diap. 28, 30, 39: filas = jerarquias completas, columnas = nodos.
  # tsaggr elimina sobrantes al principio; alineamos todos por el final.
  days <- min(floor(lengths(residuals) / (f / k)))
  if (days < 4) return(NULL)
  E <- do.call(cbind, lapply(rev(seq_along(k)), function(j)
    matrix(tail(residuals[[j]], days * f / k[j]),
           ncol = f / k[j], byrow = TRUE)))
  # El primer ciclo se descarta por la inicializacion del modelo.
  E <- E[-1, , drop = FALSE]
  E <- E[apply(E, 1, function(z) all(is.finite(z))), , drop = FALSE]
  if (nrow(E) < 3) NULL else E
}

thief_covariances <- function(E, n_hours) {
  # Diap. 35-40. W es covarianza, no su inversa.
  # Para medias: W_struct = D^-1 diag(n_hours) D^-1 = diag(1/n_hours).
  Ws <- list(THieF_struct = diag(1 / n_hours),
             THieF_OLS = diag(length(n_hours)))
  lambda <- NA_real_
  ridge <- NA_real_
  if (!is.null(E)) {
    # Aproximacion con residuos un paso de los modelos de entrenamiento.
    # No equivale a conocer la covarianza real del error a 24 pasos.
    centered <- scale(E, center = TRUE, scale = FALSE)
    n <- nrow(E)
    C <- crossprod(centered) / n
    floor_var <- max(mean(diag(C)), 1) * 1e-8
    variances <- pmax(diag(C), floor_var)
    D <- diag(variances)
    # Contraccion de correlaciones hacia cero, intensidad tipo
    # Schafer-Strimmer. Formula de referencia: hts/R/MinT.R.
    Z <- sweep(centered, 2, sqrt(variances), "/")
    R <- crossprod(Z) / n
    sampling_var <- (crossprod(Z^2) - crossprod(Z)^2 / n) / (n * (n - 1))
    off <- row(R) != col(R)
    denominator <- sum(R[off]^2)
    lambda <- if (denominator <= .Machine$double.eps) 1 else
      max(0, min(1, sum(sampling_var[off]) / denominator))
    shrunk <- (1 - lambda) * C + lambda * D
    # Regularizacion numerica explicita incluso ante residuos constantes.
    ridge <- max(0, floor_var - min(eigen(shrunk, symmetric = TRUE,
                                         only.values = TRUE)$values))
    Ws$THieF_variance <- D
    Ws$THieF_MinT <- shrunk + diag(ridge, ncol(E))
  }
  list(W = Ws, lambda = lambda, ridge = ridge)
}

thief_reconcile <- function(S, W, fbase) {
  # Diap. 33-35: G=(S' W^-1 S)^-1 S' W^-1; no invertir W explicitamente.
  W_S <- solve(W, S)
  G <- solve(crossprod(S, W_S), t(W_S))
  stopifnot(max(abs(G %*% S - diag(ncol(S)))) < 1e-6)
  list(G = G, forecast = as.numeric(S %*% G %*% fbase))
}

thief_metrics <- function(actual, predicted, train, seasonal_lag) {
  error <- actual - predicted
  scale_errors <- diff(as.numeric(train), lag = seasonal_lag)
  mae_scale <- mean(abs(scale_errors))
  mse_scale <- mean(scale_errors^2)
  c(ME = mean(error), MAE = mean(abs(error)),
    RMSE = sqrt(mean(error^2)),
    MASE = if (is.finite(mae_scale) && mae_scale > 0)
      mean(abs(error)) / mae_scale else NA_real_,
    RMSSE = if (is.finite(mse_scale) && mse_scale > 0)
      sqrt(mean(error^2) / mse_scale) else NA_real_)
}

thief_png <- function(path, width, height, draw) {
  grDevices::png(path, width = width, height = height, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}

# 2. FLUJO PRINCIPAL
# ------------------------------------------------------------
run_thief_esios <- function(config = list()) {
  defaults <- list(
    project_dir = thief_script_dir,
    input_file = "data/raw/precio_2020_2026.csv",
    mode = "holdout", last_hours = NULL,
    forecast_file = NULL, output_dir = NULL,
    calibration_days = 0L, simulations = 2000L,
    include_mapa = FALSE, arima_args = list(), seed = 2025L
  )
  unknown <- setdiff(names(config), names(defaults))
  if (length(unknown)) stop("Opciones desconocidas: ", paste(unknown, collapse = ", "))
  cfg <- utils::modifyList(defaults, config)
  cfg$mode <- match.arg(cfg$mode, c("holdout", "future"))
  if (!is.list(cfg$arima_args) || any(c("y", "x") %in% names(cfg$arima_args)))
    stop("arima_args debe ser una lista de opciones de auto.arima, sin y/x.")
  if (length(cfg$calibration_days) != 1 ||
      !is.finite(cfg$calibration_days) || cfg$calibration_days < 0 ||
      cfg$calibration_days %% 1 != 0) stop("calibration_days debe ser entero >=0.")
  stopifnot(cfg$simulations >= 100, cfg$simulations %% 1 == 0)
  for (pkg in c("MAPA", "tsutils", "forecast")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Falta el paquete ", pkg, ". Instalelo antes de ejecutar el script.")
  }
  full_path <- function(path) {
    if (grepl("^([A-Za-z]:|/|\\\\)", path)) path else file.path(cfg$project_dir, path)
  }
  input_file <- full_path(cfg$input_file)
  suffix <- if (cfg$mode == "future") "_future" else ""
  forecast_file <- full_path(if (is.null(cfg$forecast_file))
    paste0("data/processed/frc_arima_base", suffix, ".rds") else cfg$forecast_file)
  output_dir <- full_path(if (is.null(cfg$output_dir))
    paste0("output/thief_esios", suffix) else cfg$output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 3. DATOS Y CALIDAD: los bloques se anclan al ultimo instante UTC.
  # encoding marca UTF-8 sin recodificar a la configuracion regional de Windows.
  raw <- read.csv(input_file, check.names = FALSE, encoding = "UTF-8")
  required <- c("geo_name", "datetime_utc", "value")
  if (!all(required %in% names(raw))) stop("Faltan columnas: ", paste(required, collapse=", "))
  precio_es <- raw[!is.na(raw$geo_name) & raw$geo_name == "Espa\u00f1a", ]
  precio_es$datetime_utc <- as.POSIXct(precio_es$datetime_utc,
                                     format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  precio_es <- precio_es[order(precio_es$datetime_utc), ]
  if (nrow(precio_es) < 24 * 10 || anyNA(precio_es$datetime_utc) ||
      !is.numeric(precio_es$value) || any(!is.finite(precio_es$value)) ||
      anyDuplicated(precio_es$datetime_utc) ||
      any(as.numeric(diff(precio_es$datetime_utc), units = "secs") != 3600))
    stop("Se necesitan al menos 10 bloques diarios, precios finitos y horas UTC continuas.")
  if (!is.null(cfg$last_hours)) {
    if (length(cfg$last_hours) != 1 || !is.finite(cfg$last_hours) ||
        cfg$last_hours < 240 || cfg$last_hours %% 1 != 0)
      stop("last_hours debe ser NULL o un entero >=240.")
    precio_es <- tail(precio_es, cfg$last_hours)
  }
  precio_ts <- ts(precio_es$value, frequency = 24)
  message("Observaciones usadas: ", length(precio_ts), "; rango UTC: ",
          min(precio_es$datetime_utc), " / ", max(precio_es$datetime_utc))

  # 4. JERARQUIA (diap. 27-31): 1,2,3,4,6,8,12,24 horas.
  # Es un ciclo diario de 24 horas UTC, no un dia civil de 23/25 horas.
  # No modela explicitamente la estacionalidad semanal (168 h).
  f <- frequency(precio_ts)
  k <- (1:f)[f %% (1:f) == 0]
  p <- length(k)
  h <- f
  H <- h / k
  Y <- list(thief_aggregate(precio_ts, k))
  holdout <- cfg$mode == "holdout"
  n_train <- length(precio_ts) - if (holdout) h else 0
  y_train <- ts(head(precio_ts, n_train), frequency = f)
  train <- thief_aggregate(y_train, k)
  actual_hourly <- if (holdout) as.numeric(tail(precio_ts, h)) else rep(NA_real_, h)
  target_time <- if (holdout) tail(precio_es$datetime_utc, h) else
    max(precio_es$datetime_utc) + 3600 * seq_len(h)

  # 5. ARIMA POR NIVEL CON PERSISTENCIA (diap. 41).
  # Contexto compatible con la version anterior cuando se usan todos los datos.
  cache_context <- list(input_md5 = unname(tools::md5sum(input_file)),
                        k = as.numeric(k), H = H, aggregation = "mean", holdout = holdout)
  if (!is.null(cfg$last_hours)) cache_context$last_hours <- cfg$last_hours
  if (length(cfg$arima_args)) cache_context$arima_args <- cfg$arima_args
  frc <- thief_load_or_fit(forecast_file, cache_context, y_train, k, H, cfg$arima_args)
  fbase <- as.numeric(unlist(rev(frc), use.names = FALSE))

  # 6. S PARA MEDIAS, W Y G (diap. 30-40).
  S_sum <- tsutils::Sthief(f)
  n_hours <- rowSums(S_sum)
  S <- sweep(S_sum, 1, n_hours, "/")
  nodes <- do.call(rbind, lapply(rev(seq_along(k)), function(j)
    data.frame(k = k[j], bucket = seq_len(H[j]))))
  nodes$node <- paste0("AL", nodes$k, "_", nodes$bucket)
  nodes$start_utc <- format(target_time[(nodes$bucket-1)*nodes$k+1],
                            "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  nodes$end_utc <- format(target_time[nodes$bucket*nodes$k],
                          "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  rownames(S) <- rownames(S_sum) <- nodes$node
  colnames(S) <- colnames(S_sum) <- paste0("h", seq_len(h))
  stopifnot(nrow(S) == sum(H), ncol(S) == h, length(fbase) == nrow(S),
            max(abs(rowSums(S)-1)) < 1e-10)
  actual <- as.numeric(S %*% actual_hourly)
  if (holdout) {
    observed <- unlist(lapply(rev(seq_along(k)), function(j) tail(Y[[1]][[j]], H[j])))
    stopifnot(max(abs(actual - observed)) < 1e-7)
  }
  E <- thief_residual_matrix(attr(frc, "residuals"), k, f)
  if (is.null(E)) warning("La cache no contiene residuos suficientes. Se conservan ",
    "sus forecasts; W por varianzas y MinT no se estiman. Una nueva cache completa los incluira.")
  covariances <- thief_covariances(E, n_hours)
  Ws <- covariances$W
  reconciled <- lapply(Ws, function(W) thief_reconcile(S, W, fbase))
  predictions <- c(list(ARIMA_base = fbase),
                   lapply(reconciled, function(x) x$forecast))
  # Benchmarks: bottom-up y repeticion del ultimo ciclo horario.
  predictions$Bottom_up <- as.numeric(S %*% as.numeric(frc[[1]]))
  predictions$Seasonal_naive <- as.numeric(S %*% as.numeric(tail(y_train, f)))

  # 7. COMBINACION POR REPETICION (diap. 69-72).
  # Ilustracion de pesos iguales; NO es MAPA por estados ni los metodos
  # sparse no publicados de la diap. 82. Puede atenuar la estacionalidad.
  expanded <- vapply(seq_along(k), function(j)
    rep(as.numeric(frc[[j]]), each = k[j]), numeric(h))
  colnames(expanded) <- paste0("AL", k)
  predictions$MTA_equal <- as.numeric(S %*% rowMeans(expanded))

  # 8. MAPA ORIGINAL OPCIONAL (diap. 17-26).
  # Combina estados mediante el paquete del autor. Modelo aditivo para
  # admitir precios negativos; ifh=0 evita el costoso ajuste retrospectivo.
  if (isTRUE(cfg$include_mapa)) {
    mapa_file <- paste0(tools::file_path_sans_ext(forecast_file), "_mapa.rds")
    if (file.exists(mapa_file)) {
      mapa_cache <- readRDS(mapa_file)
      if (!identical(mapa_cache$context, cache_context))
        stop("Cache MAPA incompatible: ", mapa_file)
      mapa_hourly <- mapa_cache$forecast
    } else {
      mapa_hourly <- as.numeric(MAPA::mapa(
        y_train, ppy=f, fh=h, ifh=0, minimumAL=1, maximumAL=f,
        comb="mean", paral=0, display=0, outplot=0,
        hybrid=FALSE, model="AAA", type="ets")$outfor)
      stopifnot(length(mapa_hourly) == h, all(is.finite(mapa_hourly)))
      thief_save_cache(list(context=cache_context, forecast=mapa_hourly), mapa_file)
    }
    stopifnot(length(mapa_hourly) == h, all(is.finite(mapa_hourly)))
    predictions$MAPA_states <- as.numeric(S %*% mapa_hourly)
  }

  # 9. COHERENCIA Y PESOS NEGATIVOS (diap. 35, 38, 42).
  coherence <- data.frame(method = names(predictions),
    max_abs_gap = vapply(predictions, function(z)
      max(abs(z - S %*% tail(z, h))), numeric(1)))
  stopifnot(all(coherence$max_abs_gap[coherence$method != "ARIMA_base"] < 1e-6))
  G <- reconciled$THieF_struct$G
  W <- Ws$THieF_struct
  freco <- matrix(predictions$THieF_struct, ncol=1)
  freco_hourly <- tail(freco, h)
  level_weights <- sapply(rev(k), function(kk)
    rowSums(G[, nodes$k == kk, drop=FALSE]))
  colnames(level_weights) <- paste0("AL", rev(k))
  stopifnot(max(abs(level_weights - 1/p)) < 1e-7)
  # Pesos negativos permitidos: no truncamos precios ni coeficientes.

  # 10. EVALUACION POR NIVEL (diap. 66-67).
  # MAPE no es apropiado con precios cero/negativos. No promediamos
  # metricas entre niveles: las decisiones horarias y diarias son distintas.
  metric_rows <- list()
  if (holdout) for (method in names(predictions)) for (j in seq_along(k)) {
    idx <- which(nodes$k == k[j])
    values <- thief_metrics(actual[idx], predictions[[method]][idx], train[[j]], f/k[j])
    metric_rows[[length(metric_rows)+1]] <- cbind(
      data.frame(method=method, k=k[j], n=length(idx)), as.data.frame(as.list(values)))
  }
  metrics <- if (length(metric_rows)) do.call(rbind, metric_rows) else data.frame()
  if (nrow(metrics)) {
    baseline <- metrics[metrics$method == "ARIMA_base", ]
    metrics$MAE_gain_pct <- 100 * (1-metrics$MAE / baseline$MAE[match(metrics$k, baseline$k)])
    metrics$MAE_gain_pct[!is.finite(metrics$MAE_gain_pct)] <- NA_real_
    message("Evaluacion horaria del bloque reservado (no prueba de superioridad general):")
    print(metrics[metrics$k == 1, ], row.names=FALSE)
  }

  # 11. BACKTEST OPCIONAL + DISTRIBUCIONES EMPIRICAS (diap. 65-67).
  # Origenes diarios anteriores al test final. Cada modelo ve solo su pasado.
  # Se cachea cada origen; sus errores NO se usan para elegir el ganador
  # del test final. Los intervalos son marginales; las trayectorias son coherentes.
  calibration_errors <- list()
  calibration_metrics <- list()
  intervals <- data.frame()
  scenarios <- NULL
  if (cfg$calibration_days > 0) {
    origins <- n_train - f * rev(seq_len(cfg$calibration_days))
    if (min(origins) < f*10) stop("No hay suficiente historia para calibration_days.")
    cal_dir <- paste0(tools::file_path_sans_ext(forecast_file), "_calibration")
    for (origin in origins) {
      cal_context <- cache_context
      cal_context$origin <- origin
      cal_y <- ts(head(y_train, origin), frequency=f)
      cal_frc <- thief_load_or_fit(file.path(cal_dir, paste0("origin_",origin,".rds")),
        cal_context, cal_y, k, H, cfg$arima_args)
      cal_base <- as.numeric(unlist(rev(cal_frc), use.names=FALSE))
      # W estructural no usa informacion posterior al origen.
      cal_pred <- as.numeric(S %*% G %*% cal_base)
      cal_actual <- as.numeric(S %*% y_train[origin+seq_len(h)])
      calibration_errors[[length(calibration_errors)+1]] <- cal_actual-cal_pred
      cal_train <- thief_aggregate(cal_y, k)
      for (method in c("ARIMA_base","THieF_struct")) for (j in seq_along(k)) {
        idx <- which(nodes$k == k[j])
        values <- thief_metrics(cal_actual[idx],
          if (method == "ARIMA_base") cal_base[idx] else cal_pred[idx],
          cal_train[[j]], f/k[j])
        calibration_metrics[[length(calibration_metrics)+1]] <- cbind(
          data.frame(origin=origin, method=method, k=k[j]),
          as.data.frame(as.list(values)))
      }
    }
    if (length(calibration_errors) >= 20) {
      errors <- do.call(rbind, calibration_errors)
      # Muestreo de vectores diarios enteros preserva dependencia intradia.
      # Errores sin centrar: se conserva el sesgo empirico de calibracion.
      set.seed(cfg$seed)
      sampled <- sample.int(nrow(errors), cfg$simulations, replace=TRUE)
      hourly_errors <- errors[sampled, (ncol(errors)-h+1):ncol(errors), drop=FALSE]
      hourly_paths <- sweep(hourly_errors, 2, as.numeric(freco_hourly), "+")
      scenarios <- S %*% t(hourly_paths)
      stopifnot(max(abs(scenarios - S %*% tail(scenarios, h))) < 1e-6)
      q <- t(apply(scenarios, 1, stats::quantile, probs=c(.025,.1,.5,.9,.975)))
      intervals <- cbind(nodes, data.frame(lo95=q[,1],lo80=q[,2],
                          median=q[,3],hi80=q[,4],hi95=q[,5]))
      if (holdout) {
        intervals$covered80 <- actual >= intervals$lo80 & actual <= intervals$hi80
        intervals$covered95 <- actual >= intervals$lo95 & actual <= intervals$hi95
      }
    } else message("Intervalos omitidos: se requieren al menos 20 dias de calibracion.")
  }
  backtest <- if (length(calibration_metrics)) do.call(rbind, calibration_metrics) else data.frame()
  backtest_summary <- data.frame()
  if (nrow(backtest)) {
    backtest_summary <- aggregate(cbind(ME,MAE,MASE,RMSSE) ~ method+k,
                                   data=backtest, FUN=mean)
    rmse_summary <- aggregate(I(RMSE^2) ~ method+k, data=backtest, FUN=mean)
    names(rmse_summary)[3] <- "MSE"
    backtest_summary <- merge(backtest_summary,rmse_summary,by=c("method","k"))
    backtest_summary$RMSE <- sqrt(backtest_summary$MSE)
    backtest_summary$days <- cfg$calibration_days
  }

  # 12. GRAFICOS: agregacion, comparacion por nivel y pesos.
  thief_png(file.path(output_dir, "01_agregaciones.png"), 1600, 1600, function() {
    par(mfrow=c(4,2), mar=c(3,4,3,1))
    for (j in seq_along(k)) {
      yy <- tail(Y[[1]][[j]], max(14*f/k[j], H[j]+1))
      plot(as.numeric(yy), type="l", col="#333333", xlab="Bloque temporal",
           ylab="Precio medio", main=paste("Agregacion:", k[j], "horas"))
    }
  })
  thief_png(file.path(output_dir, "02_arima_vs_thief.png"), 1600, 1600, function() {
    par(mfrow=c(4,2), mar=c(3,4,3,1))
    for (j in seq_along(k)) {
      idx <- which(nodes$k == k[j])
      history <- tail(as.numeric(train[[j]]), max(7*f/k[j], 2))
      nh <- length(history)
      xx <- nh + seq_len(H[j])
      yy <- c(history, actual[idx], fbase[idx], freco[idx])
      plot(seq_len(nh), history, type="l", xlim=c(1,nh+H[j]),
           ylim=range(yy, finite=TRUE), xlab="Bloque temporal", ylab="Precio medio",
           main=paste("k =", k[j], "horas"))
      abline(v=nh+.5, lty=3, col="grey50")
      if (holdout) lines(c(nh,xx),c(tail(history,1),actual[idx]),col="black",lwd=2)
      lines(c(nh,xx),c(tail(history,1),fbase[idx]),col="#c74635",lwd=2,type="o",pch=16,cex=.5)
      lines(c(nh,xx),c(tail(history,1),freco[idx]),col="#2364aa",lwd=2,type="o",pch=16,cex=.5)
      legend("topleft", c("Observado","ARIMA","THieF estructural"),
             col=c("black","#c74635","#2364aa"),lty=1,cex=.7,bty="n")
    }
  })
  thief_png(file.path(output_dir, "03_pesos_G.png"), 1600, 700, function() {
    par(mar=c(5,4,3,2))
    bound <- max(abs(G))
    image(seq_len(ncol(G)), seq_len(nrow(G)), t(G),
          col=colorRampPalette(c("#2166ac","white","#b2182b"))(101),
          zlim=c(-bound,bound), xlab="Nodo base (AL24 ... AL1)",
          ylab="Hora reconciliada", main="G estructural: azul negativo, rojo positivo")
    abline(v=cumsum(rev(H))+.5, col="grey50", lty=3)
  })
  if (holdout) thief_png(file.path(output_dir, "04_error_por_nivel.png"), 1400, 700, function() {
    keep <- c("ARIMA_base", names(reconciled))
    error_matrix <- vapply(keep, function(m) {
      rows <- metrics[metrics$method == m, ]
      rows$MAE[match(k, rows$k)]
    }, numeric(p))
    # Cada fila es un nivel. Las series se dibujan sin promediar niveles.
    matplot(k, error_matrix, type="b", pch=seq_along(keep), lty=1,
            col=seq_along(keep), xlab="Horas agregadas (k)", ylab="MAE",
            main="Error en el bloque reservado por nivel")
    legend("topright",keep,col=seq_along(keep),pch=seq_along(keep),lty=1,cex=.8,bty="n")
  })
  if (nrow(intervals)) thief_png(file.path(output_dir,"05_intervalos_horarios.png"),1400,700,function() {
    hourly <- tail(intervals,h)
    plot(seq_len(h),as.numeric(freco_hourly),type="n",
         ylim=range(c(hourly$lo95,hourly$hi95,actual_hourly,freco_hourly),finite=TRUE),
         xlab="Hora del horizonte",ylab="Precio",
         main="THieF estructural: intervalos empiricos (calibracion retrospectiva)")
    polygon(c(seq_len(h),rev(seq_len(h))),c(hourly$lo95,rev(hourly$hi95)),
            col="#dce9f5",border=NA)
    polygon(c(seq_len(h),rev(seq_len(h))),c(hourly$lo80,rev(hourly$hi80)),
            col="#a5c8e8",border=NA)
    lines(seq_len(h),as.numeric(freco_hourly),col="#2364aa",lwd=2)
    if (holdout) lines(seq_len(h),actual_hourly,col="black",lwd=2)
    legend("topleft",c("THieF","Observado","80%","95%"),
           col=c("#2364aa","black","#a5c8e8","#dce9f5"),lty=1,lwd=c(2,2,8,8),bty="n")
  })

  # 13. EXPORTACION REPRODUCIBLE.
  forecasts <- cbind(nodes, actual=actual, as.data.frame(predictions))
  write.csv(forecasts, file.path(output_dir,"forecasts_por_nivel.csv"),row.names=FALSE)
  write.csv(coherence, file.path(output_dir,"coherencia.csv"),row.names=FALSE)
  if (nrow(metrics)) write.csv(metrics,file.path(output_dir,"metricas_por_nivel.csv"),row.names=FALSE)
  if (nrow(backtest)) write.csv(backtest,file.path(output_dir,"backtest.csv"),row.names=FALSE)
  if (nrow(backtest_summary)) write.csv(backtest_summary,file.path(output_dir,"backtest_resumen.csv"),row.names=FALSE)
  if (nrow(intervals)) write.csv(intervals,file.path(output_dir,"intervalos_empiricos.csv"),row.names=FALSE)
  if (!is.null(attr(frc,"models")))
    write.csv(attr(frc,"models"),file.path(output_dir,"modelos_arima.csv"),row.names=FALSE)
  write.csv(S,file.path(output_dir,"S_medias.csv"))
  for (method in names(reconciled)) {
    write.csv(Ws[[method]],file.path(output_dir,paste0("W_",method,".csv")))
    write.csv(reconciled[[method]]$G,file.path(output_dir,paste0("G_",method,".csv")))
  }
  write.csv(level_weights,file.path(output_dir,"pesos_por_nivel.csv"))
  capture.output(sessionInfo(),file=file.path(output_dir,"sessionInfo.txt"))
  notes <- c(
    "THieF - ESIOS. Referencia: Kourentzes, ISF 2025.",
    paste("Modo:", cfg$mode, "| Observaciones:", length(precio_ts)),
    paste("Cache:", forecast_file),
    "Medias de bloques horarios UTC anclados al final; no dias civiles locales.",
    "S: 60 x 24. Orden AL24, AL12, AL8, AL6, AL4, AL3, AL2, AL1.",
    "Diap. 27-42: agregacion, ARIMA independiente, S/W/G, coherencia y pesos.",
    "Diap. 36-40: estructural, varianzas y MinT; OLS como referencia adicional.",
    paste("Filas completas de residuos para W:", if(is.null(E)) 0 else nrow(E)),
    paste("MinT lambda:", covariances$lambda, "| regularizacion:", covariances$ridge),
    "W estimada con residuos in-sample un paso: aproximacion, no covarianza conocida a 24 pasos.",
    "Diap. 66-67: metricas separadas por nivel; no media global de la jerarquia.",
    "Un solo bloque de test no demuestra superioridad general.",
    paste("Diap. 17-26: MAPA por estados activado:", cfg$include_mapa),
    "Diap. 69-72: MTA_equal ilustra repeticion + media, no el algoritmo MAPA.",
    paste("Diap. 65: dias de calibracion retrospectiva:", cfg$calibration_days),
    "Intervalos: bootstrap empirico de dias completos; sin garantia de cobertura futura.",
    "Los escenarios son coherentes; los cuantiles marginales no tienen por que serlo.",
    "Diap. 3-16, 43-50: motivacion/resultados publicados; no se reproducen sus estudios.",
    "Diap. 51-64: demanda intermitente, ciclos de vida, cross-temporal, exogenas y ML",
    "requieren modelos/datos/disenos especificos y no son pasos obligatorios para esta serie.",
    "Diap. 73-82: resultados teoricos y propuestas no publicadas; no se atribuye al",
    "script una replica de simulaciones o algoritmos sparse sin especificacion completa.",
    "Referencia complementaria: https://otexts.com/fpp3/reconciliation.html",
    "Intensidad de shrinkage: https://github.com/earowang/hts/blob/master/R/MinT.R"
  )
  writeLines(notes,file.path(output_dir,"LEEME_resultados.txt"),useBytes=TRUE)
  results <- list(config=cfg, precio_es=precio_es, precio_ts=precio_ts, f=f,p=p,h=h,k=k, H=H, Y=Y,
    frc=frc, S_sum=S_sum, S=S, W=W, G=G, W_all=Ws, reconciled=reconciled,
    fbase=matrix(fbase,ncol=1), freco=freco, freco_hourly=freco_hourly,
    residual_matrix=E, shrinkage=covariances[c("lambda","ridge")],
    forecasts=forecasts, metrics=metrics, coherence=coherence,
    backtest=backtest, backtest_summary=backtest_summary,intervals=intervals, scenarios=scenarios)
  saveRDS(results,file.path(output_dir,"resultados_thief.rds"))
  message("Resultados guardados en: ", output_dir)
  invisible(results)
}

# 14. EJECUCION Y OBJETOS DISPONIBLES EN LA CONSOLA
if (isTRUE(getOption("thief.autorun", TRUE))) {
  resultados_thief <- run_thief_esios(getOption("thief.config", list()))
  # Mantiene los nombres usados en el taller para inspeccion interactiva.
  list2env(resultados_thief[c("precio_es","precio_ts","f","p","h","k","H","Y","frc",
    "S_sum","S","W","G","fbase","freco","freco_hourly")], envir=environment())
}
