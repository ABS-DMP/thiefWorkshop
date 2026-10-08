# Ejecutar desde la raiz del proyecto: source("tests/test_thief_esios.R")
# Pruebas con los datos locales, sin modificar data/raw ni la cache principal.
local({
  env <- new.env(parent=globalenv())
  previous <- options(thief.autorun=FALSE)
  on.exit(options(previous),add=TRUE)
  sys.source("thief_esios.R",envir=env)
  cfg <- list(project_dir=normalizePath("."),last_hours=24*45,
    forecast_file="data/processed/tests/frc_smoke.rds",
    output_dir="output/tests/verification",include_mapa=TRUE,
    calibration_days=20L,simulations=100L,
    arima_args=list(max.p=1,max.q=1,max.P=0,max.Q=0,
                    stepwise=TRUE,approximation=TRUE))
  x <- env$run_thief_esios(cfg)
  stopifnot(identical(dim(x$S),c(60L,24L)),nrow(x$metrics)==72,
    nrow(x$backtest)==320,nrow(x$backtest_summary)==16,
    nrow(x$intervals)==60,identical(dim(x$scenarios),c(60L,100L)),
    all(x$coherence$max_abs_gap[x$coherence$method!="ARIMA_base"]<1e-6),
    max(abs(x$G%*%x$S-diag(24)))<1e-8,
    max(abs(x$scenarios-x$S%*%tail(x$scenarios,24)))<1e-6,
    all(x$intervals$lo95<=x$intervals$lo80),
    all(x$intervals$lo80<=x$intervals$median),
    all(x$intervals$median<=x$intervals$hi80),
    all(x$intervals$hi80<=x$intervals$hi95))
  # Equivalencia con agregacion mediante sumas.
  D <- diag(rowSums(x$S_sum))
  sum_reco <- env$thief_reconcile(x$S_sum,D,D%*%x$fbase)$forecast
  stopifnot(max(abs(sum_reco/diag(D)-x$freco))<1e-6)
  # Cada bloque de test equivale a la media horaria observada.
  observed <- as.numeric(x$S%*%tail(x$precio_ts,24))
  stopifnot(max(abs(observed-x$forecasts$actual))<1e-8)
  # No debe entrenar al repetir, tampoco en los origenes de calibracion.
  fit_original <- env$thief_fit
  env$thief_fit <- function(...) stop("Entrenamiento inesperado")
  stamp <- file.info(cfg$forecast_file)$mtime
  y <- env$run_thief_esios(cfg)
  stopifnot(identical(x$frc,y$frc),
    identical(stamp,file.info(cfg$forecast_file)$mtime),
    identical(x$scenarios,y$scenarios))
  changed <- cfg
  changed$last_hours <- 24*44
  problem <- tryCatch({env$run_thief_esios(changed);""},
                      error=function(e) conditionMessage(e))
  stopifnot(grepl("otros datos/configuracion",problem,fixed=TRUE))
  # Historia no divisible por 24 y pronostico fuera de muestra.
  env$thief_fit <- fit_original
  future <- cfg
  future$mode <- "future"
  future$last_hours <- 24*45-1
  future$calibration_days <- 0L
  future$include_mapa <- FALSE
  future$forecast_file <- "data/processed/tests/frc_future.rds"
  future$output_dir <- "output/tests/future"
  z <- env$run_thief_esios(future)
  stopifnot(nrow(z$metrics)==0,all(is.na(z$forecasts$actual)),
    all(as.POSIXct(z$forecasts$start_utc,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC") >
          max(z$precio_es$datetime_utc)))
  # Residuos degenerados: regularizacion de W y solucion estable.
  ws <- env$thief_covariances(matrix(0,30,60),rowSums(x$S_sum))
  constant <- env$thief_reconcile(x$S,ws$W$THieF_MinT,rep(50,60))
  stopifnot(max(abs(constant$forecast-50))<1e-6)
  cat("PASS: modelos, MAPA, S/W/G, cache, fechas, backtest e intervalos.\n")
})

