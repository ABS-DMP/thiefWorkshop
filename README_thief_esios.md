# THieF sobre los precios ESIOS locales

El programa lee `data/raw/precio_2020_2026.csv` dentro de este proyecto.
No descarga datos. Los otros CSV de demanda, generacion y gas no entran en
el modelo univariante: utilizarlos como regresores requiere decidir sus
retardos y conocer su disponibilidad en cada origen de prediccion.

## Ejecucion habitual

En R:

```r
options(thief.autorun = TRUE, thief.config = list())
source("C:/Users/aleja/Documents/source/repos/thiefWorkshop/thief_esios.R")
```

Se usa todo el historico y se reservan las ultimas 24 observaciones para
evaluacion. El horizonte son 24 horas UTC consecutivas, ancladas al final
de los datos; no implica que cada bloque corresponda a un dia civil local.

La primera ejecucion entrena ARIMA en los ocho niveles. Las siguientes
cargan `data/processed/frc_arima_base.rds`. Una cache incompatible detiene
la ejecucion: se debe elegir otra ruta de cache o archivar la anterior.
Las caches nuevas incluyen residuos y descripcion de los modelos.
Una cache antigua sin residuos permite reconciliacion estructural y OLS,
pero no permite reconstruir W por varianzas/MinT sin nuevos ajustes.

Resultados en `output/thief_esios/`: forecasts con fechas por nivel,
metricas, coherencia, matrices S/W/G, modelos, graficos, notas de la
ejecucion y `resultados_thief.rds`. El objeto `resultados_thief` y los
objetos principales del taller quedan disponibles en la consola.
Los CSV de resultados opcionales solo se escriben cuando se activa la
opcion correspondiente; para comparar configuraciones use carpetas de
salida diferentes. El RDS y LEEME_resultados.txt describen la ultima ejecucion.

## Correspondencia con el PDF del taller

| Diapositivas | Implementacion |
|---|---|
| 27-31 | Divisores de 24, medias no solapadas, alineacion al final, S de 60 x 24. |
| 33-35 | Forecasts independientes, G, reconciliacion y comprobacion G S = I. |
| 36-40 | W estructural para medias, W por varianzas, MinT shrinkage; OLS adicional. Residuos agrupados en jerarquias diarias completas. |
| 38, 42 | Pesos por nivel, mapa de G, pesos negativos permitidos y prueba de coherencia. |
| 41 | Graficos de observaciones, ARIMA base y THieF por nivel. |
| 17-26 | Comparacion MAPA por estados, opcional, mediante el paquete del autor. Modelo AAA aditivo para precios que pueden ser negativos. |
| 65 | Muestreo de vectores diarios de errores retrospectivos, con trayectorias coherentes e intervalos marginales 80/95%. Opcional. |
| 66-67 | MAE, RMSE, ME, MASE y RMSSE por nivel; backtest diario opcional. |
| 69-72 | Ilustracion MTA_equal: repeticion de forecasts de cada nivel y media simple. No se presenta como MAPA por estados. |

El PDF tambien incluye motivacion teorica, resultados de estudios,
aplicaciones a demanda intermitente/ciclos de vida, jerarquias entre series,
variables exogenas, aprendizaje automatico y propuestas no publicadas de
jerarquias dispersas (diapositivas 73-82). No define un unico algoritmo
ejecutable que aplique todas esas extensiones a una serie de precios.
El script identifica ese alcance; no reproduce estudios ni inventa
implementaciones para resultados sin especificacion completa.

## Activar MAPA y la evaluacion retrospectiva

```r
options(thief.config = list(
  include_mapa = TRUE,
  calibration_days = 30L,
  output_dir = "output/thief_esios_extended"
))
source("C:/Users/aleja/Documents/source/repos/thiefWorkshop/thief_esios.R")
```

Cada origen del backtest utiliza solo datos anteriores a su horizonte.
Se guardan forecasts por origen en una carpeta de calibracion junto a
la cache base, por lo que repetir la ejecucion no vuelve a entrenarlos.
El backtest compara ARIMA con THieF estructural. El test final compara
tambien las otras variantes. Las metricas del backtest de cada nivel
no se promedian con las de otros niveles.

Con al menos 20 dias de calibracion se generan intervalos empiricos y
escenarios conjuntos. Muestrear dias completos conserva la dependencia
entre las horas. Su cobertura futura depende de que esos errores sean
representativos. Los cuantiles marginales no tienen que ser coherentes,
aunque cada escenario si lo es. W por varianzas y MinT utiliza residuos
in-sample de un paso: es una aproximacion a los errores del horizonte.

## Predecir las siguientes 24 horas, en vez de reservar un test

```r
options(thief.config = list(mode = "future"))
source("C:/Users/aleja/Documents/source/repos/thiefWorkshop/thief_esios.R")
```

Este modo entrena con todos los valores disponibles y utiliza otra cache,
`data/processed/frc_arima_base_future.rds`, y otra carpeta de resultados.
No calcula metricas frente a observaciones futuras desconocidas.

## Comprobaciones realizadas

`tests/test_thief_esios.R` prueba el flujo completo con 45 dias de datos
reales locales y una busqueda ARIMA reducida: ocho niveles, las variantes
de W, MAPA, backtest de 20 dias, intervalos, cache sin reentrenamiento,
rechazo de cache incompatible, equivalencia entre sumas/medias, historia
no divisible por 24, modo futuro y residuos degenerados.

Se ejecuta desde la raiz del proyecto:

```r
source("tests/test_thief_esios.R")
```

Las salidas de estas pruebas estan en `output/tests/` y sus caches en
`data/processed/tests/`. Son comprobaciones del funcionamiento, no los
resultados definitivos del TFM. No se ha ejecutado aqui el entrenamiento
completo con los parametros ARIMA por defecto.

Fuentes metodologicas complementarias:
- https://otexts.com/fpp3/reconciliation.html
- https://github.com/earowang/hts/blob/master/R/MinT.R

