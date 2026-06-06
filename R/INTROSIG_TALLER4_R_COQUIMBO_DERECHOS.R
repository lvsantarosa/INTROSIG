# ==============================================================================
# EVALUACIÓN PRÁCTICA: VULNERABILIDAD HÍDRICA Y ESTRÉS VEGETATIVO (COQUIMBO)
# ==============================================================================
# # Institución: Pontificia Universidad Católica de Valparaíso (PUCV)
# Prof. Lucas Vituri Santarosa
#
# Descripción: Diagnóstico espacial para evaluar la dependencia de aguas 
# subterráneas de la matriz frutícola, utilizando datos de la DGA y 
# mosaicos de teledetección estática (Landsat 8 L2SP - Sequía estival).
# ==============================================================================

# Cargar librerías necesarias
library(terra)      # Procesamiento espacial de vectores y rasters
library(dplyr)      # Manipulación de tablas de datos
library(ggplot2)    # Visualización y cartografía
library(tidyterra)  # Integración de SpatRaster/SpatVector con ggplot2
library(readxl)     # Lectura de archivos Excel de la DGA
library(janitor)    # Limpieza y estandarización de nombres de columnas

# ------------------------------------------------------------------------------
# FASE I: ANÁLISIS VECTORIAL Y PRESIÓN SUBTERRÁNEA
# ------------------------------------------------------------------------------
# Configurar el directorio de trabajo (Ajustar según corresponda al alumno)
setwd('~/03_AULAS/01_POS/01_CLASES/2026/01_INTROSIG/02_TALLER/05_TALLER3_VECTOR_RASTER_R/')

# 1. Cargar y estandarizar datos provinciales
prov <- vect("DATOS/04_REGIONES_PROVINCIAS/01_PROVINCIAS/Provincias.shp") %>% 
  project("EPSG:32719")
prov_coquimbo <- prov[prov$Region == "Región de Coquimbo", ]
prov_coquimbo$Area_km2 <- terra::expanse(prov_coquimbo, unit = "km")

# 2. Cargar Derechos de Agua Subterránea (DGA)
dga_excel <- read_excel("DATOS/Derechos_Concedidos_IV_Region.xls", skip = 6) %>%
  clean_names() %>%
  filter(naturaleza_del_agua == "Subterranea") %>%
  filter(!is.na(utm_este_captacion_m) & !is.na(utm_norte_captacion_m))

# Crear vector de puntos explícito
derechos <- vect(dga_excel, 
                 geom = c("utm_este_captacion_m", "utm_norte_captacion_m"), 
                 crs = "EPSG:32719")

# 3. Intersección espacial: Calcular Densidad de Derechos por Provincia
derechos_prov <- terra::intersect(derechos, prov_coquimbo)

tabla_derechos <- as.data.frame(derechos_prov) %>%
  group_by(Provincia) %>%
  summarise(
    Total_Derechos = n(),
    Area_Provincial_km2 = first(Area_km2)
  ) %>%
  mutate(Densidad_Derechos_km2 = Total_Derechos / Area_Provincial_km2) %>%
  arrange(desc(Densidad_Derechos_km2))

print("TABLA 1: Presión sobre el Acuífero - Densidad de Derechos por Provincia")
print(tabla_derechos)

# 4. Cuantificar Superficie Frutícola y Seleccionar Área Prioritaria
frutales <- vect("DATOS/01_FRUTALES/4__cubierta_frutal_de_pol_gonos_regi_n_de_coquimbo_a_o_2024/4__cubierta_frutal_de_pol_gonos_regi_n_de_coquimbo_a_o_2024.shp") %>% 
  project("EPSG:32719")

frutales_prov <- terra::intersect(frutales, prov_coquimbo)
frutales_prov$Area_ha <- terra::expanse(frutales_prov, unit = "ha")

tabla_fruticola <- as.data.frame(frutales_prov) %>%
  group_by(Provincia) %>%
  summarise(Total_Frutales_ha = sum(Area_ha, na.rm = TRUE)) %>%
  arrange(desc(Total_Frutales_ha))

# Identificamos programáticamente la provincia prioritaria (Mayor Área Agrícola)
prov_prioritaria_nome <- tabla_fruticola$Provincia[1]

# Subconjuntos espaciales definitivos para la Fase II
prov_prioritaria <- prov_coquimbo[prov_coquimbo$Provincia == prov_prioritaria_nome, ]
frutales_prioridad <- frutales_prov[frutales_prov$Provincia == prov_prioritaria_nome, ]

# ------------------------------------------------------------------------------
# FASE II: ANÁLISIS RÁSTER (CONSTRUCCIÓN DE MOSAICOS Y CÁLCULO DE NDMI)
# ------------------------------------------------------------------------------
# Fórmula Teórica NDMI: (NIR - SWIR) / (NIR + SWIR)
# En Landsat 8/9: NIR = Banda 5, SWIR = Banda 6

# 1. Función para procesar y crear mosaicos con PROGRAMACIÓN DEFENSIVA (tryCatch)
procesar_y_mosaicar <- function(directorio, patron, poligono) {
  archivos <- list.files(directorio, pattern = patron, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  
  if(length(archivos) == 0) {
    stop(paste("ERROR: No se encontraron archivos para el patrón:", patron))
  }
  
  print(sprintf("Se encontraron %d escenas. Iniciando pre-procesamiento...", length(archivos)))
  
  lista_rasters <- list()
  
  for(i in seq_along(archivos)) {
    # PASO A: Lectura SEGURA de la escena original usando tryCatch
    # Evita que un archivo TIF corrupto (GDAL Error 4) detenga todo el script
    r_raw <- tryCatch({
      rast(archivos[i])
    }, error = function(e) {
      warning(sprintf("\n¡ATENCIÓN! Archivo corrupto o ilegible omitido: %s\n", basename(archivos[i])))
      return(NULL) # Si falla, devuelve NULL
    })
    
    # Si el raster es NULL (hubo error), saltamos al siguiente archivo del loop
    if(is.null(r_raw)) next 
    
    # PASO B: Homogeneización del sistema de coordenadas
    r_proj <- project(r_raw, crs(poligono))
    
    # PASO C: Recorte inmediato (Crop) a la extensión (Bounding Box)
    r_crop <- crop(r_proj, ext(poligono))
    
    lista_rasters[[length(lista_rasters) + 1]] <- r_crop
  }
  
  # Verificación por si todos los archivos estaban corruptos
  if(length(lista_rasters) == 0) stop("Todos los archivos encontrados estaban corruptos.")
  
  # PASO D: Ensamblaje del mosaico (Unión espacial)
  if(length(lista_rasters) > 1) {
    coleccion <- sprc(lista_rasters) # Spatial Raster Collection
    mosaico <- terra::mosaic(coleccion, fun = "mean")
  } else {
    mosaico <- lista_rasters[[1]]
  }
  
  # PASO E: Máscara final ajustada a la geometría de la provincia
  mosaico_final <- mask(mosaico, poligono)
  
  return(mosaico_final)
}

# 2. Generación de los mosaicos por banda
print("Construyendo mosaico para la Banda 5 (NIR)...")
b5_mosaic <- procesar_y_mosaicar("DATOS/Landsat_Coquimbo", "B5\\.tif$", prov_prioritaria)

print("Construyendo mosaico para la Banda 6 (SWIR)...")
b6_mosaic <- procesar_y_mosaicar("DATOS/Landsat_Coquimbo", "B6\\.tif$", prov_prioritaria)

# 3. Calibración Radiométrica (Factor de Escala Landsat Collection 2 Nivel 2)
# Reflectancia = (DN * 0.0000275) - 0.2
b5_ref <- clamp((b5_mosaic * 0.0000275) - 0.2, lower = 0, upper = 1)
b6_ref <- clamp((b6_mosaic * 0.0000275) - 0.2, lower = 0, upper = 1)

# 4. Cálculo del NDMI (Índice de Humedad de Diferencia Normalizada)
print("Calculando NDMI del área de estudio...")
ndmi <- (b5_ref - b6_ref) / (b5_ref + b6_ref)
names(ndmi) <- "NDMI_Marzo"

# ------------------------------------------------------------------------------
# FASE III: INTEGRACIÓN VECTOR/RÁSTER Y ESTADÍSTICAS ZONALES
# ------------------------------------------------------------------------------

# 1. Extracción de métricas de humedad
print("Extrayendo estadísticas zonales (Media y Mínimo)...")
ext_mean <- terra::extract(ndmi, frutales_prioridad, fun = mean, na.rm = TRUE)
ext_min  <- terra::extract(ndmi, frutales_prioridad, fun = min, na.rm = TRUE)

frutales_prioridad$NDMI_Medio <- ext_mean$NDMI_Marzo
frutales_prioridad$NDMI_Minimo <- ext_min$NDMI_Marzo

# Filtrar polígonos inválidos (menores que el tamaño del píxel)
df_resultados <- as.data.frame(frutales_prioridad) %>% 
  filter(!is.na(NDMI_Medio))

# 2. Ranking de Resiliencia Hídrica por Especie
ranking_resiliencia <- df_resultados %>%
  group_by(ESPECIE) %>%
  summarise(
    Humedad_Promedio_NDMI = mean(NDMI_Medio, na.rm = TRUE),
    NDMI_Mas_Critico = min(NDMI_Minimo, na.rm = TRUE),
    Area_Evaluada_ha = sum(Area_ha, na.rm = TRUE)
  ) %>%
  filter(Area_Evaluada_ha > 10) %>% # Excluir cultivos residuales/anecdóticos
  arrange(desc(Humedad_Promedio_NDMI))

print("TABLA 2: Ranking de Resiliencia Frutícola frente a Sequía (Marzo)")
print(ranking_resiliencia)

# 3. Mapa de Anomalías y Estrés Severo
# Clasificación categórica: Huertos bajo el umbral crítico de humedad (NDMI < 0)
frutales_prioridad$Estado_Hidrico <- ifelse(frutales_prioridad$NDMI_Medio < 0, 
                                            "Estrés Hídrico Severo (< 0)", 
                                            "Humedad de Dosel Estable (>= 0)")

# Generación cartográfica
mapa_anomalias <- ggplot() +
  geom_spatvector(data = prov_prioritaria, fill = "gray95", color = "black", linewidth = 0.5) +
  geom_spatvector(data = frutales_prioridad, aes(fill = Estado_Hidrico), color = NA) +
  scale_fill_manual(values = c("Estrés Hídrico Severo (< 0)" = "firebrick", 
                               "Humedad de Dosel Estable (>= 0)" = "forestgreen")) +
  labs(
    title = paste("Anomalías de Humedad Frutícola - Provincia de", prov_prioritaria_nome),
    subtitle = "Identificación de huertos con estrés severo mediante NDMI (Landsat 8)",
    fill = "Estado Fisiológico",
    caption = "Fuente: Catastro Frutícola ODEPA e Imágenes Landsat L2SP"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(mapa_anomalias)

# ==============================================================================
# FIN DEL DIAGNÓSTICO
# ==============================================================================
