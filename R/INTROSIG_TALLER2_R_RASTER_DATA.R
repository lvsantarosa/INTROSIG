# ==============================================================================
# SCRIPT INTEGRACIÓN VECTOR-RASTER Y CÁLCULO DE ÍNDICES ESPECTRALES
# ==============================================================================
# Autor: Prof. Lucas Vituri Santarosa
# Laboratório de Recursos Hídricos y Geotecnologias
# Institución: Pontificia Universidad Católica de Valparaíso (PUCV)
#
# Descripción: Este script toma como base los datos vectoriales provinciales y 
# de cultivos frutales, y los integra con imágenes satelitales Landsat 8 (Nivel 2). 
# El flujo metodológico incluye el recorte a la Provincia de Quillota, la 
# conversión de valores de píxel a reflectancia de superficie, el cálculo 
# de índices exploratorios (NDVI, NDWI, NDMI) y la extracción de estadísticas 
# zonales para comparar el estado espectral de las diferentes especies frutales.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# PASO 1: CARGA DE LIBRERÍAS
# ------------------------------------------------------------------------------
library(terra)      # Procesamiento espacial de raster y vectores
library(dplyr)      # Manipulación de datos tabulares
library(ggplot2)    # Visualización gráfica
library(tidyr)      # Reestructuración de datos para gráficos
library(tidyterra)

# ------------------------------------------------------------------------------
# PASO 2: PREPARACIÓN DE DATOS VECTORIALES (RECORTE A QUILLOTA)
# ------------------------------------------------------------------------------
# Como buena práctica, no cargaremos toda la región si solo analizaremos Quillota.
# Cargamos el shapefile de la cubierta frutal.
# Nota: Ajusta la ruta según dónde guardaste tus archivos.
frut <- vect("~/03_AULAS/01_POS/01_CLASES/2026/01_INTROSIG/02_TALLER/02_TALLER1_VECTOR/00_DATOS/00_VECTOR/01_FRUTALES/5__cubierta_frutal_de_pol_gonos_regi_n_de_valpara_so_2025/5__cubierta_frutal_de_pol_gonos_regi_n_de_valpara_so_2025.shp")
prov <- vect('~/03_AULAS/01_POS/01_CLASES/2026/01_INTROSIG/02_TALLER/02_TALLER1_VECTOR/00_DATOS/00_VECTOR/04_REGIONES_PROVINCIAS/01_PROVINCIAS/Provincias.shp')

# 2.1 Cargar y filtrar la Provincia de Quillota
quillota_prov <- prov %>%
  filter(Provincia == "Quillota") %>% # Filtro explícito para la zona de estudio
  project("EPSG:32719")               # Proyección UTM 19S WGS84

# 2.2 Cargar la cobertura frutal y recortarla (clip) con el límite de Quillota
frut_proj <- project(frut, "EPSG:32719")

# Usamos 'intersect' para mantener solo los huertos dentro de Quillota
frut_quillota <- terra::intersect(quillota_prov, frut_proj) 

# ------------------------------------------------------------------------------
# PASO 3: PREPARACIÓN DE DATOS RASTER (LANDSAT 8 L2SP)
# ------------------------------------------------------------------------------
# Imagen base: LC08_L2SP_233083_20250718_20250726_02_T1
# Identificamos las bandas necesarias para vegetación y agua:
# Banda 3 (Green), Banda 4 (Red), Banda 5 (NIR), Banda 6 (SWIR 1)

# Asumimos que los archivos TIF están en una carpeta de trabajo
ruta_img <- "C:/Users/lucas/Documents/03_AULAS/01_POS/01_CLASES/2026/01_INTROSIG/02_TALLER/03_TALLER2_RASTER/LC08_L2SP_233083_20250718_20250726_02_T1/"

# 3.1 Cargar bandas individuales a la memoria temporal
b3_raw <- rast(paste0(ruta_img, "LC08_L2SP_233083_20250718_20250726_02_T1_SR_B3.tif")) %>% project("EPSG:32719") 
b4_raw <- rast(paste0(ruta_img, "LC08_L2SP_233083_20250718_20250726_02_T1_SR_B4.tif")) %>% project("EPSG:32719")
b5_raw <- rast(paste0(ruta_img, "LC08_L2SP_233083_20250718_20250726_02_T1_SR_B5.tif")) %>% project("EPSG:32719")
b6_raw <- rast(paste0(ruta_img, "LC08_L2SP_233083_20250718_20250726_02_T1_SR_B6.tif")) %>% project("EPSG:32719")

# 3.2 Recortar (Crop) los rasters al área de Quillota antes del procesamiento matemático
# Esto ahorra drásticamente el uso de memoria RAM.
b3_ref <- crop(b3_raw, quillota_prov) 
b4_ref <- crop(b4_raw, quillota_prov)
b5_ref <- crop(b5_raw, quillota_prov)
b6_ref <- crop(b6_raw, quillota_prov)

# ------------------------------------------------------------------------------
# PASO 4: CÁLCULO DE ÍNDICES ESPECTRALES EXPLORATORIOS
# ------------------------------------------------------------------------------
# Se calculan por separado para permitir la depuración de errores.

# NDVI: Normalized Difference Vegetation Index (Vigor de la vegetación)
# Fórmula: (NIR - Red) / (NIR + Red)
ndvi <- (b5_ref - b4_ref) / (b5_ref + b4_ref)
names(ndvi) <- "NDVI"

# NDWI: Normalized Difference Water Index (McFeeters, 1996 - Cuerpos de agua)
# Fórmula: (Green - NIR) / (Green + NIR)
ndwi <- (b3_ref - b5_ref) / (b3_ref + b5_ref)
names(ndwi) <- "NDWI"

# NDMI: Normalized Difference Moisture Index (Gao, 1996 - Humedad del dosel)
# Fórmula: (NIR - SWIR1) / (NIR + SWIR1)
ndmi <- (b5_ref - b6_ref) / (b5_ref + b6_ref)
names(ndmi) <- "NDMI"

# Agrupamos los índices en un único objeto RasterStack multicapa
indices_stack <- c(ndvi, ndwi, ndmi)
plot(indices_stack, main = "Índices Espectrales - Quillota")

# ------------------------------------------------------------------------------
# PASO 5: EXTRACCIÓN DE VALORES (ZONAL STATISTICS)
# ------------------------------------------------------------------------------
# Extraemos el valor promedio (mean) de cada índice para cada polígono de cultivo.
# na.rm = TRUE asegura que píxeles vacíos o bordes no generen errores.

extraccion <- terra::extract(indices_stack, frut_quillota, fun = mean, na.rm = TRUE)

# Acoplamos los resultados extraídos
frut_quillota$NDVI_mean <- extraccion$NDVI
frut_quillota$NDWI_mean <- extraccion$NDWI
frut_quillota$NDMI_mean <- extraccion$NDMI

# ------------------------------------------------------------------------------
# PASO 6: ANÁLISIS ESTADÍSTICO Y GRÁFICO COMPARATIVO
# ------------------------------------------------------------------------------

# 6.1 Convertimos a DataFrame para trabajar cómodamente en ggplot2
df_resultados <- as.data.frame(frut_quillota) %>%
  filter(!is.na(NDVI_mean)) %>% # Eliminamos polígonos demasiado pequeños para el píxel Landsat (30m)
  select(ESPECIE, NDVI_mean, NDWI_mean, NDMI_mean)

# 6.2 Calculamos los promedios globales por Especie
resumen_especies <- df_resultados %>%
  group_by(ESPECIE) %>%
  summarise(
    Promedio_NDVI = mean(NDVI_mean, na.rm = TRUE),
    Promedio_NDWI = mean(NDWI_mean, na.rm = TRUE),
    Promedio_NDMI = mean(NDMI_mean, na.rm = TRUE),
    Conteo_Poligonos = n()
  ) %>%
  filter(Conteo_Poligonos > 10) # Filtramos especies representativas (>10 huertos)

print("Resumen de Índices por Especie Frutal en Quillota:")
print(resumen_especies)

# 6.3 Reestructuración de datos (Pivot Longer) para el gráfico de barras múltiples
df_grafico <- resumen_especies %>%
  select(ESPECIE, Promedio_NDVI, Promedio_NDWI, Promedio_NDMI) %>%
  pivot_longer(cols = starts_with("Promedio"), 
               names_to = "Indice", 
               values_to = "Valor")

# Limpieza de los nombres de los índices para la leyenda
df_grafico$Indice <- gsub("Promedio_", "", df_grafico$Indice)

# 6.4 Generación del Gráfico Comparativo de Índices
grafico_comparativo <- ggplot(df_grafico, aes(x = reorder(ESPECIE, Valor), y = Valor, fill = Indice)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.3) +
  scale_fill_manual(values = c("NDVI" = "forestgreen", "NDMI" = "dodgerblue", "NDWI" = "navy")) +
  coord_flip() + # Volteamos los ejes para facilitar la lectura de las especies
  labs(
    title = "Comparación de Índices Espectrales por Especie Frutal",
    subtitle = "Provincia de Quillota - Landsat 8 L2SP (Julio 2025)",
    x = "Especie Frutal",
    y = "Valor Promedio del Índice",
    fill = "Índice Espectral"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 10)
  )

# Mostramos el gráfico
print(grafico_comparativo)

