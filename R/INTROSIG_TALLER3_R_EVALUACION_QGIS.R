# ==============================================================================
# EVALUACIÓN PRÁCTICA: ESTRUCTURA FRUTÍCOLA Y DINÁMICA FENOLÓGICA
# ==============================================================================
# Institución: Pontificia Universidad Católica de Valparaíso (PUCV)
# Prof. Lucas Vituri Santarosa
#
# Descripción: Diagnóstico integral combinando vectores y rasters.
# Fase I: Densidad de riego y superficie frutícola provincial.
# Fase II & III: Análisis de NDVI estacional y estadísticas zonales.
# ==============================================================================

# Cargar librerías necesarias
library(terra)      # Motor espacial principal
library(dplyr)      # Manejo de tablas de atributos
library(ggplot2)    # Creación de mapas de calidad cartográfica
library(tidyterra)  # Integración de terra con ggplot2

# ------------------------------------------------------------------------------
# FASE I: ANÁLISIS VECTORIAL Y SELECCIÓN DEL ÁREA DE ESTUDIO
# ------------------------------------------------------------------------------

# 1. Cargar y estandarizar datos vectoriales (Proyección EPSG:32719)

setwd('03_AULAS/01_POS/01_CLASES/2026/01_INTROSIG/02_TALLER/05_TALLER3_VECTOR_RASTER_R/')

prov <- vect("DATOS/04_REGIONES_PROVINCIAS/01_PROVINCIAS/Provincias.shp") %>% project("EPSG:32719")
canales <- vect("DATOS/02_CANALES_RIEGO/canales_nacional_final.shp") %>% project("EPSG:32719")
frutales <- vect("DATOS/01_FRUTALES/5__cubierta_frutal_de_pol_gonos_regi_n_de_valpara_so_2025/5__cubierta_frutal_de_pol_gonos_regi_n_de_valpara_so_2025.shp") %>% project("EPSG:32719")

# 2. Área Provincial (Base para los cálculos de densidad y porcentaje)
# Separamos geometrías por seguridad y calculamos área en km² y hectáreas.

prov <- terra::disagg(prov) %>% 
  filter(Region == "Región de Valparaíso") %>% 
  mutate(Area_km2 = terra::expanse(.[], unit = "km")) %>% 
  mutate(Area_ha = terra::expanse(.[], unit = "ha")) %>% 
  filter(Area_ha > 100000)
plot(prov)

# --- TAREA A: Densidad de Riego (km/km²) ---
# Intersección espacial explícita para evitar pérdida de geometría
canales_prov <- terra::intersect(canales, prov)
canales_prov$Longitud_km <- terra::perim(canales_prov) / 1000 # perim() calcula longitud en líneas

# Agregación tabular de densidad
tabla_densidad_riego <- as.data.frame(canales_prov) %>%
  group_by(Provincia) %>%
  summarise(
    Total_Canales_km = sum(Longitud_km, na.rm = TRUE),
    Area_Provincial_km2 = first(Area_km2)
  ) %>%
  mutate(Densidad_km_km2 = Total_Canales_km / Area_Provincial_km2) %>%
  arrange(desc(Densidad_km_km2))

print("TABLA 1: Densidad de Riego por Provincia")
print(tabla_densidad_riego)

# --- TAREA B: Superficie Frutícola (%) ---
frutales_prov <- terra::intersect(frutales, prov)
frutales_prov$Area_Frutal_ha <- terra::expanse(frutales_prov, unit = "ha")

tabla_tabla_tabla_superficie_fruticola <- as.data.frame(frutales_prov) %>%
  group_by(Provincia) %>%
  summarise(
    Total_Frutales_ha = sum(Area_Frutal_ha, na.rm = TRUE),
    Area_Provincial_ha = first(Area_ha)
  ) %>%
  mutate(Porcentaje_Fruticola = (Total_Frutales_ha / Area_Provincial_ha) * 100) %>%
  arrange(desc(Total_Frutales_ha))

print("TABLA 2: Superficie Frutícola por Provincia")
print(tabla_superficie_fruticola)

# --- TAREA C: Selección del Área Prioritaria y Mapa ---
# Identificamos la provincia con mayor superficie plantada (Top 1 de la tabla anterior)
nombre_prov_prioridad <- tabla_superficie_fruticola$Provincia[1]

# Subconjuntos espaciales usando indexación nativa (¡A prueba de errores!)
prov_prioridad <- prov[prov$Provincia == nombre_prov_prioridad, ]
frutales_prioridad <- frutales_prov[frutales_prov$Provincia == nombre_prov_prioridad, ]

# MAPA 1: Cobertura de Frutales en el área prioritaria
mapa_frutales <- ggplot() +
  geom_spatvector(data = prov_prioridad, fill = "gray95", color = "black") +
  geom_spatvector(data = frutales_prioridad, aes(fill = ESPECIE), color = NA) +
  labs(
    title = paste("Distribución Frutícola - Provincia de", nombre_prov_prioridad),
    subtitle = "Fase I: Área Prioritaria",
    fill = "Especie"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(mapa_frutales)

# ------------------------------------------------------------------------------
# FASE II & III: ANÁLISIS RÁSTER, ESTADÍSTICAS ZONALES Y DELTA NDVI
# ------------------------------------------------------------------------------

# 1. Cargar imágenes NDVI descargadas de GEE
# (Asumiendo que las imágenes ya vienen con el factor de escala aplicado de GEE)
ndvi_verano <- rast("DATOS/Landsat-8_Quillota_ÍndiceNDVI_Verano/Landsat-8_Quillota_ÍndiceNDVI_VERANO.tif") %>% project("EPSG:32719")
ndvi_invierno <- rast("DATOS/Landsat-8_Quillota_ÍndiceNDVI_Invierno/Landsat-8_Quillota_ÍndiceNDVI_INVIERNO.tif") %>% project("EPSG:32719")

# Renprov_prioridad# Renombrar capas
names(ndvi_verano) <- "NDVI_Verano"
names(ndvi_invierno) <- "NDVI_Invierno"
ndvi_stack <- c(ndvi_verano, ndvi_invierno)

# --- TAREA D: Estadísticas Zonales ---
# Extraemos media, desviación estándar y máximo. 
# Lo hacemos por separado para evitar funciones complejas que saturen la RAM.
ext_mean <- terra::extract(ndvi_stack, frutales_prioridad, fun = mean, na.rm = TRUE)
ext_sd   <- terra::extract(ndvi_stack, frutales_prioridad, fun = sd, na.rm = TRUE)
ext_max  <- terra::extract(ndvi_stack, frutales_prioridad, fun = max, na.rm = TRUE)

# Incorporamos los datos al vector original
frutales_prioridad$NDVI_Verano_Mean <- ext_mean$NDVI_Verano
frutales_prioridad$NDVI_Verano_SD   <- ext_sd$NDVI_Verano
frutales_prioridad$NDVI_Verano_Max  <- ext_max$NDVI_Verano

frutales_prioridad$NDVI_Invier_Mean <- ext_mean$NDVI_Invierno
frutales_prioridad$NDVI_Invier_SD   <- ext_sd$NDVI_Invierno
frutales_prioridad$NDVI_Invier_Max  <- ext_max$NDVI_Invierno

# Tabla final de diagnóstico
tabla_estadisticas_zonales <- as.data.frame(frutales_prioridad) %>%
  select(ESPECIE, Area_Frutal_ha, starts_with("NDVI")) %>%
  filter(!is.na(NDVI_Verano_Mean)) # Limpiar polígonos sin datos válidos

print("TABLA 3: Estadísticas Zonales por Huerto (Muestra)")
print(head(tabla_estadisticas_zonales))

# --- TAREA E: Análisis de Cambio (Delta NDVI) ---
# Cálculo del álgebra de mapas
delta_ndvi <- ndvi_verano - ndvi_invierno
names(delta_ndvi) <- "Delta_NDVI"

quillota <- prov %>% filter(Provincia == "Quillota")

# MAPA 2: Mapa Delta NDVI
# Valores positivos (verde) indican mayor vigor en verano; negativos (marrón) mayor en invierno.
mapa_delta <- ggplot() +
  geom_spatraster(data = delta_ndvi) +
  geom_spatvector(data = quillota, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_whitebox_c(
    palette = "muted", 
    n.breaks = 7, 
    limits = c(-0.5, 0.5), # Ajustar límites según la realidad de la región
    labels = scales::label_number(accuracy = 0.1)
  ) +
  labs(
    title = paste("Dinámica Fenológica (Delta NDVI) -", nombre_prov_prioridad),
    subtitle = "Variación Estacional: Verano - Invierno",
    fill = "Δ NDVI"
  ) +
  theme_minimal()

print(mapa_delta)

# ==============================================================================
# FIN DEL DIAGNÓSTICO
# ==============================================================================
