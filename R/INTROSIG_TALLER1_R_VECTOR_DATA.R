# ==============================================================================
# SCRIPT DE ANÁLISIS ESPACIAL CON R 
# ==============================================================================
# Autor: Prof. Lucas Vituri Santarosa
# Laboratório de Recursos Hídricos y Geotecnologias
# Institución: Pontificia Universidad Católica de Valparaíso (PUCV)
#
# El objetivo principal es cargar coberturas de uso de suelo (frutales), 
# procesar las geometrías de las provincias de la Región de Valparaíso, 
# aplicar geoprocesos de intersección y finalmente extraer la estadística 
# descriptiva de las áreas cultivadas.
# ==============================================================================

# ------------------------------------------------------------------------------
# PASO 1: CONFIGURACIÓN DEL ENTORNO Y LIBRERÍAS
# ------------------------------------------------------------------------------

# Instalación de las librerías necesarias. 
# (Descomentar estas líneas si se ejecuta el script por primera vez)
# install.packages("terra")
# install.packages("tidyterra")
# install.packages("dplyr")
# install.packages("ggplot2")

# Cargamos las bibliotecas.
# 'terra' es la librería principal y más moderna para procesar datos espaciales.
library(terra) 
library(tidyterra)
library(dplyr)    
library(ggplot2)  

# ------------------------------------------------------------------------------
# PASO 2: CARGA Y EXPLORACIÓN DE DATOS ESPACIALES
# ------------------------------------------------------------------------------
# Siempre analicen críticamente sus inputs antes de operar con ellos.

# Cargamos el shapefile de la cubierta frutal.
# Nota: Ajusta la ruta según dónde guardaste tus archivos.
frut <- vect("~/03_AULAS/01_POS/01_CLASES/2026/01_INTROSIG/02_TALLER/02_TALLER1_VECTOR/00_DATOS/00_VECTOR/01_FRUTALES/5__cubierta_frutal_de_pol_gonos_regi_n_de_valpara_so_2025/5__cubierta_frutal_de_pol_gonos_regi_n_de_valpara_so_2025.shp")

# Verificamos la proyección, los nombres de atributos y la extensión
print(frut)

# Cargamos el shapefile que contiene el límite de las Provincias
prov <- vect('~/03_AULAS/01_POS/01_CLASES/2026/01_INTROSIG/02_TALLER/02_TALLER1_VECTOR/00_DATOS/00_VECTOR/04_REGIONES_PROVINCIAS/01_PROVINCIAS/Provincias.shp')
print(prov)

# ------------------------------------------------------------------------------
# PASO 3: PROCESAMIENTO Y FILTRADO ESPACIAL INTERMEDIO
# ------------------------------------------------------------------------------
# Es recomendable ir paso a paso, no intenten simplificar todo en una línea si 
# la lógica aún no está clara. Desmembrar el procesamiento ayuda a la depuración.

# Creamos un objeto específico para la Región de Valparaíso
valp_prov <- prov %>%
  # disagg() rompe polígonos multiparte en partes individuales para evitar errores
  terra::disagg() %>%
  # Filtramos usando dplyr para quedarnos solo con la región de interés
  filter(Region == "Región de Valparaíso") %>%
  # Reproyectamos al sistema WGS 84/UTM zona 19S 
  project("EPSG:32719")

# Calculamos el área de cada provincia en hectáreas (ha) 
valp_prov$area_ha <- terra::expanse(valp_prov, unit = "ha")

# Filtramos polígonos continentales mayores a 100.000 ha.
# Análisis crítico: Esto nos ayuda a eliminar polígonos que corresponden a islas
# pequeñas, rocas o artefactos topológicos que ensuciarán el resultado.
valp_prov_cont <- valp_prov %>%
  filter(area_ha > 100000)

# Verificamos el tamaño de las áreas que sobrevivieron al filtro
print(valp_prov_cont$area_ha)

# Vemos la estructura del nuevo polígono limplio
print(valp_prov_cont)

# Visualizamos las provincias continentales de Valparaíso para validar
plot(valp_prov_cont, main = "Provincias Continentales - Región de Valparaíso")

# ------------------------------------------------------------------------------
# PASO 4: INTERSECCIÓN ESPACIAL
# ------------------------------------------------------------------------------
# Extraeremos la información de la cubierta frutal que se encuentre ESTRICTAMENTE
# dentro de las provincias que filtramos.

inter_fr_pr <- terra::intersect(valp_prov_cont, frut)

# Validamos que los nombres de los atributos de ambos shapefiles se hayan combinado
print(names(inter_fr_pr))

# Visualizamos la intersección. Esto nos dará una idea de la densidad agrícola.
plot(inter_fr_pr, main = "Intersección: Cubierta Frutal por Provincias")

# ------------------------------------------------------------------------------
# PASO 5: LIMPIEZA DE ATRIBUTOS Y RE-CÁLCULO GEOMÉTRICO
# ------------------------------------------------------------------------------
# Después de un corte (intersect), los polígonos frutales en los bordes de la
# provincia cambiaron de tamaño. Las hectáreas anteriores ya no son reales. 

# Dejamos solo los datos que necesitamos en la tabla para no saturar la memoria
frut_prov <- select(inter_fr_pr, 'Provincia', 'ESPECIE', 'ANO_PLANT')

# Estandarizamos el nombre de la columna para evitar confusiones de sintaxis
frut_prov$PROVINCIA <- frut_prov$Provincia
frut_prov <- select(frut_prov, -'Provincia')

# Recalculamos la nueva área post-intersección de cada polígono frutal
frut_prov$AREA_HA <- terra::expanse(frut_prov, unit = "ha")

# Exploramos el objeto final
print(frut_prov)

# ------------------------------------------------------------------------------
# PASO 6: ANÁLISIS ESTADÍSTICO Y GENERACIÓN DE RESULTADOS
# ------------------------------------------------------------------------------

# Para los resúmenes estadísticos, convertimos los datos espaciales a un 
# data.frame normal. Esto hace el procesamiento mucho más rápido y comprensible.

# 1. Tabla resumen del área total por provincia
tabela_provincia <- frut_prov %>%
  as.data.frame() %>%
  group_by(PROVINCIA) %>%
  summarise(TOTAL_AREA_HA = sum(AREA_HA, na.rm = TRUE)) %>%
  arrange(desc(TOTAL_AREA_HA)) # Orden jerárquico descendente

print("Área Total de Frutales por Provincia (ha):")
print(tabela_provincia)

# 2. Tabla resumen del área total por especie frutal
tabela_especie <- frut_prov %>%
  as.data.frame() %>%
  group_by(ESPECIE) %>%
  summarise(TOTAL_AREA_HA = sum(AREA_HA, na.rm = TRUE)) %>%
  arrange(desc(TOTAL_AREA_HA))

print("Top 10: Área Total por Especie Frutal (ha):")
print(head(tabela_especie, 10))

# 3. Gráfico de Barras 
# Utilizamos las funciones base de R para levantar un gráfico exploratorio
barplot(tabela_provincia$TOTAL_AREA_HA,
        names.arg = tabela_provincia$PROVINCIA,
        las = 2,                 # Permite que los nombres en el eje X giren 90 grados
        col = "steelblue",       # Color a elección
        main = "Área Total de Frutas por Provincia",
        ylab = "Área (ha)")

# ------------------------------------------------------------------------------
# 7. MAPA TEMÁTICO: FRUTA PREDOMINANTE POR PROVINCIA
# ------------------------------------------------------------------------------

# Paso 7.1: Identificar la fruta predominante en cada provincia
# Agrupamos por Provincia y Especie para sumar las áreas
resumen_especies <- frut_prov %>%
  as.data.frame %>% 
  group_by(PROVINCIA, ESPECIE) %>%
  summarise(AREA_TOTAL = sum(AREA_HA, na.rm = TRUE), .groups = 'drop')

# Para cada provincia, filtramos y nos quedamos solo con la especie de mayor área
fruta_predominante <- resumen_especies %>%
  group_by(PROVINCIA) %>%
  slice_max(order_by = AREA_TOTAL, n = 1) %>%
  rename(FRUTA_PRINCIPAL = ESPECIE) # Renombramos para mayor claridad

print(fruta_predominante)

# Paso 7.2: Unir los datos tabulares con los polígonos espaciales
# Hacemos un 'merge' (unión) entre el vector espacial de las provincias y nuestra tabla.
# 'by.x' es el nombre de la columna en el mapa y 'by.y' en la tabla.
valp_prov_mapa <- merge(valp_prov_cont, 
                        fruta_predominante, 
                        by.x = "Provincia", 
                        by.y = "PROVINCIA", 
                        all.x = TRUE)

# Paso 7.3: Dibujar el mapa utilizando ggplot2 y tidyterra
# geom_spatvector es una función de tidyterra que permite usar shapefiles en ggplot
mapa_final <- ggplot() +
  geom_spatvector(data = valp_prov_mapa, 
                  aes(fill = FRUTA_PRINCIPAL), # El color de relleno dependerá de la fruta
                  color = "black",             # Color de los bordes de la provincia
                  linewidth = 0.5) +           # Grosor del borde
  scale_fill_brewer(palette = "Set2", 
                    name = "Especie Predominante", 
                    na.value = "white") +      # Colores categóricos para las frutas
  labs(
    title = "Especie Frutal Predominante por Provincia",
    caption = "Fuente: ODEPA (2025)"
  ) +
  theme_minimal() +                            # Tema limpio y profesional
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right"
  )

# Mostramos el mapa generado
print(mapa_final)

