## Identifying structural and individual-level correlates of sociodemographic inequalities on CVD mortality in Mexico				
## Data Analysis: Mario Cesar Torres Chaves and Neftali Eduardo Antonio-Villa
## Latest version of Analysis: January-2026
## Any question regarding analysis, please contact Neftali Eduardo Antonio-Villa (neftalivilla@comunidad.unam.mx)


#######Codigo necesario para analisis posteriores ######

#---- Cargar paquetes ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rlang,foreign, knitr, dplyr, ggplot2, sf, readr, janitor, stringr, spdep, tmap, readxl,tidyr,purrr,rlang,classInt, glmmTMB, broom,ggsci )

###### Parametros modificables ######

# Se guarda el Directorio en dónde se almancenaron las bases de datos
#Se estandariza la ruta del directorio para su uso sin contratiempos en R

#ruta_base <- normalizePath("C:/Users/enriq/OneDrive/MORTALIDAD INEGI/Bases de datos", winslash = "/", mustWork = TRUE)
ruta_base <- normalizePath("/Users/ctorres92/Library/CloudStorage/OneDrive-Personal/Lown Scholar/MORTALIDAD INEGI/Bases de datos", winslash = "/", mustWork = TRUE)
ruta_figuras <- normalizePath("/Users/ctorres92/Library/CloudStorage/OneDrive-Personal/Lown Scholar/Cardiovascular mortality/Figures", winslash = "/", mustWork = TRUE)

#Años que se desean evaluar
anios <- 1990:2024

#Codigos CIE de causa de defunción de interes
causa_def_interes <- list(
  icd9  = as.character(390:459),                      # 390–459
  icd10 = paste0("I", sprintf("%02d", 10:99))          # I00–I99
)

#Comprobar existencia de capetas y archivos necesarios
comprobar_carpetas <- function(ruta_base) {
  
  #En caso de que no exista un directorio en ruta_base se detiene el codigo y arroja un mensaje
  
  if (!dir.exists(ruta_base)) {
    warning ("No se encontraró la ruta ingresada'\n") 
    return(NULL)}
  
  #Se almacenan los archivos que coincidan con "defunciones_base_datos" que es el nombre que se utiliza en las bases de datos de registros de defunciones IENGI
  carpetas <- list.files(ruta_base, pattern = "defunciones_base_datos", full.names = TRUE)
  
  #En caso de que no existan carpetas en la ruta_base se arroja un mensaje
  if (length(carpetas) == 0 ) {
    warning("No se encontraron carpetas de registro de defunciones en la ruta ingresada\n")
    return(NULL)}
  
  #Se elige al azar una de las carpetas que contienen los registros de defuncion para mostrar los archivos que contiene
  random_files <- list.files(sample(list.files(ruta_base, pattern = "defunciones_base_datos", full.names = TRUE), 1))
  
  #Se genera texto en la consola para obsservar 1) ruta_base 2) Carpetas almacenadas en la ruta base y 3) los archivos almancenados en una carpeta aleatoria
  cat ("Ruta ingresada:\n")
  print(ruta_base)
  cat ("\n\n")
  
  cat ("Carpetas de registros de defunciones encotradas:\n")
  print(basename(carpetas))
  cat ("\n\n")
  
  cat("Los archivos almacenados dentro de una carpeta aleatoria son:\n")
  print(random_files)
  cat ("\n\n")
  
}

comprobar_carpetas(ruta_base = ruta_base)


###### Generar las rutas para los registros de defunciones ######

generacion_rutas_defun <-  function(anios, ruta_base) {
  
  #Se almacenan las rutas completas de los registros de defunciones
  carpetas <- list.files(ruta_base, pattern = "defunciones_base_datos", full.names = TRUE)
  
  #Funcion que genera a partir de los años de las carpetas identificadas un patron DEFUN(dos útlimos digitos del año) 
  #DEFUNaño es el nombre de los archivos que contienen las bases de datos de los regstros de defunción
  ruta_anio <- lapply(anios, function(anio) {
    
    yy <- sprintf("%02d", anio %% 100)
    patron <- sprintf("DEFUN%s", yy)
    
    #Funcion interna para generar un vector que contiene las rutas usando el patron DEFUNaño
    rutas_defun <- unlist(lapply(carpetas, function (defun) {
      list.files(defun, pattern = patron, full.names = TRUE,ignore.case = TRUE)}))
    
    #Mensaje de Warning en caso de que no existan rutas con el patron DEFUNaño para cada año en especifico  
    if (length(rutas_defun) == 0) {
      warning(sprintf("No se encontraron archivos %s para el año %s",patron, anio))
    }
    return(rutas_defun)
  }
  )
  return(ruta_anio)
}

#Se genera una lista usando la función anterior dónde se encuentran las rutas los registros de defuncion de cada archivo por año
rutas_defun <- generacion_rutas_defun (anios = anios, ruta_base = ruta_base)


#---- Cargar los dataframe que contienen los registros de defuncion INEGI en una lista ----
all_deaths <- setNames(lapply(rutas_defun, foreign::read.dbf, as.is = FALSE), anios) 

#---- Función para clasificar muertes ----
clasificar_cvd <- function(anios, all_deaths, causa_def_interes) {
  
  resultado <- lapply(anios, function(anio) {
    
    registros <- all_deaths[[as.character(anio)]]
    if (is.null(registros)) return(NULL)
    
    # Limpieza de CAUSA_DEF
    registros$CAUSA_DEF <- trimws(gsub("\\.", "", as.character(registros$CAUSA_DEF)))
    
    # Patrón por sistema de codificación según año
    if (anio <= 1997) {
      patron <- paste0("^(", paste(causa_def_interes$icd9, collapse = "|"), ")")
    } else {
      patron <- paste0("^(", paste(causa_def_interes$icd10, collapse = "|"), ")")
    }
    
    # Variable indicadora
    registros$muerte_cvd <- ifelse(grepl(patron, registros$CAUSA_DEF), 1L, 0L)
    
    registros
  })
  
  names(resultado) <- as.character(anios)
  resultado
}

#---- Aplicar clasificación ----
classified_deaths <- clasificar_cvd(
  anios = anios,
  all_deaths = all_deaths,
  causa_def_interes = causa_def_interes
)

#---- Cargar datos población CONAPO ----

ruta_poblacion_conapo <- function(ruta_base) {
  
  #Almacenar la ruta de los archivos de la carpeta 00_Republica_mexicana que contienene las proyecciones de poblacion de CONAPO
  carpeta <- list.files(ruta_base, pattern = "00_Republica_mexicana", full.names = TRUE, ignore.case = TRUE)
  
  
  if(!file.exists(carpeta)) {
    warning("No se encontró ninguna carpeta con el nombre '00_Republica_mexicana' ")
    return(NULL)}
  
  #Almacenar la ruta de los archivos de la carpeta 3_indicadores_Dem..." que contiene el DF de las poblaciones
  archivo <- list.files(carpeta, pattern = "3_Indicadores_Dem_00_RM" , full.names = TRUE, ignore.case = TRUE)
  
  if (length(archivo) == 0) {
    warning("No se encontró ningun archivo dentro de '00_Republica_mexicana'")
    return(NULL)}
  else { 
    archivos <- list.files(carpeta, full.names = TRUE)
    cat("La carpeta '00_Republica_mexicana' contiene los siguientes archivos:\n\n", paste(archivos, collapse ="\n"),"\n")}
  
  return(archivo[1])
  
}

#Almacenar la ruta
archivo_pob <- ruta_poblacion_conapo(ruta_base = ruta_base)

#Cargar el archivo
poblacion_conapo <- read_excel(archivo_pob, sheet = 1)


#---- Cargar datos población CONAPO quinquenal ----
ruta_poblacion_quinquenal <- function(ruta_base) {
  
  #Almacenar la ruta de los archivos de la carpeta 00_Republica_mexicana que contienene las proyecciones de poblacion de CONAPO
  carpeta <- list.files(ruta_base, pattern = "00_Republica_mexicana", full.names = TRUE, ignore.case = TRUE)
  
  
  if(!file.exists(carpeta)) {
    warning("No se encontró ninguna carpeta con el nombre '00_Republica_mexicana' ")
    return(NULL)}
  
  #Almacenar la ruta de los archivos de la carpeta que contiene la población segmentada por grupos quinquenales que contiene el DF de las poblaciones
  archivo <- list.files(carpeta, pattern = "^1_Grupo_Quinq_00_RM\\.xlsx$" , full.names = TRUE, ignore.case = TRUE)
  
  if (length(archivo) == 0) {
    warning("No se encontró ningun archivo dentro de '1_Grupo_Quinq_00_RM'")
    return(NULL)}
  else {
    archivos <- list.files(carpeta, full.names = TRUE)
    cat("La carpeta '00_Republica_mexicana' contiene los siguientes archivos:\n\n", paste(archivos, collapse ="\n"),"\n")}
  
  return(archivo)
  
}
#Almanenar ruta

archivo_pop_quinquenal <- ruta_poblacion_quinquenal(ruta_base = ruta_base)

poblacion_conapo_quinquenal <- read_excel(archivo_pop_quinquenal, sheet = 1)

poblacion_conapo_quinquenal <- poblacion_conapo_quinquenal %>%
  dplyr::mutate(
    cvegeo = sprintf("%05d", suppressWarnings(as.integer(as.character(CLAVE))))
  )


#---- Generar CLAVE Municipal en registros de defuncion filtrados----
generacion_clave_mun <- function(anios, registros_lista) {
  
  
  clave_anio <- lapply(anios, function(anio){
    
    #Almacenar los registros de defuncion con el año como character
    regis <- registros_lista[[as.character(anio)]]
    
    #Generar una columna con el nombre cvegeo usando los codigos Estatales (2 dígitos) y Municipales (3 digitos)
    regis <- regis %>% mutate(ENT_RESID = sprintf("%02d", as.integer(ENT_RESID)), MUN_RESID = sprintf("%03d", as.integer(MUN_RESID)), cvegeo = paste0(ENT_RESID, MUN_RESID))
    
    return(regis)
  })
  
  names(clave_anio) <- as.character(anios)
  
  return(clave_anio)
}

classified_deaths_id <- generacion_clave_mun (anios = anios, registros_lista = classified_deaths)


####################
#Estatdísticos descriptivos
####################

#---- Tasas de Mortalidad  ----
mortalidad_nacional <- function(anios, registros_defun_clave_mun, poblacion_conapo) {
  
  out <- lapply(anios, function(anio){
    
    df <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(df)) return(NULL)
    
    # Defunciones
    def_total <- nrow(df)
    def_cvd   <- sum(df$muerte_cvd, na.rm = TRUE)
    
    # Población nacional (suma municipal CONAPO)
    pop_total <- poblacion_conapo %>%
      dplyr::filter(AÑO == anio) %>%
      dplyr::summarise(poblacion_total = sum(POB_MIT_MUN, na.rm = TRUE)) %>%
      dplyr::pull(poblacion_total)
    
    # Tasas por 100,000
    tasa_total <- def_total / pop_total * 100000
    tasa_cvd   <- def_cvd   / pop_total * 100000
    
    dplyr::tibble(
      anio = anio,
      def_total = def_total,
      def_cvd = def_cvd,
      poblacion = pop_total,
      tasa_total = tasa_total,
      tasa_cvd = tasa_cvd,
      prop_cvd = def_cvd / def_total
    )
  })
  
  dplyr::bind_rows(out)
}

mortality_national_crude <- mortalidad_nacional(
  anios = anios,
  registros_defun_clave_mun = classified_deaths_id,   
  poblacion_conapo = poblacion_conapo
)


# CVD deaths
ggplot(mortality_national_crude, 
       aes(x = anio, y = tasa_cvd)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Cardiovascular Mortality Rate in Mexico",
    subtitle = "Crude rate per 100,000 population",
    x = "Year",
    y = "Cardiovascular mortality rate (per 100,000 population)"
  ) +
  theme_minimal()



#State mortality rates
mortality_state_crude <- function(anios, registros_defun_clave_mun, poblacion_conapo) {
  
  mortalidad_anio <- lapply(anios, function(anio) {
    
    df <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(df)) return(NULL)
    
    # Defunciones CVD por estado (residencia)
    defun_estatal <- df %>%
      mutate(ENT_RESID = sprintf("%02d", suppressWarnings(as.integer(as.character(ENT_RESID))))) %>%
      group_by(ENT_RESID) %>%
      summarise(defunciones_cvd = sum(muerte_cvd, na.rm = TRUE), .groups = "drop")
    
    # Población estatal CONAPO
    poblacion_estatal <- poblacion_conapo %>%
      filter(AÑO == anio) %>%
      mutate(CLAVE_ENT = sprintf("%02d", suppressWarnings(as.integer(as.character(CLAVE_ENT))))) %>%
      group_by(CLAVE_ENT) %>%
      summarise(
        poblacion_total = sum(POB_MIT_MUN, na.rm = TRUE),
        estado = first(NOM_ENT),
        .groups = "drop"
      ) %>%
      rename(ENT_RESID = CLAVE_ENT)
    
    # Tasa CVD por 100,000
    df_estatal <- defun_estatal %>%
      inner_join(poblacion_estatal, by = "ENT_RESID") %>%
      mutate(
        tasa_cvd = defunciones_cvd / poblacion_total * 100000,
        anio = anio
      )
    
    df_estatal
  })
  
  dplyr::bind_rows(mortalidad_anio)
}

mortality_state_crude <- mortality_state_crude(
  anios = anios,
  registros_defun_clave_mun = classified_deaths_id,
  poblacion_conapo = poblacion_conapo
)


#Generar gráfico tasa mortalidad estatal 
ggplot(mortality_state_crude, 
       aes(x = anio, y = tasa_cvd, group = estado, color = estado)) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Cardiovascular Mortality Rate Trends by State in Mexico",
    subtitle = "Crude rate per 100,000 population",
    x = "Year",
    y = "Cardiovascular mortality rate (per 100,000 population)",
    color = "State"
  ) +
  theme_minimal()

#Tasas de mortalidad por macroregion
mortalidad_macroregion_crude <- function(anios, registros_defun_clave_mun, poblacion_conapo) {
  
  asignar_macro <- function(ent) {
    dplyr::case_when(
      ent %in% c(2, 3, 5, 8, 19, 25, 26, 28) ~ 1,  # North
      ent %in% c(1, 6, 10, 11, 14, 16, 18, 24, 32) ~ 2,  # Central-West
      ent %in% c(9, 13, 15, 17, 21, 22, 29) ~ 3,  # Central
      ent %in% c(4, 7, 12, 20, 23, 27, 30, 31) ~ 4,  # South-Southeast
      TRUE ~ NA_real_
    )
  }
  
  mortalidad_anio <- lapply(anios, function(anio) {
    
    df <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(df)) return(NULL)
    
    # CVD deaths by macro-region (ENT_OCURR)
    defun_macro <- df %>%
      mutate(
        ent = suppressWarnings(as.integer(as.character(ENT_OCURR))),
        macro_cod = asignar_macro(ent)
      ) %>%
      filter(!is.na(macro_cod)) %>%
      group_by(macro_cod) %>%
      summarise(defunciones_cvd = sum(muerte_cvd, na.rm = TRUE), .groups = "drop")
    
    # Population by macro-region (CLAVE_ENT)
    pob_macro <- poblacion_conapo %>%
      filter(AÑO == anio) %>%
      mutate(
        ent = suppressWarnings(as.integer(as.character(CLAVE_ENT))),
        macro_cod = asignar_macro(ent)
      ) %>%
      filter(!is.na(macro_cod)) %>%
      group_by(macro_cod) %>%
      summarise(poblacion = sum(POB_MIT_MUN, na.rm = TRUE), .groups = "drop")
    
    pob_macro %>%
      left_join(defun_macro, by = "macro_cod") %>%
      mutate(
        defunciones_cvd = ifelse(is.na(defunciones_cvd), 0L, defunciones_cvd),
        tasa_cvd = defunciones_cvd / poblacion * 100000,
        anio = anio
      ) %>%
      arrange(macro_cod)
  })
  
  dplyr::bind_rows(mortalidad_anio)
}

mortalidad_macroregion_crude <- mortalidad_macroregion_crude(
  anios = anios,
  registros_defun_clave_mun = classified_deaths_id,  # tu lista con muerte_cvd
  poblacion_conapo = poblacion_conapo
)

#Generar gráfico tasa mortalidad macroregion 
mortalidad_macroregion_crude %>%
  mutate(
    macroregion = factor(
      macro_cod,
      levels = c(1, 2, 3, 4),
      labels = c("North", "Central-West", "Central", "South-Southeast")
    )
  ) %>%
  ggplot(aes(
    x = anio,
    y = tasa_cvd,
    group = macroregion,
    color = macroregion
  )) +
  geom_line(alpha = 0.8, linewidth = 0.9) +
  labs(
    title = "Cardiovascular Mortality Rate Trends by Macro-Region",
    subtitle = "Mexico",
    x = "Year",
    y = "Cardiovascular mortality rate (per 100,000 population)",
    color = "Macro-region"
  ) +
  theme_minimal()


#Tasas de mortalidad Municipal
mortality_municipal_crude <- function(anios, registros_defun_clave_mun, poblacion_conapo) {
  
  mortalidad_anio <- lapply(anios, function(anio) {
    
    # Población municipal del año
    poblacion_mun <- poblacion_conapo %>%
      dplyr::filter(AÑO == anio) %>%
      dplyr::mutate(cvegeo = sprintf("%05d", suppressWarnings(as.integer(as.character(CLAVE))))) %>%
      dplyr::group_by(cvegeo) %>%
      dplyr::summarise(
        poblacion_total = sum(POB_MIT_MUN, na.rm = TRUE),
        municipio = dplyr::first(NOM_MUN),
        estado = dplyr::first(NOM_ENT),
        .groups = "drop"
      )
    
    # Registros de defunción del año
    regis <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(regis)) return(NULL)
    
    # Defunciones CVD por municipio
    defun_mun <- regis %>%
      dplyr::mutate(cvegeo = sprintf("%05d", suppressWarnings(as.integer(as.character(cvegeo))))) %>%
      dplyr::group_by(cvegeo) %>%
      dplyr::summarise(defunciones_cvd = sum(muerte_cvd, na.rm = TRUE), .groups = "drop")
    
    # Unir y calcular tasa CVD por 100,000
    df_municipal <- poblacion_mun %>%
      dplyr::left_join(defun_mun, by = "cvegeo") %>%
      dplyr::mutate(
        defunciones_cvd = dplyr::coalesce(defunciones_cvd, 0L),
        tasa_cvd = defunciones_cvd / poblacion_total * 100000,
        anio = anio
      )
    
    df_municipal
  })
  
  dplyr::bind_rows(mortalidad_anio)
}

mortality_municipal_crude <- mortality_municipal_crude(
  anios = anios,
  registros_defun_clave_mun = classified_deaths_id,  # tu lista con muerte_cvd
  poblacion_conapo = poblacion_conapo
)


# ---- Tasas de mortalidad ajustada por edad ----
clasificacion_edad_quinq <- function(registros_defun_clave_mun, anios) {
  
  out <- lapply(anios, function(anio){
    
    df <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(df)) return(NULL)
    
    # Ensure EDAD_AGRU is 2-char string ("01"..."30")
    df <- df %>%
      mutate(
        EDAD_AGRU = sprintf("%02d", suppressWarnings(as.integer(as.character(EDAD_AGRU)))),
        grupo_edad = case_when(
          EDAD_AGRU %in% c("01","02","03","04","05") ~ "POB_00_04",
          EDAD_AGRU == "06" ~ "POB_05_09",
          EDAD_AGRU == "07" ~ "POB_10_14",
          EDAD_AGRU == "08" ~ "POB_15_19",
          EDAD_AGRU == "09" ~ "POB_20_24",
          EDAD_AGRU == "10" ~ "POB_25_29",
          EDAD_AGRU == "11" ~ "POB_30_34",
          EDAD_AGRU == "12" ~ "POB_35_39",
          EDAD_AGRU == "13" ~ "POB_40_44",
          EDAD_AGRU == "14" ~ "POB_45_49",
          EDAD_AGRU == "15" ~ "POB_50_54",
          EDAD_AGRU == "16" ~ "POB_55_59",
          EDAD_AGRU == "17" ~ "POB_60_64",
          EDAD_AGRU == "18" ~ "POB_65_69",
          EDAD_AGRU == "19" ~ "POB_70_74",
          EDAD_AGRU == "20" ~ "POB_75_79",
          EDAD_AGRU == "21" ~ "POB_80_84",
          EDAD_AGRU %in% c("22","23","24","25","26","27","28","29") ~ "POB_85_mm",
          TRUE ~ NA_character_  # includes "30" no especificada
        )
      )
    
    df
  })
  
  names(out) <- as.character(anios)
  out
}

# ---- Main function: age-adjusted rates ----
clasificacion_edad_quinq <- function(registros_defun_clave_mun, anios) {
  
  out <- lapply(anios, function(anio){
    
    df <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(df)) return(NULL)
    
    # Ensure EDAD_AGRU is 2-char string ("01"..."30")
    df <- df %>%
      mutate(
        EDAD_AGRU = sprintf("%02d", suppressWarnings(as.integer(as.character(EDAD_AGRU)))),
        grupo_edad = case_when(
          EDAD_AGRU %in% c("01","02","03","04","05") ~ "POB_00_04",
          EDAD_AGRU == "06" ~ "POB_05_09",
          EDAD_AGRU == "07" ~ "POB_10_14",
          EDAD_AGRU == "08" ~ "POB_15_19",
          EDAD_AGRU == "09" ~ "POB_20_24",
          EDAD_AGRU == "10" ~ "POB_25_29",
          EDAD_AGRU == "11" ~ "POB_30_34",
          EDAD_AGRU == "12" ~ "POB_35_39",
          EDAD_AGRU == "13" ~ "POB_40_44",
          EDAD_AGRU == "14" ~ "POB_45_49",
          EDAD_AGRU == "15" ~ "POB_50_54",
          EDAD_AGRU == "16" ~ "POB_55_59",
          EDAD_AGRU == "17" ~ "POB_60_64",
          EDAD_AGRU == "18" ~ "POB_65_69",
          EDAD_AGRU == "19" ~ "POB_70_74",
          EDAD_AGRU == "20" ~ "POB_75_79",
          EDAD_AGRU == "21" ~ "POB_80_84",
          EDAD_AGRU %in% c("22","23","24","25","26","27","28","29") ~ "POB_85_mm",
          TRUE ~ NA_character_  # includes "30" no especificada
        )
      )
    
    df
  })
  
  names(out) <- as.character(anios)
  out
}

# ---- Main function: age-adjusted rates ----
tasa_mortalidad_ajustada_edad_total_y_cvd <- function(registros_defun_clave_mun, anios, poblacion_conapo_quinquenal) {
  
  niveles_edad <- c(
    "POB_00_04","POB_05_09","POB_10_14","POB_15_19","POB_20_24",
    "POB_25_29","POB_30_34","POB_35_39","POB_40_44","POB_45_49",
    "POB_50_54","POB_55_59","POB_60_64","POB_65_69","POB_70_74",
    "POB_75_79","POB_80_84","POB_85_mm"
  )
  
  # ---- Standard weights (Mexico 2010 from CONAPO quinquennial pop) ----
  pesos_2010 <- poblacion_conapo_quinquenal %>%
    filter(AÑO == 2010) %>%
    pivot_longer(cols = starts_with("POB_"), names_to = "grupo_edad", values_to = "poblacion") %>%
    filter(grupo_edad %in% niveles_edad) %>%
    group_by(grupo_edad) %>%
    summarise(pob_std = sum(poblacion, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE),
      peso = pob_std / sum(pob_std, na.rm = TRUE)
    ) %>%
    dplyr::select(grupo_edad, peso)
  
  ajuste_directo <- function(df) {
    df %>%
      left_join(pesos_2010, by = "grupo_edad") %>%
      mutate(tasa_esp = ifelse(poblacion > 0, defunciones / poblacion, NA_real_)) %>%
      summarise(tasa_ajustada = sum(peso * tasa_esp, na.rm = TRUE) * 1e5, .groups = "drop")
  }
  
  # ---- Add grupo_edad to each death record ----
  registros_edad_quinq <- clasificacion_edad_quinq(registros_defun_clave_mun, anios)
  
  # ---- Loop by year ----
  out <- purrr::map(setNames(as.list(anios), as.character(anios)), function(anio) {
    
    # Population long (municipal x age group)
    pob_year <- poblacion_conapo_quinquenal %>%
      filter(AÑO == anio) %>%
      mutate(
        cvegeo    = sprintf("%05d", suppressWarnings(as.integer(as.character(CLAVE)))),
        ENT_RESID = sprintf("%02d", suppressWarnings(as.integer(as.character(CLAVE_ENT)))),
        ANIO      = as.integer(AÑO)
      ) %>%
      pivot_longer(
        cols = starts_with("POB_"),
        names_to = "grupo_edad",
        values_to = "poblacion"
      ) %>%
      filter(grupo_edad %in% niveles_edad) %>%
      mutate(grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE)) %>%
      group_by(cvegeo, ENT_RESID, ANIO, grupo_edad) %>%
      summarise(poblacion = sum(poblacion, na.rm = TRUE), .groups = "drop")
    
    # Deaths aggregated (municipal x age group): total + CVD
    def_year <- registros_edad_quinq[[as.character(anio)]] %>%
      mutate(
        cvegeo    = sprintf("%05d", suppressWarnings(as.integer(as.character(cvegeo)))),
        ENT_RESID = substr(cvegeo, 1, 2),
        ANIO      = as.integer(anio)
      ) %>%
      filter(!is.na(grupo_edad)) %>%
      group_by(cvegeo, ENT_RESID, ANIO, grupo_edad) %>%
      summarise(
        def_total = n(),
        def_cvd   = sum(muerte_cvd, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Merge pop + deaths
    df_age_mun <- pob_year %>%
      left_join(def_year, by = c("cvegeo","ENT_RESID","ANIO","grupo_edad")) %>%
      mutate(
        def_total = coalesce(def_total, 0L),
        def_cvd   = coalesce(def_cvd, 0L)
      ) %>%
      filter(poblacion > 0)
    
    # ---- Municipal age-adjusted (total + CVD) ----
    mun_adj_total <- df_age_mun %>%
      mutate(defunciones = def_total) %>%
      group_by(cvegeo, ENT_RESID, ANIO) %>%
      group_modify(~ ajuste_directo(.x)) %>%
      ungroup() %>%
      rename(tasa_ajustada_edad_mun_total = tasa_ajustada)
    
    mun_adj_cvd <- df_age_mun %>%
      mutate(defunciones = def_cvd) %>%
      group_by(cvegeo, ENT_RESID, ANIO) %>%
      group_modify(~ ajuste_directo(.x)) %>%
      ungroup() %>%
      rename(tasa_ajustada_edad_mun_cvd = tasa_ajustada)
    
    # ---- Municipal crude (total + CVD) ----
    mun_cruda <- df_age_mun %>%
      group_by(cvegeo, ENT_RESID, ANIO) %>%
      summarise(
        def_mun_total = sum(def_total),
        def_mun_cvd   = sum(def_cvd),
        pob_mun_total = sum(poblacion),
        .groups = "drop"
      ) %>%
      mutate(
        tasa_cruda_total = ifelse(pob_mun_total > 0, def_mun_total / pob_mun_total * 1e5, NA_real_),
        tasa_cruda_cvd   = ifelse(pob_mun_total > 0, def_mun_cvd   / pob_mun_total * 1e5, NA_real_)
      )
    
    # ---- State age-adjusted (total + CVD) ----
    edo_adj_total <- df_age_mun %>%
      group_by(ENT_RESID, ANIO, grupo_edad) %>%
      summarise(defunciones = sum(def_total), poblacion = sum(poblacion), .groups = "drop") %>%
      group_by(ENT_RESID, ANIO) %>%
      group_modify(~ ajuste_directo(.x)) %>%
      ungroup() %>%
      rename(tasa_ajustada_edad_edo_total = tasa_ajustada)
    
    edo_adj_cvd <- df_age_mun %>%
      group_by(ENT_RESID, ANIO, grupo_edad) %>%
      summarise(defunciones = sum(def_cvd), poblacion = sum(poblacion), .groups = "drop") %>%
      group_by(ENT_RESID, ANIO) %>%
      group_modify(~ ajuste_directo(.x)) %>%
      ungroup() %>%
      rename(tasa_ajustada_edad_edo_cvd = tasa_ajustada)
    
    # ---- National age-adjusted (total + CVD) ----
    pais_adj_total <- df_age_mun %>%
      group_by(ANIO, grupo_edad) %>%
      summarise(defunciones = sum(def_total), poblacion = sum(poblacion), .groups = "drop") %>%
      group_by(ANIO) %>%
      group_modify(~ ajuste_directo(.x)) %>%
      ungroup() %>%
      rename(tasa_ajustada_pais_total = tasa_ajustada)
    
    pais_adj_cvd <- df_age_mun %>%
      group_by(ANIO, grupo_edad) %>%
      summarise(defunciones = sum(def_cvd), poblacion = sum(poblacion), .groups = "drop") %>%
      group_by(ANIO) %>%
      group_modify(~ ajuste_directo(.x)) %>%
      ungroup() %>%
      rename(tasa_ajustada_pais_cvd = tasa_ajustada)
    
    # ---- Join all yearly rates ----
    df_rates_year <- mun_cruda %>%
      left_join(mun_adj_total, by = c("cvegeo","ENT_RESID","ANIO")) %>%
      left_join(mun_adj_cvd,   by = c("cvegeo","ENT_RESID","ANIO")) %>%
      left_join(edo_adj_total, by = c("ENT_RESID","ANIO")) %>%
      left_join(edo_adj_cvd,   by = c("ENT_RESID","ANIO")) %>%
      left_join(pais_adj_total, by = "ANIO") %>%
      left_join(pais_adj_cvd,   by = "ANIO")
    
    df_rates_year
  })
  
  out
}

mortality_adjusted <- tasa_mortalidad_ajustada_edad_total_y_cvd(
  registros_defun_clave_mun = classified_deaths_id,
  anios = anios,
  poblacion_conapo_quinquenal = poblacion_conapo_quinquenal
)
mortality_adjusted[["2010"]] %>% head()


#Funcion para generar un resumen de las tablas
tabla_tasas_cvd <- function(mortality_adjusted, anios, usar_ponderacion = TRUE) {
  
  purrr::map_dfr(anios, function(anio) {
    
    df <- mortality_adjusted[[as.character(anio)]]
    if (is.null(df)) return(NULL)
    
    # 1 fila por municipio
    df_mun <- df %>%
      transmute(
        ANIO,
        cvegeo,
        def_mun_cvd,
        pob_mun_total,
        tasa_cruda_cvd,
        tasa_ajustada_edad_mun_cvd,
        tasa_ajustada_edad_edo_cvd,
        tasa_ajustada_pais_cvd
      ) %>%
      distinct(ANIO, cvegeo, .keep_all = TRUE)
    
    summarise(
      df_mun,
      ANIO = first(ANIO),
      
      # Cruda nacional CVD (a partir de suma municipal)
      tasa_cruda_pais_cvd =
        sum(def_mun_cvd, na.rm = TRUE) / sum(pob_mun_total, na.rm = TRUE) * 1e5,
      
      # Ajustada nacional CVD (viene constante por año dentro del df)
      tasa_ajustada_pais_cvd = first(na.omit(tasa_ajustada_pais_cvd)),
      
      # Promedio municipal ajustado (ponderado por población, opcional)
      tasa_ajustada_mun_media_cvd =
        if (usar_ponderacion) stats::weighted.mean(tasa_ajustada_edad_mun_cvd, w = pob_mun_total, na.rm = TRUE)
      else mean(tasa_ajustada_edad_mun_cvd, na.rm = TRUE),
      
      # Promedio estatal ajustado (replicado por municipios; ponderación por población municipal)
      tasa_ajustada_edo_media_cvd =
        if (usar_ponderacion) stats::weighted.mean(tasa_ajustada_edad_edo_cvd, w = pob_mun_total, na.rm = TRUE)
      else mean(tasa_ajustada_edad_edo_cvd, na.rm = TRUE),
      
      defunciones_cvd = sum(def_mun_cvd, na.rm = TRUE),
      
      .groups = "drop"
    )
  })
}

tabla_anual_cvd <- tabla_tasas_cvd(mortality_adjusted, anios = anios, usar_ponderacion = TRUE)
tabla_anual_cvd



#Plot
plot_tasas_en_tiempo_cvd <- function(df_ts) {
  
  df_ts <- df_ts %>% arrange(ANIO)
  
  ggplot(df_ts, aes(x = ANIO)) +
    geom_line(aes(y = tasa_cruda_pais_cvd, color = "Crude rate"), linewidth = 1) +
    geom_point(aes(y = tasa_cruda_pais_cvd, color = "Crude rate"), size = 2) +
    geom_line(aes(y = tasa_ajustada_pais_cvd, color = "Age-adjusted rate"), linewidth = 1) +
    geom_point(aes(y = tasa_ajustada_pais_cvd, color = "Age-adjusted rate"), size = 2) +
    scale_x_continuous(
      breaks = seq(min(df_ts$ANIO, na.rm = TRUE), max(df_ts$ANIO, na.rm = TRUE), by = 2)
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(
      title = "Cardiovascular Mortality Rates in Mexico",
      subtitle = "Crude vs age-adjusted (per 100,000 population)",
      x = "Year",
      y = "Rate per 100,000 population",
      color = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

p_cvd <- plot_tasas_en_tiempo_cvd(tabla_anual_cvd)
p_cvd

plot_tasa_ajustada_cvd <- function(df_ts) {
  
  df_ts <- df_ts %>% arrange(ANIO)
  
  p <- ggplot(df_ts, aes(x = ANIO, y = tasa_ajustada_pais_cvd)) +
    
    geom_line(linewidth = 1.1, color = "#1B4F72") +
    
    scale_x_continuous(
      breaks = seq(
        floor(min(df_ts$ANIO, na.rm = TRUE) / 5) * 5,
        ceiling(max(df_ts$ANIO, na.rm = TRUE) / 5) * 5,
        by = 5
      ),
      minor_breaks = seq(
        min(df_ts$ANIO, na.rm = TRUE),
        max(df_ts$ANIO, na.rm = TRUE),
        by = 1
      )
    ) +
    
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 8)
    ) +
    
    labs(
      x = "Year",
      y = "Age-adjusted rate per 100,000"
    ) +
    
    theme_classic(base_size = 13) +
    
    theme(
      axis.title = element_text(face = "bold"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(3, "pt"),
      
      # Activar solo líneas horizontales
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  return(p)
}

p_adj <- plot_tasa_ajustada_cvd(tabla_anual_cvd)
p_adj


# Exportar
ggsave(
  filename = file.path(ruta_figuras, "tasa_ajustada_cvd.svg"),
  plot = p_adj,
  device = "svg",
  width = 9, height = 5.5, units = "in"
)
ggsave(
  filename = file.path(ruta_figuras, "tasa_ajustada_cvd.png"),
  plot = p_adj,
  device = "png",
  width = 9, height = 5.5, units = "in",
  dpi = 600
)

#Tasas estratificadas por sexo
tasa_cvd_nacional_total_sexo <- function(registros_defun_clave_mun, anios, poblacion_conapo_quinquenal) {
  
  niveles_edad <- c(
    "POB_00_04","POB_05_09","POB_10_14","POB_15_19","POB_20_24",
    "POB_25_29","POB_30_34","POB_35_39","POB_40_44","POB_45_49",
    "POB_50_54","POB_55_59","POB_60_64","POB_65_69","POB_70_74",
    "POB_75_79","POB_80_84","POB_85_mm"
  )
  
  # Pesos estándar 2010 (sumando H+M)
  pesos_2010 <- poblacion_conapo_quinquenal %>%
    dplyr::filter(AÑO == 2010) %>%
    tidyr::pivot_longer(cols = starts_with("POB_"), names_to = "grupo_edad", values_to = "poblacion") %>%
    dplyr::filter(grupo_edad %in% niveles_edad) %>%
    dplyr::group_by(grupo_edad) %>%
    dplyr::summarise(pob_std = sum(poblacion, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE),
      peso = pob_std / sum(pob_std, na.rm = TRUE)
    ) %>%
    dplyr::select(grupo_edad, peso)
  
  ajuste_directo <- function(df) {
    df %>%
      dplyr::left_join(pesos_2010, by = "grupo_edad") %>%
      dplyr::mutate(tasa_esp = dplyr::if_else(poblacion > 0, defunciones / poblacion, NA_real_)) %>%
      dplyr::summarise(tasa_ajustada = sum(peso * tasa_esp, na.rm = TRUE) * 1e5, .groups = "drop")
  }
  
  registros_edad_quinq <- clasificacion_edad_quinq(registros_defun_clave_mun, anios)
  
  res <- purrr::map_dfr(as.character(anios), function(anio) {
    
    # Población nacional por edad x sexo
    pob_nat_sexo_age <- poblacion_conapo_quinquenal %>%
      dplyr::filter(AÑO == as.integer(anio)) %>%
      dplyr::mutate(
        ANIO = as.integer(AÑO),
        SEXO = dplyr::case_when(
          SEXO == "HOMBRES" ~ 1L,
          SEXO == "MUJERES" ~ 2L,
          TRUE ~ NA_integer_
        )
      ) %>%
      tidyr::pivot_longer(cols = starts_with("POB_"), names_to = "grupo_edad", values_to = "poblacion") %>%
      dplyr::filter(grupo_edad %in% niveles_edad, !is.na(SEXO)) %>%
      dplyr::mutate(grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE)) %>%
      dplyr::group_by(ANIO, SEXO, grupo_edad) %>%
      dplyr::summarise(poblacion = sum(poblacion, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(poblacion > 0)
    
    # Defunciones nacionales por edad x sexo (CVD)
    def_nat_sexo_age <- registros_edad_quinq[[as.character(anio)]] %>%
      dplyr::mutate(
        ANIO = as.integer(anio),
        SEXO = suppressWarnings(as.integer(as.character(SEXO)))
      ) %>%
      dplyr::filter(!is.na(grupo_edad), SEXO %in% c(1,2)) %>%
      dplyr::group_by(ANIO, SEXO, grupo_edad) %>%
      dplyr::summarise(def_cvd = sum(muerte_cvd, na.rm = TRUE), .groups = "drop")
    
    df_nat_sexo_age <- pob_nat_sexo_age %>%
      dplyr::left_join(def_nat_sexo_age, by = c("ANIO","SEXO","grupo_edad")) %>%
      dplyr::mutate(def_cvd = dplyr::coalesce(def_cvd, 0L))
    
    # Conteos por sexo
    conteos_sexo <- df_nat_sexo_age %>%
      dplyr::group_by(ANIO, SEXO) %>%
      dplyr::summarise(
        poblacion = sum(poblacion, na.rm = TRUE),
        defunciones_cvd = sum(def_cvd, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Conteos total
    conteos_total <- conteos_sexo %>%
      dplyr::group_by(ANIO) %>%
      dplyr::summarise(
        poblacion = sum(poblacion, na.rm = TRUE),
        defunciones_cvd = sum(defunciones_cvd, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(SEXO = 0L)
    
    conteos <- dplyr::bind_rows(conteos_total, conteos_sexo) %>%
      dplyr::mutate(SEXO = dplyr::recode(as.character(SEXO), `0`="TOTAL", `1`="HOMBRES", `2`="MUJERES"))
    
    # Tasas ajustadas por sexo
    tasas_sexo <- df_nat_sexo_age %>%
      dplyr::mutate(defunciones = def_cvd) %>%
      dplyr::group_by(ANIO, SEXO) %>%
      dplyr::group_modify(~ ajuste_directo(.x)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(tasa_ajustada_pais_cvd = tasa_ajustada) %>%
      dplyr::mutate(SEXO = dplyr::recode(as.character(SEXO), `1`="HOMBRES", `2`="MUJERES"))
    
    # Tasas ajustadas total (sumando sexos por edad)
    df_nat_total_age <- df_nat_sexo_age %>%
      dplyr::group_by(ANIO, grupo_edad) %>%
      dplyr::summarise(
        defunciones = sum(def_cvd, na.rm = TRUE),
        poblacion   = sum(poblacion, na.rm = TRUE),
        .groups = "drop"
      )
    
    tasas_total <- df_nat_total_age %>%
      dplyr::group_by(ANIO) %>%
      dplyr::group_modify(~ ajuste_directo(.x)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(tasa_ajustada_pais_cvd = tasa_ajustada) %>%
      dplyr::mutate(SEXO = "TOTAL")
    
    tasas <- dplyr::bind_rows(tasas_total, tasas_sexo)
    
    dplyr::left_join(tasas, conteos, by = c("ANIO","SEXO"))
  })
  
  list(
    tasas_nacional_total = res %>% dplyr::filter(SEXO=="TOTAL") %>% dplyr::select(ANIO, tasa_ajustada_pais_cvd) %>% dplyr::arrange(ANIO),
    tasas_nacional_sexo  = res %>% dplyr::filter(SEXO %in% c("HOMBRES","MUJERES")) %>% dplyr::select(ANIO, SEXO, tasa_ajustada_pais_cvd) %>% dplyr::arrange(ANIO, SEXO),
    conteos_nacional     = res %>% dplyr::select(ANIO, SEXO, poblacion, defunciones_cvd) %>% dplyr::distinct() %>% dplyr::arrange(ANIO, SEXO)
  )
}

res_nat <- tasa_cvd_nacional_total_sexo(
  registros_defun_clave_mun = classified_deaths_id,
  anios = anios,
  poblacion_conapo_quinquenal = poblacion_conapo_quinquenal
)

tasas_nacional_sexo  <- res_nat$tasas_nacional_sexo
defunciones_sexo     <- res_nat$conteos_nacional


# Juntar ambas tablas 
# 1) Asegurar nombre estándar de la columna de tasa (por si viene como tasa_ajustada_pais_cvd o tasa_ajustada_pais_cvd_sexo)
tasas_nacional_sexo2 <- tasas_nacional_sexo %>%
  rename(tasa_ajustada_pais_cvd = dplyr::any_of(c("tasa_ajustada_pais_cvd", "tasa_ajustada_pais_cvd_sexo", "tasa_ajustada_pais_cvd"))) %>%
  select(ANIO, SEXO, tasa_ajustada_pais_cvd)

# 2) Wide: columnas para HOMBRES y MUJERES
tasas_sexo_wide <- tasas_nacional_sexo2 %>%
  mutate(SEXO = toupper(SEXO)) %>%
  tidyr::pivot_wider(
    names_from  = SEXO,
    values_from = tasa_ajustada_pais_cvd,
    names_prefix = "tasa_ajustada_pais_cvd_"
  ) %>%
  # nombres más cortos opcional
  rename(
    tasa_ajustada_pais_cvd_hombres = tasa_ajustada_pais_cvd_HOMBRES,
    tasa_ajustada_pais_cvd_mujeres = tasa_ajustada_pais_cvd_MUJERES
  )

tasas_estratificadas <- tabla_anual_cvd %>%
  left_join(tasas_sexo_wide, by = "ANIO")

tasas_estratificadas


#Plot tasas estratificadas por sexo 
plot_tasas_ajustadas_sexo_cvd <- function(df_ts) {
  
  df_ts <- df_ts %>% arrange(ANIO)
  
  p <- ggplot(df_ts, aes(x = ANIO)) +
    
    geom_line(aes(y = tasa_ajustada_pais_cvd_hombres, color = "Men"),
              linewidth = 1.1) +
    
    geom_line(aes(y = tasa_ajustada_pais_cvd_mujeres, color = "Women"),
              linewidth = 1.1) +
    
    scale_color_manual(
      values = c("Men" = "#1B4F72",
                 "Women" = "#A93226")
    ) +
    
    scale_x_continuous(
      breaks = seq(
        floor(min(df_ts$ANIO, na.rm = TRUE) / 5) * 5,
        ceiling(max(df_ts$ANIO, na.rm = TRUE) / 5) * 5,
        by = 5
      ),
      minor_breaks = seq(
        min(df_ts$ANIO, na.rm = TRUE),
        max(df_ts$ANIO, na.rm = TRUE),
        by = 1
      )
    ) +
    
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 8)
    ) +
    
    labs(
      x = "Year",
      y = "Age-adjusted rate per 100,000",
      color = NULL
    ) +
    
    theme_classic(base_size = 13) +
    
    theme(
      axis.title = element_text(face = "bold"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(3, "pt"),
      
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      legend.position = "top"
    )
  
  return(p)
}

p_sexo <- plot_tasas_ajustadas_sexo_cvd(tasas_estratificadas)
p_sexo

#Exportar
ggsave(
  filename = file.path(ruta_figuras, "Figure1A.svg"),
  plot = p_sexo,
  device = "svg",
  width = 9, height = 5.5, units = "in"
)
ggsave(
  filename = file.path(ruta_figuras, "Figure1A.png"),
  plot = p_sexo,
  device = "png",
  width = 9, height = 5.5, units = "in",
  dpi = 600
)

#Tasas ajustadas por macroregion
mortalidad_macroregion_age_adj_cvd <- function(anios, registros_defun_clave_mun, poblacion_conapo_quinquenal) {
  
  asignar_macro <- function(ent) {
    dplyr::case_when(
      ent %in% c(2, 3, 5, 8, 19, 25, 26, 28) ~ 1L,  # North
      ent %in% c(1, 6, 10, 11, 14, 16, 18, 24, 32) ~ 2L,  # Central-West
      ent %in% c(9, 13, 15, 17, 21, 22, 29) ~ 3L,  # Central
      ent %in% c(4, 7, 12, 20, 23, 27, 30, 31) ~ 4L,  # South-Southeast
      TRUE ~ NA_integer_
    )
  }
  
  niveles_edad <- c(
    "POB_00_04","POB_05_09","POB_10_14","POB_15_19","POB_20_24",
    "POB_25_29","POB_30_34","POB_35_39","POB_40_44","POB_45_49",
    "POB_50_54","POB_55_59","POB_60_64","POB_65_69","POB_70_74",
    "POB_75_79","POB_80_84","POB_85_mm"
  )
  
  # Pesos estándar 2010 (nacional) sumando todo
  pesos_2010 <- poblacion_conapo_quinquenal %>%
    dplyr::filter(AÑO == 2010) %>%
    tidyr::pivot_longer(cols = starts_with("POB_"), names_to = "grupo_edad", values_to = "poblacion") %>%
    dplyr::filter(grupo_edad %in% niveles_edad) %>%
    dplyr::group_by(grupo_edad) %>%
    dplyr::summarise(pob_std = sum(poblacion, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE),
      peso = pob_std / sum(pob_std, na.rm = TRUE)
    ) %>%
    dplyr::select(grupo_edad, peso)
  
  ajuste_directo <- function(df) {
    df %>%
      dplyr::left_join(pesos_2010, by = "grupo_edad") %>%
      dplyr::mutate(tasa_esp = dplyr::if_else(poblacion > 0, defunciones / poblacion, NA_real_)) %>%
      dplyr::summarise(tasa_ajustada = sum(peso * tasa_esp, na.rm = TRUE) * 1e5, .groups = "drop")
  }
  
  # 1) Mapear edades en defunciones a quinquenios CONAPO
  registros_edad_quinq <- clasificacion_edad_quinq(registros_defun_clave_mun, anios)
  
  purrr::map_dfr(anios, function(anio) {
    
    df <- registros_edad_quinq[[as.character(anio)]]
    if (is.null(df)) return(NULL)
    
    # Defunciones CVD por macroregion x edad
    def_macro_age <- df %>%
      dplyr::mutate(
        anio = as.integer(anio),
        ent  = suppressWarnings(as.integer(as.character(ENT_OCURR))),
        macro_cod = asignar_macro(ent)
      ) %>%
      dplyr::filter(!is.na(macro_cod), !is.na(grupo_edad)) %>%
      dplyr::group_by(anio, macro_cod, grupo_edad) %>%
      dplyr::summarise(defunciones = sum(muerte_cvd, na.rm = TRUE), .groups = "drop")
    
    # Población por macroregion x edad (agregando entidades)
    pob_macro_age <- poblacion_conapo_quinquenal %>%
      dplyr::filter(AÑO == as.integer(anio)) %>%
      dplyr::mutate(
        anio = as.integer(AÑO),
        ent  = suppressWarnings(as.integer(as.character(CLAVE_ENT))),
        macro_cod = asignar_macro(ent)
      ) %>%
      dplyr::filter(!is.na(macro_cod)) %>%
      tidyr::pivot_longer(cols = starts_with("POB_"), names_to = "grupo_edad", values_to = "poblacion") %>%
      dplyr::filter(grupo_edad %in% niveles_edad) %>%
      dplyr::mutate(grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE)) %>%
      dplyr::group_by(anio, macro_cod, grupo_edad) %>%
      dplyr::summarise(poblacion = sum(poblacion, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(poblacion > 0)
    
    # Merge + tasas crudas por macro
    df_macro_age <- pob_macro_age %>%
      dplyr::left_join(def_macro_age, by = c("anio","macro_cod","grupo_edad")) %>%
      dplyr::mutate(defunciones = dplyr::coalesce(defunciones, 0L))
    
    crude_macro <- df_macro_age %>%
      dplyr::group_by(anio, macro_cod) %>%
      dplyr::summarise(
        defunciones_cvd = sum(defunciones, na.rm = TRUE),
        poblacion = sum(poblacion, na.rm = TRUE),
        tasa_cvd_cruda = defunciones_cvd / poblacion * 1e5,
        .groups = "drop"
      )
    
    # Ajustada por edad por macro
    adj_macro <- df_macro_age %>%
      dplyr::group_by(anio, macro_cod) %>%
      dplyr::group_modify(~ ajuste_directo(.x)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(tasa_cvd_ajustada = tasa_ajustada)
    
    dplyr::left_join(crude_macro, adj_macro, by = c("anio","macro_cod"))
  })
}

macro_rates <- mortalidad_macroregion_age_adj_cvd(
  anios = anios,
  registros_defun_clave_mun = classified_deaths_id,
  poblacion_conapo_quinquenal = poblacion_conapo_quinquenal
)

macro_rates_plot <- macro_rates %>%
  dplyr::mutate(
    macroregion = factor(
      macro_cod,
      levels = c(1,2,3,4),
      labels = c("North","Central-West","Central","South-Southeast")
    )
  )

p_macro <- ggplot(macro_rates_plot,
                  aes(x = anio,
                      y = tasa_cvd_ajustada,
                      group = macroregion,
                      color = macroregion)) +
  
  geom_line(linewidth = 1.1) +
  
  scale_x_continuous(
    breaks = seq(
      floor(min(macro_rates_plot$anio, na.rm = TRUE) / 5) * 5,
      ceiling(max(macro_rates_plot$anio, na.rm = TRUE) / 5) * 5,
      by = 5
    ),
    minor_breaks = seq(
      min(macro_rates_plot$anio, na.rm = TRUE),
      max(macro_rates_plot$anio, na.rm = TRUE),
      by = 1
    )
  ) +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 8)
  ) +
  
  labs(
    x = "Year",
    y = "Age-adjusted CVD mortality rate (per 100,000)",
    color = "Macro-region"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(3, "pt"),
    
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    legend.position = "right"
  )

p_macro

#Exportar
ggsave(
  filename = file.path(ruta_figuras, "Figure1B.svg"),
  plot = p_macro,
  device = "svg",
  width = 9, height = 5.5, units = "in"
)
ggsave(
  filename = file.path(ruta_figuras, "Figure1B.png"),
  plot = p_macro,
  device = "png",
  width = 9, height = 5.5, units = "in",
  dpi = 600
)


#Integrar a tasas estratificadas
tasas_macro_wide <- macro_rates %>%
  transmute(
    ANIO = as.integer(anio),
    macro = dplyr::case_when(
      macro_cod == 1L ~ "north",
      macro_cod == 2L ~ "centralwest",
      macro_cod == 3L ~ "central",
      macro_cod == 4L ~ "southsoutheast",
      TRUE ~ NA_character_
    ),
    tasa = tasa_cvd_ajustada
  ) %>%
  filter(!is.na(macro)) %>%
  tidyr::pivot_wider(
    names_from = macro,
    values_from = tasa,
    names_prefix = "tasa_adj_cvd_"
  ) %>%
  arrange(ANIO)

tasas_estratificadas <- tasas_estratificadas %>%
  left_join(tasas_macro_wide, by = "ANIO")

#---- Generar rutas y cargar IRS-----
generar_rutas_irs <- function(ruta_base) {
  
  ruta_base <- normalizePath(ruta_base, winslash = "/", mustWork = TRUE)
  
  carpeta <- list.files(ruta_base, pattern = "IRS_ent_mun", full.names = TRUE)
  
  if (!dir.exists(carpeta)) {
    warning("No se encontró la carpeta ", carpeta)
    return(NULL)
  }
  
  archivos <- list.files(carpeta, full.names = TRUE)
  
  if (length(archivos) == 0) {
    warning("No se encontró ningun archivo dentro de la carpeta ", carpeta)
    return(NULL)
  }
  
  ruta <- lapply(archivos, function(irs) {
    
    irs_pre <- readxl::read_excel(irs, sheet = "Municipios", skip = 2)
    
    irs_post <- irs_pre[!is.na(irs_pre[[1]]), ]
    
    irs_post
  })
  
  names(ruta) <- gsub("[^0-9]", "", basename(archivos))
  return(ruta)
}

irs__mun <- generar_rutas_irs (ruta_base = ruta_base)

#---- Ajustar df de IRS -----
renombrar_columnas_irs <- function(lista_irs) {
  
  nuevos_6_16 <- c(
    "Población de 15 años o más analfabeta",
    "Población de 6 a 14 años que no asiste a la escuela",
    "Población de 15 años y más con educación básica incompleta",
    "Población sin derechohabiencia a servicios de salud",
    "Viviendas con piso de tierra",
    "Viviendas que no disponen de excusado o sanitario",
    "Viviendas que no disponen de agua entubada de la red pública",
    "Viviendas que no disponen de drenaje",
    "Viviendas que no disponen de energía eléctrica",
    "Viviendas que no disponen de lavadora",
    "Viviendas que no disponen de refrigerador")
  
  lista_ajustada <- lapply(lista_irs, function(df) {
    
    nombres <- colnames(df)
    
    if (length(nombres) >= 16) {
      nombres[6:16] <- nuevos_6_16
    }
    
    colnames(df) <- nombres
    df
  })
  
  names(lista_ajustada) <- names(lista_irs)
  
  return(lista_ajustada)
}

irs__mun <- renombrar_columnas_irs(irs__mun)

# ---- Tasas ajustadas por IRS ----
# asignar año IRS más cercano (>=2000)
map_anio_irs_cercano <- function(y, irs_years){
  y <- as.integer(y)
  irs_years <- sort(as.integer(irs_years))
  if (is.na(y) || y < 2000) return(NA_integer_)
  irs_years[which.min(abs(irs_years - y))]
}

# preparar IRS: lista -> df largo (anio_irs, cvegeo, estrato_irs)
prep_irs_long <- function(irs__mun){
  
  irs_long <- purrr::imap_dfr(irs__mun, function(df, yr){
    
    df2 <- tibble::as_tibble(df)
    
    # Clave municipal (5 dígitos; puede tener 0 inicial)
    if (!"Clave municipio" %in% names(df2)) {
      stop("No encontré 'Clave municipio' en IRS año ", yr,
           ". Columnas: ", paste(names(df2), collapse = ", "))
    }
    
    # Estrato: usar "Grado de rezago social" (categoría)
    estr_col <- c("Grado de rezago social", "Grado de rezago")[
      c("Grado de rezago social", "Grado de rezago") %in% names(df2)
    ][1]
    
    if (is.na(estr_col)) {
      stop("No encontré 'Grado de rezago social' en IRS año ", yr,
           ". Columnas: ", paste(names(df2), collapse = ", "))
    }
    
    df2 %>%
      dplyr::transmute(
        anio_irs = as.integer(yr),
        cvegeo = sprintf("%05d", suppressWarnings(as.integer(`Clave municipio`))),
        estrato_irs = as.character(.data[[estr_col]])
      ) %>%
      dplyr::filter(!is.na(cvegeo), !is.na(estrato_irs), estrato_irs != "")
  })
  
  irs_long
}

# Función principal: tasas ajustadas por edad estratificadas por estrato IRS
tasa_cvd_adj_por_irs <- function(anios,
                                 registros_defun_clave_mun,
                                 poblacion_conapo_quinquenal,
                                 irs__mun){
  
  niveles_edad <- c(
    "POB_00_04","POB_05_09","POB_10_14","POB_15_19","POB_20_24",
    "POB_25_29","POB_30_34","POB_35_39","POB_40_44","POB_45_49",
    "POB_50_54","POB_55_59","POB_60_64","POB_65_69","POB_70_74",
    "POB_75_79","POB_80_84","POB_85_mm"
  )
  
  # Pesos estándar 2010 (nacional)
  pesos_2010 <- poblacion_conapo_quinquenal %>%
    dplyr::filter(AÑO == 2010) %>%
    tidyr::pivot_longer(cols = starts_with("POB_"), names_to = "grupo_edad", values_to = "poblacion") %>%
    dplyr::filter(grupo_edad %in% niveles_edad) %>%
    dplyr::group_by(grupo_edad) %>%
    dplyr::summarise(pob_std = sum(poblacion, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE),
      peso = pob_std / sum(pob_std, na.rm = TRUE)
    ) %>%
    dplyr::select(grupo_edad, peso)
  
  ajuste_directo <- function(df){
    
    df2 <- df %>%
      dplyr::mutate(
        grupo_edad = factor(
          as.character(grupo_edad),
          levels = niveles_edad,
          ordered = TRUE
        )
      ) %>%
      dplyr::left_join(pesos_2010, by = "grupo_edad") %>%
      dplyr::mutate(
        tasa_esp = if_else(poblacion > 0, defunciones / poblacion, 0)
      )
    
    # DEBUG opcional:
    # print(sum(df2$peso, na.rm=TRUE))
    
    tibble(
      tasa_ajustada = sum(df2$peso * df2$tasa_esp, na.rm = TRUE) * 1e5
    )
  }
  
  # IRS long
  irs_long <- prep_irs_long(irs__mun)
  irs_years <- sort(unique(irs_long$anio_irs))
  
  # Defunciones con grupo_edad
  registros_edad_quinq <- clasificacion_edad_quinq(registros_defun_clave_mun, anios)
  
  purrr::map_dfr(anios, function(anio){
    
    anio <- as.integer(anio)
    anio_irs <- map_anio_irs_cercano(anio, irs_years)
    if (is.na(anio_irs)) return(NULL)  # 1999 y anteriores fuera
    
    # IRS del año asignado
    irs_y <- irs_long %>%
      dplyr::filter(anio_irs == !!anio_irs) %>%
      dplyr::select(cvegeo, estrato_irs)
    
    # Población municipal por edad (año)  <<< AQUÍ VA EL CAMBIO CLAVE >>>
    # Construimos cvegeo (5 dígitos) desde entidad + municipio (municipio puede venir en 4 dígitos)
    # Población municipal por edad (año)
    pob_year <- poblacion_conapo_quinquenal %>%
      dplyr::filter(AÑO == anio) %>%
      dplyr::mutate(
        # CLAVE es el cvegeo municipal pero a veces viene sin 0 inicial -> estandarizar a 5 dígitos
        cvegeo = sprintf("%05d", suppressWarnings(as.integer(as.character(CLAVE)))),
        ANIO   = as.integer(AÑO)
      ) %>%
      tidyr::pivot_longer(
        cols = starts_with("POB_"),
        names_to = "grupo_edad",
        values_to = "poblacion"
      ) %>%
      dplyr::filter(grupo_edad %in% niveles_edad) %>%
      dplyr::mutate(grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE)) %>%
      dplyr::group_by(cvegeo, ANIO, grupo_edad) %>%
      dplyr::summarise(poblacion = sum(poblacion, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(poblacion > 0) %>%
      dplyr::left_join(irs_y, by = "cvegeo") %>%
      dplyr::filter(!is.na(estrato_irs))
    
    # Defunciones municipal por edad (año)
    def_year <- registros_edad_quinq[[as.character(anio)]] %>%
      dplyr::mutate(
        cvegeo = sprintf("%05d", suppressWarnings(as.integer(as.character(cvegeo)))),
        ANIO   = as.integer(anio)
      ) %>%
      dplyr::filter(!is.na(grupo_edad)) %>%
      dplyr::group_by(cvegeo, ANIO, grupo_edad) %>%
      dplyr::summarise(defunciones = sum(muerte_cvd, na.rm = TRUE), .groups = "drop") %>%
      dplyr::left_join(irs_y, by = "cvegeo") %>%
      dplyr::filter(!is.na(estrato_irs))
    
    # Merge edad-específico por estrato IRS
    df_age_irs <- pob_year %>%
      dplyr::left_join(def_year, by = c("cvegeo","ANIO","grupo_edad","estrato_irs")) %>%
      dplyr::mutate(defunciones = dplyr::coalesce(defunciones, 0L))
    
    # Conteos por estrato
    conteos <- df_age_irs %>%
      dplyr::group_by(ANIO, estrato_irs) %>%
      dplyr::summarise(
        poblacion = sum(poblacion, na.rm = TRUE),
        defunciones_cvd = sum(defunciones, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Ajuste directo por estrato
    # 1) Agregar a nivel estrato x edad (CRÍTICO)
    df_strato_age <- df_age_irs %>%
      dplyr::group_by(ANIO, estrato_irs, grupo_edad) %>%
      dplyr::summarise(
        defunciones = sum(defunciones, na.rm = TRUE),
        poblacion   = sum(poblacion, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 2) Ajuste directo por estrato (ya con 1 fila por edad)
    tasas_adj <- df_strato_age %>%
      dplyr::group_by(ANIO, estrato_irs) %>%
      dplyr::group_modify(~ ajuste_directo(.x)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(tasa_ajustada = tasa_ajustada)
    
    # 3) Unir tasas + conteos y devolver
    dplyr::left_join(tasas_adj, conteos, by = c("ANIO", "estrato_irs")) %>%
      dplyr::mutate(anio_irs_asignado = anio_irs)
  })
}
# ---- Ejecutar ----
tasas_adj_irs <- tasa_cvd_adj_por_irs(
  anios = anios,
  registros_defun_clave_mun = classified_deaths_id,
  poblacion_conapo_quinquenal = poblacion_conapo_quinquenal,
  irs__mun = irs__mun)

summary(tasas_adj_irs$tasa_ajustada)

tasas_adj_irs_wide <- tasas_adj_irs %>%
  dplyr::select(ANIO, estrato_irs, tasa_ajustada) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from = estrato_irs,
    values_from = tasa_ajustada,
    names_prefix = "tasa_irs_"
  ) %>%
  dplyr::arrange(ANIO)

#Integrar a tasas_estratificadas
tasas_estratificadas <- tasas_estratificadas %>%
  dplyr::left_join(tasas_adj_irs_wide, by = "ANIO")

#Crear tasas estratificadas por grupos grandes de edad
clasificar_grupo_amplio <- function(grupo_edad){
  
  dplyr::case_when(
    grupo_edad %in% c("POB_00_04","POB_05_09","POB_10_14",
                      "POB_15_19","POB_20_24") ~ "<25",
    
    grupo_edad %in% c("POB_25_29","POB_30_34","POB_35_39",
                      "POB_40_44") ~ "25-45",
    
    grupo_edad %in% c("POB_45_49","POB_50_54","POB_55_59",
                      "POB_60_64") ~ "45-65",
    
    grupo_edad %in% c("POB_65_69","POB_70_74","POB_75_79",
                      "POB_80_84","POB_85_mm") ~ ">65",
    
    TRUE ~ NA_character_
  )
}

tasas_grupo_edad <- purrr::map_dfr(anios, function(anio){
  
  # Defunciones por edad
  df_def <- clasificacion_edad_quinq(classified_deaths_id, anios)[[as.character(anio)]] %>%
    dplyr::filter(!is.na(grupo_edad)) %>%
    dplyr::mutate(grupo_amplio = clasificar_grupo_amplio(grupo_edad)) %>%
    dplyr::group_by(grupo_amplio) %>%
    dplyr::summarise(defunciones = sum(muerte_cvd, na.rm = TRUE), .groups="drop")
  
  # Población por edad
  df_pob <- poblacion_conapo_quinquenal %>%
    dplyr::filter(AÑO == as.integer(anio)) %>%
    tidyr::pivot_longer(cols = starts_with("POB_"),
                        names_to = "grupo_edad",
                        values_to = "poblacion") %>%
    dplyr::mutate(grupo_amplio = clasificar_grupo_amplio(grupo_edad)) %>%
    dplyr::group_by(grupo_amplio) %>%
    dplyr::summarise(poblacion = sum(poblacion, na.rm = TRUE), .groups="drop")
  
  df_pob %>%
    dplyr::left_join(df_def, by="grupo_amplio") %>%
    dplyr::mutate(
      ANIO = as.integer(anio),
      defunciones = dplyr::coalesce(defunciones,0),
      tasa_cvd = defunciones / poblacion * 1e5
    )
})

tasas_grupo_edad_wide <- tasas_grupo_edad %>%
  dplyr::select(ANIO, grupo_amplio, tasa_cvd) %>%
  tidyr::pivot_wider(
    names_from = grupo_amplio,
    values_from = tasa_cvd,
    names_prefix = "tasa_cvd_edad_"
  )


#Integrar en tasas estratificadas
tasas_estratificadas <- tasas_estratificadas %>%
  dplyr::left_join(tasas_grupo_edad_wide, by="ANIO")


#Crear tabla por quinquenios 
tasas_estratificadas_quinquenios <- tasas_estratificadas %>%
  mutate(
    quinquenio = case_when(
      ANIO >= 1990 & ANIO <= 1994 ~ "1990-1994",
      ANIO >= 1995 & ANIO <= 1999 ~ "1995-1999",
      ANIO >= 2000 & ANIO <= 2004 ~ "2000-2004",
      ANIO >= 2005 & ANIO <= 2009 ~ "2005-2009",
      ANIO >= 2010 & ANIO <= 2014 ~ "2010-2014",
      ANIO >= 2015 & ANIO <= 2019 ~ "2015-2019",
      ANIO >= 2020 & ANIO <= 2024 ~ "2020-2024",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(quinquenio)) %>%
  group_by(quinquenio) %>%
  summarise(
    across(
      -ANIO,
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

#Crear flex table
library(flextable)
library(officer)
ft <- flextable(tasas_estratificadas_quinquenios)

ft <- ft %>%
  autofit() %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header")

ft

#Exportar a word
Tabla_tasas_quinquenios <- read_docx() %>%
  body_add_par(
    "Age-standardized cardiovascular mortality rates by sex, 1990–2024 (5-year periods)",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft)

archivo_salida <- file.path(
  ruta_figuras,
  "Table1_adjusted_rates_quinquenios.docx"
)

# Guardar
print(Tabla_tasas_quinquenios, target = archivo_salida)


# ---- Crear incremento por año
tasas_estratificadas_quinquenios <- tasas_estratificadas_quinquenios %>%
  dplyr::mutate(
    quinquenio_dummy = as.integer(factor(quinquenio, levels = sort(unique(quinquenio))))
  )

#Hombres
incremento_hombres <-MASS:: glm.nb(tasa_ajustada_pais_cvd_hombres~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_hombres, exp=T)

#Mujeres
incremento_mujeres <-MASS:: glm.nb(tasa_ajustada_pais_cvd_mujeres~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_mujeres, exp=T)

#Norte
incremento_norte <- m_hombres<-MASS:: glm.nb(tasa_adj_cvd_north~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_norte, exp=T)

#Centro-oeste
incremento_co<- m_hombres<-MASS:: glm.nb(tasa_adj_cvd_centralwest~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_hombres, exp=T)

#Central
incremento_c<- m_hombres<-MASS:: glm.nb(tasa_adj_cvd_central~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_c, exp=T)

#South east
incremento_se<- m_hombres<-MASS:: glm.nb(tasa_adj_cvd_southsoutheast~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_se, exp=T)

#<25 years
incremento_25<- m_hombres<-MASS:: glm.nb(tasa_cvd_edad_less25~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_25, exp=T)

#25-45
incremento_25_45<- m_hombres<-MASS:: glm.nb(tasa_cvd_edad_25_45~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_25_45, exp=T)


#45-65
incremento_45_65<- m_hombres<-MASS:: glm.nb(tasa_cvd_edad_45_65~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_45_65, exp=T)

#>65 years
incremento_65<- m_hombres<-MASS:: glm.nb(tasa_cvd_edad_65~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_65, exp=T)

#Very high SLI
incremento_very_high<- m_hombres<-MASS:: glm.nb(tasa_muy_alto~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_very_high, exp=T)

#High SLI
incremento_high<- m_hombres<-MASS:: glm.nb(tasa_irs_Alto~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_high, exp=T)

#Medium SLI
incremento_medium<- m_hombres<-MASS:: glm.nb(tasa_irs_Medio~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_medium, exp=T)

#Low SLI
incremento_low<- m_hombres<-MASS:: glm.nb(tasa_irs_Bajo~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_low, exp=T)

#Very low SLI
incremento_ver_low<- m_hombres<-MASS:: glm.nb(tasa_muy_bajo~quinquenio_dummy, data = tasas_estratificadas_quinquenios)
jtools::summ(incremento_ver_low, exp=T)

#Change names
tasas_estratificadas_quinquenios <- tasas_estratificadas_quinquenios %>%
  rename(tasa_muy_bajo = "tasa_irs_Muy bajo")
  
class(tasas_estratificadas_quinquenios$tasa_cvd_edad_65)
  


#Por quinquenio
tasas_estratificadas_quinquenios <- tasas_estratificadas_quinquenios %>%
  dplyr::mutate(
    year_dummy = as.integer(factor(quinquenio, levels = sort(unique(quinquenio))))
  )

modelo_quinquenio <- glm(
  log(tasa_ajustada_pais_cvd) ~ year_dummy,
  data = tasas_estratificadas_quinquenios
)
summary(modelo_quinquenio)

exp(coef(modelo_quinquenio)["year_dummy"]) - 1

#Comparacion periodos
tasas_estratificadas_quinquenios %>%
  mutate(periodo = factor(quinquenio))

first <- tasas_estratificadas_quinquenios$tasa_ajustada_pais_cvd[1]
last  <- tail(tasas_estratificadas_quinquenios$tasa_ajustada_pais_cvd,1)

((last/first) - 1) * 100


# ---- Crear DF de solo muertes CV ----
classified_deaths_id_cvd <- purrr::map(
  classified_deaths_id,
  ~ dplyr::filter(.x, muerte_cvd == 1)
)


####### Funciones para estandarización de variables ######
#----Educacion----
recodificacion_educacion <- function(registros_defun_clave_mun, anios){
  
  df_anio <- lapply(anios, function(anio){
    
    regis_anio <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(regis_anio)) return(NULL)
    
    if (!("ESCOLARIDA" %in% names(regis_anio))) {
      warning("No existe la columna ESCOLARIDA en el año: ", anio)
      regis_anio$ANIO <- as.integer(anio)
      regis_anio$NIVEL_ESCOLARIDAD <- NA_character_
      return(regis_anio)
    }
    
    esc <- suppressWarnings(as.integer(as.character(regis_anio$ESCOLARIDA)))
    
    df_final <- regis_anio %>%
      dplyr::mutate(
        ANIO = as.integer(anio),
        ESCOLARIDA_int = esc,
        NIVEL_ESCOLARIDAD = dplyr::case_when(
          
          ANIO <= 2011 & ESCOLARIDA_int %in% c(1, 2, 3, 4, 5) ~ "low_education",
          ANIO <= 2011 & ESCOLARIDA_int == 6               ~ "medium_education",
          ANIO <= 2011 & ESCOLARIDA_int == 7        ~ "high_education",
          ANIO <= 2011 & ESCOLARIDA_int %in% c(8, 9)       ~ NA_character_,
          
          ANIO >= 2012 & ESCOLARIDA_int %in% c(1, 2, 3, 4, 5,6) ~ "low_education",
          ANIO >= 2012 & ESCOLARIDA_int %in% c(7, 8)         ~ "medium_education",
          ANIO >= 2012 & ESCOLARIDA_int %in% c(9, 10)     ~ "high_education",
          ANIO >= 2012 & ESCOLARIDA_int %in% c(88, 99)       ~ NA_character_,
          
          TRUE ~ NA_character_
        )
      )
    
    df_final
  })
  
  names(df_anio) <- as.character(anios)
  df_anio
}

edu_classified_deaths_id <- recodificacion_educacion(classified_deaths_id, anios)
table(edu_classified_deaths_id[["2024"]]$NIVEL_ESCOLARIDAD, useNA = "ifany")


#----Derechohabiencia----
recodificacion_derechohabiencia <- function(registros_defun_clave_mun, anios){
  
  out <- lapply(anios, function(anio){
    
    regis <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(regis)) return(NULL)
    
    if (!("DERECHOHAB" %in% names(regis))) {
      warning("No existe la columna DERECHOHAB en el año: ", anio)
      regis$ANIO <- as.integer(anio)
      regis$seguro <- NA_integer_
      return(regis)
    }
    
    dh <- suppressWarnings(as.integer(as.character(regis$DERECHOHAB)))
    
    regis %>%
      dplyr::mutate(
        ANIO = as.integer(anio),
        DERECHOHAB_int = dh,
        seguro = dplyr::case_when(
          
          # <= 2003
          ANIO <= 2003 & DERECHOHAB_int %in% c(1,6) ~ 0L,
          ANIO <= 2003 & DERECHOHAB_int %in% c(2,3,4,5,23,24,25,26,34,35,36,45,46,56) ~ 1L,
          ANIO <= 2003 & DERECHOHAB_int == 99 ~ NA_integer_,
          
          # 2004–2011
          ANIO >= 2004 & ANIO <= 2011 & DERECHOHAB_int %in% c(1,8) ~ 0L,
          ANIO >= 2004 & ANIO <= 2011 & DERECHOHAB_int %in% c(2,3,4,5,6,7) ~ 1L,
          ANIO >= 2004 & ANIO <= 2011 & DERECHOHAB_int == 99 ~ NA_integer_,
          
          # 2012–2017
          ANIO >= 2012 & ANIO <= 2017 & DERECHOHAB_int %in% c(1,8) ~ 0L,
          ANIO >= 2012 & ANIO <= 2017 & DERECHOHAB_int %in% c(2,3,4,5,6,7,9) ~ 1L,
          ANIO >= 2012 & ANIO <= 2017 & DERECHOHAB_int == 99 ~ NA_integer_,
          
          # 2018–2021 
          ANIO >= 2018 & ANIO <= 2021 & DERECHOHAB_int %in% c(1,8) ~ 0L,
          ANIO >= 2018 & ANIO <= 2021 & DERECHOHAB_int %in% c(2,3,4,5,6,7,9) ~ 1L,
          ANIO >= 2018 & ANIO <= 2021 & DERECHOHAB_int == 99 ~ NA_integer_,
          
          # 2022-2023
          ANIO >= 2022 & ANIO <= 2023 & DERECHOHAB_int %in% c(1,8) ~ 0L,
          ANIO >= 2022 & ANIO <= 2023 & DERECHOHAB_int %in% c(2,3,4,5,6,7,9,10,12) ~ 1L,
          ANIO >= 2022 & ANIO <= 2023 & DERECHOHAB_int == 99 ~ NA_integer_,
          
          # 2024
          ANIO == 2024 & DERECHOHAB_int %in% c(1,8) ~ 0L,
          ANIO == 2024 & DERECHOHAB_int %in% c(2,3,4,5,6,7,9,10,12) ~ 1L,
          ANIO == 2024 & DERECHOHAB_int == 99 ~ NA_integer_,
          
          
          TRUE ~ NA_integer_
        )
      )
  })
  
  names(out) <- as.character(anios)
  out
}

ss_classified_deaths_id <- recodificacion_derechohabiencia(classified_deaths_id, anios)
df2010 <- ss_classified_deaths_id[["2010"]]
table(df2010$seguro, useNA = "ifany")

#----Lengua indigena----
recodificacion_indigena <- function(registros_defun_clave_mun, anios){
  
  df_anio <- lapply(anios, function(anio){
    
    regis_anio <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(regis_anio)) return(NULL)
    
    regis_anio <- regis_anio %>% dplyr::mutate(ANIO = as.integer(anio))
    
    # Si no existe LENGUA (p.ej. <2012), crear NA
    if (!("LENGUA" %in% names(regis_anio))) {
      regis_anio$LENGUA_INDIGENA <- NA_integer_
      return(regis_anio)
    }
    
    lengua_int <- suppressWarnings(as.integer(as.character(regis_anio$LENGUA)))
    
    regis_anio %>%
      dplyr::mutate(
        LENGUA_int = lengua_int,
        LENGUA_INDIGENA = dplyr::case_when(
          LENGUA_int == 2 ~ 0L,
          LENGUA_int == 1 ~ 1L,
          LENGUA_int %in% c(9, 99) ~ NA_integer_,
          TRUE ~ NA_integer_
        )
      )
  })
  
  names(df_anio) <- as.character(anios)
  df_anio
}

indigenous_classified_deaths_id <- recodificacion_indigena(classified_deaths_id, anios)
table(indigenous_classified_deaths_id[["2012"]]$LENGUA_INDIGENA, useNA="ifany")


table(all_deaths$`2022`$OCUPACION)
#----Ocupacion----
recodificacion_ocupacion <- function(registros_defun_clave_mun, anios){
  
  out <- lapply(anios, function(anio){
    
    regis <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(regis)) return(NULL)
    
    regis <- regis %>% dplyr::mutate(ANIO = as.integer(anio))
    
    if (!("OCUPACION" %in% names(regis))) {
      warning("No existe la columna OCUPACION en el año: ", anio)
      regis$NIVEL_OCUPACION <- NA_character_
      return(regis)
    }
    
    # preparar ocupación en dos formatos
    regis <- regis %>%
      dplyr::mutate(
        OCUP_INT = suppressWarnings(as.integer(as.character(OCUPACION))),
        OCUP_3   = stringr::str_pad(stringr::str_trim(as.character(OCUPACION)),
                                    width = 3, pad = "0")
      )
    
    regis %>%
      dplyr::mutate(
        NIVEL_OCUPACION = dplyr::case_when(
          
          # 1990–1991
          ANIO %in% c(1990,1991) & OCUP_INT == 2 ~ "Unemployed",
          ANIO %in% c(1990,1991) & OCUP_INT %in% c(84) ~ "Armed forces",
          
          ANIO %in% c(1990,1991) & OCUP_INT %in% c(11,12,13,14,21,22,31,41,51) ~ "High",
          ANIO %in% c(1990,1991) & OCUP_INT %in% c(42,43,52,53,61,83) ~ "Medium",
          ANIO %in% c(1990,1991) & OCUP_INT %in% c(71,72,81,82) ~ "Low",
          ANIO %in% c(1990,1991) & OCUP_INT %in% c(98,99) ~ NA_character_,
          
          # 1992–2001 (2 dígitos)
          ANIO >= 1992 & ANIO <= 2001 & OCUP_INT %in% c(11,12,13,14,21) ~ "High",
          ANIO >= 1992 & ANIO <= 2001 & OCUP_INT %in% c(41,51,52,53,54,55,61,62,71) ~ "Medium",
          ANIO >= 1992 & ANIO <= 2001 & OCUP_INT %in% c(72,81,82) ~ "Low",
          ANIO >= 1992 & ANIO <= 2001 & OCUP_INT %in% c(2) ~ "Unemployed",
          ANIO >= 1992 & ANIO <= 2001 & OCUP_INT %in% c(83) ~ "Armed forces",
          ANIO >= 1992 & ANIO <= 2001 & OCUP_INT %in% c(98,99) ~ NA_character_,
          
          # 2002–2012 (2 dígitos)
          ANIO >= 2002 & ANIO <= 2012 & OCUP_INT %in% c(11,12,13,14,21) ~ "High",
          ANIO >= 2002 & ANIO <= 2012 & OCUP_INT %in% c(41,51,52,53,54,55,61,62,71) ~ "Medium",
          ANIO >= 2002 & ANIO <= 2012 & OCUP_INT %in% c(72,81,82) ~ "Low",
          ANIO >= 2002 & ANIO <= 2012 & OCUP_INT %in% c(2) ~ "Unemployed",
          ANIO >= 2002 & ANIO <= 2012 & OCUP_INT %in% c(83) ~ "Armed forces",
          ANIO >= 2002 & ANIO <= 2012 & OCUP_INT %in% c(97,98,99) ~ NA_character_,
          
          # 2013–2021 (1–11, 97–99)
          ANIO >= 2013 & ANIO <= 2021 & OCUP_INT %in% c(1,2) ~ "High",
          ANIO >= 2013 & ANIO <= 2021 & OCUP_INT %in% c(3,4,5,6,7,8) ~ "Medium",
          ANIO >= 2013 & ANIO <= 2021 & OCUP_INT %in% c(9) ~ "Low",
          ANIO >= 2013 & ANIO <= 2021 & OCUP_INT %in% c(10,11) ~ "Unemployed",
          ANIO >= 2013 & ANIO <= 2021 & OCUP_INT %in% c(97,98,99) ~ NA_character_,
          
          # 2022–2024 (3 dígitos; usa OCUP_3)
          ANIO >= 2022 & ANIO <= 2024 & OCUP_3 %in% c(
            "011","012","013","014","015","016","017","019",
            "021","022","023","024","025","026","027","028","029"
          ) ~ "High",
          
          ANIO >= 2022 & ANIO <= 2024 & OCUP_3 %in% c(
            "031","032","039","041","042","043","049",
            "051","052","053","061","062","063","069",
            "071","072","073","074","075","076","079",
            "081","082","083","089",
            "095"
          ) ~ "Medium",
          
          ANIO >= 2022 & ANIO <= 2024 & OCUP_3 %in% c(
            "091","092","093","094","096","097","098"
          ) ~ "Low",
          
          ANIO >= 2022 & ANIO <= 2024 & OCUP_3 %in% c("100","110") ~ "Unemployed",
          ANIO >= 2022 & ANIO <= 2024 & OCUP_3 %in% c("054") ~ "Armed forces",
          ANIO >= 2022 & ANIO <= 2024 & OCUP_3 %in% c("099","997","998","999") ~ NA_character_,
          
          TRUE ~ NA_character_
        )
      )
  })
  
  names(out) <- as.character(anios)
  out
}

work_classified_deaths_id <- recodificacion_ocupacion(classified_deaths_id, anios)
table(work_classified_deaths_id[["2012"]]$NIVEL_OCUPACION, useNA="ifany")


#----Asistencia-----
 recodificacion_asistencia <- function(registros_defun_clave_mun, anios){
  
  df_anio <- lapply(anios, function(anio){
    
    regis_anio <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(regis_anio)) return(NULL)
    
    regis_anio <- regis_anio %>% dplyr::mutate(ANIO = as.integer(anio))
    
    if (!("ASIST_MEDI" %in% names(regis_anio))) {
      regis_anio$ASISTENCIA <- NA_integer_
      return(regis_anio)
    }
    
    asist_int <- suppressWarnings(as.integer(as.character(regis_anio$ASIST_MEDI)))
    
    regis_anio %>%
      dplyr::mutate(
        ASIST_MEDI_int = asist_int,
        ASISTENCIA = dplyr::case_when(
          ASIST_MEDI_int == 1 ~ 1L,                 # sí
          ASIST_MEDI_int == 2 ~ 0L,                 # no
          ASIST_MEDI_int %in% c(9, 99) ~ NA_integer_,
          TRUE ~ NA_integer_
        )
      )
  })
  
  names(df_anio) <- as.character(anios)
  df_anio
}

assistance_classified_deaths_id <- recodificacion_asistencia(classified_deaths_id, anios)
table(assistance_classified_deaths_id[["2012"]]$ASISTENCIA, useNA="ifany")

#----Area urbana rural----
recodificacion_area <- function(registros_defun_clave_mun, anios){
  
  out <- lapply(anios, function(anio){
    
    regis <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(regis)) return(NULL)
    
    regis <- regis %>% dplyr::mutate(ANIO = as.integer(anio))
    
    # Antes de 2002: no existe la variable
    if (anio < 2002) {
      return(regis %>% dplyr::mutate(area = NA_character_))
    }
    
    # Si la variable no existe en el DF
    if (!("AREA_UR" %in% names(regis))) {
      return(regis %>% dplyr::mutate(area = NA_character_))
    }
    
    regis %>%
      dplyr::mutate(
        AREA_UR_num = suppressWarnings(as.integer(as.character(AREA_UR))),
        area = dplyr::case_when(
          AREA_UR_num == 1 ~ "Urbana",
          AREA_UR_num == 2 ~ "Rural",
          AREA_UR_num == 9 ~ NA_character_,
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::select(-AREA_UR_num)
  })
  
  names(out) <- as.character(anios)
  out
}

area_classified_deaths_id <- recodificacion_area(classified_deaths_id, anios)

table(area_classified_deaths_id[["2001"]]$area, useNA="ifany")  
table(area_classified_deaths_id[["2010"]]$area, useNA="ifany")


#----Contexto defuncion----
recodificacion_sitio_ocurr <- function(registros_defun_clave_mun, anios){
  
  out <- lapply(anios, function(anio){
    
    regis <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(regis)) return(NULL)
    
    regis <- regis %>% dplyr::mutate(ANIO = as.integer(anio))
    
    # detectar nombre variable
    var_sitio <- names(regis)[tolower(names(regis)) == "sitio_ocur"][1]
    
    if (is.na(var_sitio)) {
      regis$LUGAR_MUERTE <- NA_character_
      return(regis)
    }
    
    sitio_num <- suppressWarnings(as.integer(as.character(regis[[var_sitio]])))
    
    regis %>%
      dplyr::mutate(
        LUGAR_MUERTE = dplyr::case_when(
          # 1990–2003
          ANIO <= 2003 & sitio_num %in% c(1, 2) ~ "Hospital",
          ANIO <= 2003 & sitio_num %in% c(3, 4) ~ "Ambulatory",
          
          # 2004–2024
          ANIO >= 2004 & sitio_num %in% c(1,2,3,4,5,6,7,8,9) ~ "Hospital",
          ANIO >= 2004 & sitio_num %in% c(10,11,12) ~ "Ambulatory",
          
          # no especificado
          sitio_num %in% c(99) ~ NA_character_,
          
          TRUE ~ NA_character_
        )
      )
  })
  
  names(out) <- as.character(anios)
  out
}

contexto_classified_deaths_id <- recodificacion_sitio_ocurr(classified_deaths_id, anios)
table(contexto_classified_deaths_id[["2000"]]$LUGAR_MUERTE, useNA="ifany")
table(contexto_classified_deaths_id[["2010"]]$LUGAR_MUERTE, useNA="ifany")


#FUNCIÓN DE CLASIFICACION DE EDAD QUINQUENAL: Genera una columna adicional llamada "grupo_edad" dónde se clasifica su grupo quinqienal
## Requisitos: DF que contenga registros de defunciones con la columna "EDAD_AGRU" que codifica el grupo quinquenal según INEGI
clasificacion_edad_quinq <- function(registros_defun_clave_mun, anios) {
  niveles_edad <- c(
    "POB_00_04","POB_05_09","POB_10_14","POB_15_19","POB_20_24",
    "POB_25_29","POB_30_34","POB_35_39","POB_40_44","POB_45_49",
    "POB_50_54","POB_55_59","POB_60_64","POB_65_69","POB_70_74",
    "POB_75_79","POB_80_84","POB_85_mm")
  
  df_anio <- lapply(anios, function(anio){
    
    df <- registros_defun_clave_mun[[as.character(anio)]]
    
    df2 <- df %>%
      dplyr::mutate(
        ANIO = as.integer(anio),
        edad_agru = suppressWarnings(as.integer(as.character(EDAD_AGRU))),
        grupo_edad = dplyr::case_when(
          edad_agru %in% 1:5   ~ "POB_00_04",
          edad_agru == 6       ~ "POB_05_09",
          edad_agru == 7       ~ "POB_10_14",
          edad_agru == 8       ~ "POB_15_19",
          edad_agru == 9       ~ "POB_20_24",
          edad_agru == 10      ~ "POB_25_29",
          edad_agru == 11      ~ "POB_30_34",
          edad_agru == 12      ~ "POB_35_39",
          edad_agru == 13      ~ "POB_40_44",
          edad_agru == 14      ~ "POB_45_49",
          edad_agru == 15      ~ "POB_50_54",
          edad_agru == 16      ~ "POB_55_59",
          edad_agru == 17      ~ "POB_60_64",
          edad_agru == 18      ~ "POB_65_69",
          edad_agru == 19      ~ "POB_70_74",
          edad_agru == 20      ~ "POB_75_79",
          edad_agru == 21      ~ "POB_80_84",
          edad_agru %in% 22:29 ~ "POB_85_mm",
          TRUE                 ~ NA_character_), grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE))
    
    return(df2)
  })
  
  names(df_anio) <- as.character(anios)
  return(df_anio)
}

classified_deaths_id <- clasificacion_edad_quinq(classified_deaths_id, anios)

#FUNCIÓN PARA SEPARAR POR SEXO: Genera dos columnas adicionales "SEXO" (Hombre/Mujer) y "poblacion" para agregar la poblacion de cada sexo por municipio por año
recodificacion_sexo <- function(registros_defun_clave_mun, anios, poblacion_conapo_quinquenal){
  
  out <- lapply(anios, function(anio){
    
    # ---- Defunciones por municipio-año-sexo (totales y CV) ----
    reg <- registros_defun_clave_mun[[as.character(anio)]]
    if (is.null(reg)) return(NULL)
    
    sexo_int  <- suppressWarnings(as.integer(as.character(reg$SEXO)))
    muerte_cv <- suppressWarnings(as.integer(as.character(reg$muerte_cv)))
    
    def <- reg %>%
      dplyr::mutate(
        ANIO   = as.integer(anio),
        cvegeo = sprintf("%05d", as.integer(cvegeo)),
        sexo = dplyr::case_when(
          sexo_int == 1 ~ "HOMBRES",
          sexo_int == 2 ~ "MUJERES",
          TRUE ~ NA_character_
        ),
        muerte_cv = muerte_cv
      ) %>%
      dplyr::filter(!is.na(sexo), !is.na(cvegeo)) %>%
      dplyr::group_by(cvegeo, ANIO, sexo) %>%
      dplyr::summarise(
        defunciones_totales = dplyr::n(),
        defunciones_cv      = sum(muerte_cv == 1, na.rm = TRUE),
        .groups = "drop"
      )
    
    # ---- Población por municipio-año-sexo (CONAPO quinquenal) ----
    pob_year <- poblacion_conapo_quinquenal %>%
      dplyr::filter(AÑO == anio)
    
    pob_cols <- names(pob_year)[grepl("^POB_", names(pob_year))]
    if (length(pob_cols) == 0) stop("No encontré columnas POB_ en poblacion_conapo_quinquenal")
    
    pob <- pob_year %>%
      dplyr::transmute(
        cvegeo = sprintf("%05d", as.integer(CLAVE)),
        ANIO   = as.integer(AÑO),
        sexo   = dplyr::case_when(
          toupper(as.character(SEXO)) %in% c("HOMBRES","MUJERES") ~ toupper(as.character(SEXO)),
          toupper(as.character(SEXO)) %in% c("H","HOMBRE","MASCULINO") ~ "HOMBRES",
          toupper(as.character(SEXO)) %in% c("M","MUJER","FEMENINO")   ~ "MUJERES",
          TRUE ~ NA_character_
        ),
        poblacion = rowSums(dplyr::across(dplyr::all_of(pob_cols)), na.rm = TRUE),
        NOM_ENT, NOM_MUN
      ) %>%
      dplyr::filter(!is.na(sexo), !is.na(cvegeo))
    
    # ---- Join por municipio-año-sexo ----
    long <- pob %>%
      dplyr::left_join(def, by = c("cvegeo","ANIO","sexo")) %>%
      dplyr::mutate(
        defunciones_totales = dplyr::coalesce(defunciones_totales, 0L),
        defunciones_cv      = dplyr::coalesce(defunciones_cv, 0L)
      )
    
    # ---- Pasar a wide: hombres/mujeres en columnas ----
    wide <- long %>%
      tidyr::pivot_wider(
        names_from  = sexo,
        values_from = c(poblacion, defunciones_totales, defunciones_cv),
        values_fill = 0
      ) %>%
      dplyr::rename(
        poblacion_hombres             = poblacion_HOMBRES,
        poblacion_mujeres             = poblacion_MUJERES,
        defunciones_totales_hombres   = defunciones_totales_HOMBRES,
        defunciones_totales_mujeres   = defunciones_totales_MUJERES,
        defunciones_cv_hombres        = defunciones_cv_HOMBRES,
        defunciones_cv_mujeres        = defunciones_cv_MUJERES
      )
    
    wide
  })
  
  dplyr::bind_rows(out)
}

prueba_1 <- recodificacion_sexo(classified_deaths_id, anios, poblacion_conapo_quinquenal)

#Mortalidad DF
muertes_mun <- purrr::map_dfr(mortality_adjusted, ~ .x %>%
                                dplyr::transmute(
                                  cvegeo = sprintf("%05d", as.integer(cvegeo)),
                                  ANIO   = as.integer(ANIO),
                                  
                                  def_total = as.integer(def_mun_total),
                                  def_cvd   = as.integer(def_mun_cvd),
                                  def_noncvd = as.integer(def_mun_total - def_mun_cvd),
                                  
                                  pob_total = as.numeric(pob_mun_total),
                                  
                                  tasa_cruda_total = as.numeric(tasa_cruda_total),
                                  tasa_cruda_cvd   = as.numeric(tasa_cruda_cvd),
                                  
                                  tasa_adj_mun_total = as.numeric(tasa_ajustada_edad_mun_total),
                                  tasa_adj_mun_cvd   = as.numeric(tasa_ajustada_edad_mun_cvd)
                                ) %>%
                                dplyr::distinct(cvegeo, ANIO, .keep_all = TRUE)
)

mortalidad_sexo <- prueba_1 %>%
  dplyr::mutate(
    cvegeo = sprintf("%05d", as.integer(cvegeo)),
    ANIO   = as.integer(ANIO),
    
    poblacion_total_sexo = poblacion_hombres + poblacion_mujeres,
    
    defunciones_totales_sexo =
      defunciones_totales_hombres + defunciones_totales_mujeres,
    
    defunciones_cv_sexo =
      defunciones_cv_hombres + defunciones_cv_mujeres
  )

mortalidad_mun <- mortalidad_sexo %>%
  dplyr::left_join(muertes_mun, by = c("cvegeo","ANIO"))


#Agregar variables sociodemograficas a df
#Educacion
regis_edu_all <- recodificacion_educacion(classified_deaths_id, anios)

edu_mun_all <- purrr::map_dfr(regis_edu_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO)
    ) %>%
    dplyr::filter(!is.na(NIVEL_ESCOLARIDAD)) %>%
    dplyr::count(cvegeo, ANIO, NIVEL_ESCOLARIDAD, name = "n_all")
}) %>%
  tidyr::pivot_wider(
    names_from  = NIVEL_ESCOLARIDAD,
    values_from = n_all,
    values_fill = 0,
    names_glue  = "{NIVEL_ESCOLARIDAD}_all"
  )

edu_mun_cv <- purrr::map_dfr(regis_edu_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO),
      muerte_cvd = suppressWarnings(as.integer(as.character(muerte_cvd)))
    ) %>%
    dplyr::filter(muerte_cvd == 1, !is.na(NIVEL_ESCOLARIDAD)) %>%
    dplyr::count(cvegeo, ANIO, NIVEL_ESCOLARIDAD, name = "n_cv")
}) %>%
  tidyr::pivot_wider(
    names_from  = NIVEL_ESCOLARIDAD,
    values_from = n_cv,
    values_fill = 0,
    names_glue  = "{NIVEL_ESCOLARIDAD}_cv"
  )

mortalidad_mun <- mortalidad_mun %>%
  dplyr::left_join(edu_mun_all, by = c("cvegeo","ANIO")) %>%
  dplyr::left_join(edu_mun_cv,  by = c("cvegeo","ANIO")) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        low_education_all, medium_education_all, high_education_all,
        low_education_cv,  medium_education_cv,  high_education_cv
      ),
      ~ dplyr::coalesce(.x, 0L)
    )
  )

#Derechohabiencia
regis_dh_all <- recodificacion_derechohabiencia(classified_deaths_id, anios)

dh_mun_all <- purrr::map_dfr(regis_dh_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO)
    ) %>%
    dplyr::filter(!is.na(seguro)) %>%
    dplyr::count(cvegeo, ANIO, seguro, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = seguro,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::rename(
    no_derechohabiencia_all = `0`,
    derechohabiencia_all    = `1`
  )

dh_mun_cv <- purrr::map_dfr(regis_dh_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO),
      muerte_cvd = suppressWarnings(as.integer(as.character(muerte_cvd)))
    ) %>%
    dplyr::filter(muerte_cvd == 1, !is.na(seguro)) %>%
    dplyr::count(cvegeo, ANIO, seguro, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = seguro,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::rename(
    no_derechohabiencia_cv = `0`,
    derechohabiencia_cv    = `1`
  )

mortalidad_mun <- mortalidad_mun %>%
  dplyr::left_join(dh_mun_all, by = c("cvegeo","ANIO")) %>%
  dplyr::left_join(dh_mun_cv,  by = c("cvegeo","ANIO")) %>%
  dplyr::mutate(
    no_derechohabiencia_all = dplyr::coalesce(no_derechohabiencia_all, 0L),
    derechohabiencia_all    = dplyr::coalesce(derechohabiencia_all, 0L),
    no_derechohabiencia_cv  = dplyr::coalesce(no_derechohabiencia_cv, 0L),
    derechohabiencia_cv     = dplyr::coalesce(derechohabiencia_cv, 0L)
  )

#Lengua
regis_lengua_all <- recodificacion_indigena(classified_deaths_id, anios)

lengua_mun_all <- purrr::map_dfr(regis_lengua_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO)
    ) %>%
    dplyr::filter(!is.na(LENGUA_INDIGENA)) %>%
    dplyr::count(cvegeo, ANIO, LENGUA_INDIGENA, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = LENGUA_INDIGENA,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::rename(
    no_indigena_all = `0`,
    indigena_all    = `1`
  )

lengua_mun_cv <- purrr::map_dfr(regis_lengua_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO),
      muerte_cvd = suppressWarnings(as.integer(as.character(muerte_cvd)))
    ) %>%
    dplyr::filter(muerte_cvd == 1, !is.na(LENGUA_INDIGENA)) %>%
    dplyr::count(cvegeo, ANIO, LENGUA_INDIGENA, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = LENGUA_INDIGENA,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::rename(
    no_indigena_cv = `0`,
    indigena_cv    = `1`
  )


mortalidad_mun <- mortalidad_mun %>%
  dplyr::left_join(lengua_mun_all, by = c("cvegeo","ANIO")) %>%
  dplyr::left_join(lengua_mun_cv,  by = c("cvegeo","ANIO")) %>%
  dplyr::mutate(
    no_indigena_all = dplyr::coalesce(no_indigena_all, 0L),
    indigena_all    = dplyr::coalesce(indigena_all, 0L),
    no_indigena_cv  = dplyr::coalesce(no_indigena_cv, 0L),
    indigena_cv     = dplyr::coalesce(indigena_cv, 0L)
  )

#Ocupacion
regis_ocup_all <- recodificacion_ocupacion(classified_deaths_id, anios)

ocup_mun_all <- purrr::map_dfr(regis_ocup_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO)
    ) %>%
    dplyr::filter(!is.na(NIVEL_OCUPACION)) %>%
    dplyr::count(cvegeo, ANIO, NIVEL_OCUPACION, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = NIVEL_OCUPACION,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::rename(
    ocu_high_all       = High,
    ocu_medium_all     = Medium,
    ocu_low_all        = Low,
    ocu_unemployed_all = Unemployed,
    ocu_armed_all      = `Armed forces`
  )

ocup_mun_cv <- purrr::map_dfr(regis_ocup_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO),
      muerte_cvd = suppressWarnings(as.integer(as.character(muerte_cvd)))
    ) %>%
    dplyr::filter(muerte_cvd == 1, !is.na(NIVEL_OCUPACION)) %>%
    dplyr::count(cvegeo, ANIO, NIVEL_OCUPACION, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = NIVEL_OCUPACION,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::rename(
    ocu_high_cv       = High,
    ocu_medium_cv     = Medium,
    ocu_low_cv        = Low,
    ocu_unemployed_cv = Unemployed,
    ocu_armed_cv      = `Armed forces`
  )

mortalidad_mun <- mortalidad_mun %>%
  dplyr::left_join(ocup_mun_all, by = c("cvegeo","ANIO")) %>%
  dplyr::left_join(ocup_mun_cv,  by = c("cvegeo","ANIO")) %>%
  dplyr::mutate(
    ocu_high_all       = dplyr::coalesce(ocu_high_all, 0L),
    ocu_medium_all     = dplyr::coalesce(ocu_medium_all, 0L),
    ocu_low_all        = dplyr::coalesce(ocu_low_all, 0L),
    ocu_unemployed_all = dplyr::coalesce(ocu_unemployed_all, 0L),
    ocu_armed_all      = dplyr::coalesce(ocu_armed_all, 0L),
    
    ocu_high_cv        = dplyr::coalesce(ocu_high_cv, 0L),
    ocu_medium_cv      = dplyr::coalesce(ocu_medium_cv, 0L),
    ocu_low_cv         = dplyr::coalesce(ocu_low_cv, 0L),
    ocu_unemployed_cv  = dplyr::coalesce(ocu_unemployed_cv, 0L),
    ocu_armed_cv       = dplyr::coalesce(ocu_armed_cv, 0L)
  )

#Asistencia
regis_asist_all <- recodificacion_asistencia(classified_deaths_id, anios)

asist_mun_all <- purrr::map_dfr(regis_asist_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo     = sprintf("%05d", as.integer(cvegeo)),
      ANIO       = as.integer(ANIO),
      ASISTENCIA = suppressWarnings(as.integer(as.character(ASISTENCIA)))
    ) %>%
    dplyr::filter(!is.na(ASISTENCIA)) %>%
    dplyr::count(cvegeo, ANIO, ASISTENCIA, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = ASISTENCIA,
    values_from = n,
    values_fill = 0
  )

if (!("0" %in% names(asist_mun_all))) asist_mun_all$`0` <- 0L
if (!("1" %in% names(asist_mun_all))) asist_mun_all$`1` <- 0L

asist_mun_all <- asist_mun_all %>%
  dplyr::rename(
    asistencia_no_all = `0`,
    asistencia_si_all = `1`
  ) %>%
  dplyr::select(cvegeo, ANIO, asistencia_no_all, asistencia_si_all)

asist_mun_cv <- purrr::map_dfr(regis_asist_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo     = sprintf("%05d", as.integer(cvegeo)),
      ANIO       = as.integer(ANIO),
      muerte_cvd  = suppressWarnings(as.integer(as.character(muerte_cvd))),
      ASISTENCIA = suppressWarnings(as.integer(as.character(ASISTENCIA)))
    ) %>%
    dplyr::filter(muerte_cvd == 1, !is.na(ASISTENCIA)) %>%
    dplyr::count(cvegeo, ANIO, ASISTENCIA, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = ASISTENCIA,
    values_from = n,
    values_fill = 0
  )

if (!("0" %in% names(asist_mun_cv))) asist_mun_cv$`0` <- 0L
if (!("1" %in% names(asist_mun_cv))) asist_mun_cv$`1` <- 0L

asist_mun_cv <- asist_mun_cv %>%
  dplyr::rename(
    asistencia_no_cv = `0`,
    asistencia_si_cv = `1`
  ) %>%
  dplyr::select(cvegeo, ANIO, asistencia_no_cv, asistencia_si_cv)

mortalidad_mun <- mortalidad_mun %>%
  dplyr::mutate(
    cvegeo = sprintf("%05d", as.integer(cvegeo)),
    ANIO   = as.integer(ANIO)
  ) %>%
  # si ya existen columnas asistencia_*, esto evita que se conviertan en .x/.y
  dplyr::left_join(asist_mun_all, by = c("cvegeo","ANIO"), suffix = c("", "_tmp")) %>%
  dplyr::left_join(asist_mun_cv,  by = c("cvegeo","ANIO"), suffix = c("", "_tmp"))

# Consolidar (por si quedaron sufijos) y asegurar 0L
if (!("asistencia_no_all" %in% names(mortalidad_mun)) && ("asistencia_no_all_tmp" %in% names(mortalidad_mun))) {
  mortalidad_mun$asistencia_no_all <- mortalidad_mun$asistencia_no_all_tmp
}
if (!("asistencia_si_all" %in% names(mortalidad_mun)) && ("asistencia_si_all_tmp" %in% names(mortalidad_mun))) {
  mortalidad_mun$asistencia_si_all <- mortalidad_mun$asistencia_si_all_tmp
}
if (!("asistencia_no_cv" %in% names(mortalidad_mun)) && ("asistencia_no_cv_tmp" %in% names(mortalidad_mun))) {
  mortalidad_mun$asistencia_no_cv <- mortalidad_mun$asistencia_no_cv_tmp
}
if (!("asistencia_si_cv" %in% names(mortalidad_mun)) && ("asistencia_si_cv_tmp" %in% names(mortalidad_mun))) {
  mortalidad_mun$asistencia_si_cv <- mortalidad_mun$asistencia_si_cv_tmp
}

if (!("asistencia_no_all" %in% names(mortalidad_mun))) mortalidad_mun$asistencia_no_all <- 0L
if (!("asistencia_si_all" %in% names(mortalidad_mun))) mortalidad_mun$asistencia_si_all <- 0L
if (!("asistencia_no_cv"  %in% names(mortalidad_mun))) mortalidad_mun$asistencia_no_cv  <- 0L
if (!("asistencia_si_cv"  %in% names(mortalidad_mun))) mortalidad_mun$asistencia_si_cv  <- 0L

mortalidad_mun <- mortalidad_mun %>%
  dplyr::mutate(
    asistencia_no_all = dplyr::coalesce(asistencia_no_all, 0L),
    asistencia_si_all = dplyr::coalesce(asistencia_si_all, 0L),
    asistencia_no_cv  = dplyr::coalesce(asistencia_no_cv, 0L),
    asistencia_si_cv  = dplyr::coalesce(asistencia_si_cv, 0L)
  ) %>%
  dplyr::select(-dplyr::any_of(c("asistencia_no_all_tmp","asistencia_si_all_tmp","asistencia_no_cv_tmp","asistencia_si_cv_tmp")))

#Urbano rural
regis_area_all <- recodificacion_area(classified_deaths_id, anios)

area_mun_all <- purrr::map_dfr(regis_area_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO)
    ) %>%
    dplyr::filter(!is.na(area)) %>%
    dplyr::count(cvegeo, ANIO, area, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = area,
    values_from = n,
    values_fill = 0
  )

if (!("Urbana" %in% names(area_mun_all))) area_mun_all$Urbana <- 0L
if (!("Rural"  %in% names(area_mun_all))) area_mun_all$Rural  <- 0L

area_mun_all <- area_mun_all %>%
  dplyr::rename(
    urbana_all = Urbana,
    rural_all  = Rural
  ) %>%
  dplyr::select(cvegeo, ANIO, urbana_all, rural_all)

area_mun_cv <- purrr::map_dfr(regis_area_all, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO),
      muerte_cvd = suppressWarnings(as.integer(as.character(muerte_cvd)))
    ) %>%
    dplyr::filter(muerte_cvd == 1, !is.na(area)) %>%
    dplyr::count(cvegeo, ANIO, area, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = area,
    values_from = n,
    values_fill = 0
  )

if (!("Urbana" %in% names(area_mun_cv))) area_mun_cv$Urbana <- 0L
if (!("Rural"  %in% names(area_mun_cv))) area_mun_cv$Rural  <- 0L

area_mun_cv <- area_mun_cv %>%
  dplyr::rename(
    urbana_cv = Urbana,
    rural_cv  = Rural
  ) %>%
  dplyr::select(cvegeo, ANIO, urbana_cv, rural_cv)

mortalidad_mun <- mortalidad_mun %>%
  dplyr::left_join(area_mun_all, by = c("cvegeo","ANIO")) %>%
  dplyr::left_join(area_mun_cv,  by = c("cvegeo","ANIO")) %>%
  dplyr::mutate(
    urbana_all = dplyr::if_else(ANIO < 2002, NA_integer_, dplyr::coalesce(urbana_all, 0L)),
    rural_all  = dplyr::if_else(ANIO < 2002, NA_integer_, dplyr::coalesce(rural_all, 0L)),
    urbana_cv  = dplyr::if_else(ANIO < 2002, NA_integer_, dplyr::coalesce(urbana_cv, 0L)),
    rural_cv   = dplyr::if_else(ANIO < 2002, NA_integer_, dplyr::coalesce(rural_cv, 0L))
  )

#Contexto defuncion
lugar_mun_all <- purrr::map_dfr(contexto_classified_deaths_id, ~{
  .x %>%
    dplyr::mutate(
      cvegeo = sprintf("%05d", as.integer(cvegeo)),
      ANIO   = as.integer(ANIO)
    ) %>%
    dplyr::filter(!is.na(LUGAR_MUERTE)) %>%
    dplyr::count(cvegeo, ANIO, LUGAR_MUERTE, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = LUGAR_MUERTE,
    values_from = n,
    values_fill = 0
  )

if (!("Hospital"    %in% names(lugar_mun_all))) lugar_mun_all$Hospital    <- 0L
if (!("Ambulatory"  %in% names(lugar_mun_all))) lugar_mun_all$Ambulatory  <- 0L

lugar_mun_all <- lugar_mun_all %>%
  dplyr::rename(
    lugar_hospital_all   = Hospital,
    lugar_ambulatory_all = Ambulatory
  ) %>%
  dplyr::select(cvegeo, ANIO, lugar_hospital_all, lugar_ambulatory_all)

lugar_mun_cv <- purrr::map_dfr(contexto_classified_deaths_id, ~{
  .x %>%
    dplyr::mutate(
      cvegeo    = sprintf("%05d", as.integer(cvegeo)),
      ANIO      = as.integer(ANIO),
      muerte_cvd = suppressWarnings(as.integer(as.character(muerte_cvd)))
    ) %>%
    dplyr::filter(muerte_cvd == 1, !is.na(LUGAR_MUERTE)) %>%
    dplyr::count(cvegeo, ANIO, LUGAR_MUERTE, name = "n")
}) %>%
  tidyr::pivot_wider(
    names_from  = LUGAR_MUERTE,
    values_from = n,
    values_fill = 0
  )

if (!("Hospital"    %in% names(lugar_mun_cv))) lugar_mun_cv$Hospital    <- 0L
if (!("Ambulatory"  %in% names(lugar_mun_cv))) lugar_mun_cv$Ambulatory  <- 0L

lugar_mun_cv <- lugar_mun_cv %>%
  dplyr::rename(
    lugar_hospital_cv   = Hospital,
    lugar_ambulatory_cv = Ambulatory
  ) %>%
  dplyr::select(cvegeo, ANIO, lugar_hospital_cv, lugar_ambulatory_cv)

# 2) Join a mortalidad_mun
mortalidad_mun <- mortalidad_mun %>%
  dplyr::left_join(lugar_mun_all, by = c("cvegeo","ANIO")) %>%
  dplyr::left_join(lugar_mun_cv,  by = c("cvegeo","ANIO")) %>%
  dplyr::mutate(
    lugar_hospital_all   = dplyr::coalesce(lugar_hospital_all, 0L),
    lugar_ambulatory_all = dplyr::coalesce(lugar_ambulatory_all, 0L),
    lugar_hospital_cv    = dplyr::coalesce(lugar_hospital_cv, 0L),
    lugar_ambulatory_cv  = dplyr::coalesce(lugar_ambulatory_cv, 0L)
  )

#---- Cargar datos población INEGI ----
dir_inegi <- file.path(ruta_base, "Mario", "INEGI")

# -------- 1) Lector universal para tus CSV INEGI/INAFED --------
leer_csv_inegi <- function(path, encoding = "Latin1"){
  
  # detectar delimitador ( ; vs , )
  linea1 <- readLines(path, n = 1, warn = FALSE)
  delim <- ifelse(str_count(linea1, ";") > str_count(linea1, ","), ";", ",")
  
  df <- readr::read_delim(
    file = path,
    delim = delim,
    locale = locale(encoding = encoding),
    show_col_types = FALSE,
    guess_max = 200000
  )
  
  df %>%
    # quitar columnas totalmente vacías
    select(where(~ !all(is.na(.))))
}

# -------- 2) Limpieza estándar: grupo_edad + quitar basura final --------
limpiar_tabla_inegi <- function(df){
  
  df %>%
    rename(grupo_edad = 1) %>%
    mutate(grupo_edad = str_squish(as.character(grupo_edad))) %>%
    filter(!is.na(grupo_edad)) %>%
    filter(!str_detect(grupo_edad, regex("^fuente", ignore_case = TRUE)))
}

# -------- 3) Wrapper: leer + limpiar en un paso --------
leer_y_limpiar_csv_inegi <- function(path){
  df <- leer_csv_inegi(path)
  limpiar_tabla_inegi(df)
}

# -------- 4) Aplicarlo a TODOS tus archivos (recursivo) --------
files_csv <- list.files(dir_inegi,
                        pattern = "\\.csv$",
                        full.names = TRUE,
                        recursive = TRUE)

# Carga: lista nombrada por archivo (basename)
inegi_csv_limpios <- setNames(files_csv, basename(files_csv)) %>%
  lapply(leer_y_limpiar_csv_inegi)

#1990
INEGI_1990_edu <- inegi_csv_limpios[["INEGI_1990_edu.csv"]]
INEGI_1990_lengua <- inegi_csv_limpios[["INEGI_1990_lengua.csv"]]

# Lengua 1995 (inafed)
INEGI_1995_lengua <- inegi_csv_limpios[["inafed_lengua_1995.csv"]]

#2000
INEGI_2000_ss <- inegi_csv_limpios[["INEGI_2000_ss.csv"]]
INEGI_2000_lengua <- inegi_csv_limpios[["INEGI_2000_lengua.csv"]]
INEGI_2000_edu <- inegi_csv_limpios[["INEGI_2000_edu.csv"]]

#2005
INEGI_2005_ss <- inegi_csv_limpios[["INEGI_2005_ss.csv"]]
INEGI_2005_lengua <- inegi_csv_limpios[["INEGI_2005_lengua.csv"]]
INEGI_2005_edu <- inegi_csv_limpios[["INEGI_2005_edu.csv"]]

#2010
INEGI_2010_ss <- inegi_csv_limpios[["INEGI_2010_ss.csv"]]
INEGI_2010_lengua <- inegi_csv_limpios[["INEGI_2010_lengua.csv"]]
INEGI_2010_edu <- inegi_csv_limpios[["INEGI_2010_edu.csv"]]

#2020
INEGI_2020_ss <- inegi_csv_limpios[["INEGI_2020_ss.csv"]]
INEGI_2020_lengua <- inegi_csv_limpios[["INEGI_2020_lengua.csv"]]
INEGI_2020_edu <- inegi_csv_limpios[["INEGI_2020_edu.csv"]]

#Cargar bases de empleo
library(readxl)

dir_empleo <- file.path(ruta_base, "Mario", "INEGI", "Empleo")

# renombrar todos los .xls -> .xlsx
files_xls  <- list.files(dir_empleo, pattern="\\.xls$", full.names=TRUE)
files_xlsx <- sub("\\.xls$", ".xlsx", files_xls)

file.rename(files_xls, files_xlsx)

# cargar todos los .xlsx
files_xlsx <- list.files(dir_empleo, pattern="\\.xlsx$", full.names=TRUE)

empleo_ok <- setNames(files_xlsx, basename(files_xlsx)) |>
  lapply(readxl::read_excel)

# objetos individuales (si quieres)
ENE_empleo     <- empleo_ok[["ENE_empleo.xlsx"]]
ENE_desempleo  <- empleo_ok[["ENE_desempleo.xlsx"]]
ENOE_empleo    <- empleo_ok[["ENOE_empleo.xlsx"]]
ENOE_desempleo <- empleo_ok[["ENOE_desempleo.xlsx"]]

#Limpiar
path <- file.path(dir_empleo, "ENE_empleo.xlsx")

ENE_empleo <- read_excel(path, skip = 4)
ENE_empleo <- ENE_empleo %>%
  rename(periodo = 1) %>%
  filter(!is.na(periodo)) %>%
  filter(!str_detect(periodo, regex("^fuente", ignore_case = TRUE)))

path <- file.path(dir_empleo, "ENE_desempleo.xlsx")
ENE_desempleo <- read_excel(path, skip = 4)
ENE_desempleo <- ENE_desempleo %>%
  rename(periodo = 1) %>%
  mutate(periodo = str_squish(as.character(periodo))) %>%
  filter(!is.na(periodo)) %>%
  filter(!str_detect(periodo, regex("^fuente", ignore_case = TRUE))) %>%
  filter(!str_detect(periodo, regex("^notas",  ignore_case = TRUE))) %>%
  filter(!str_detect(periodo, regex("^El 27",  ignore_case = TRUE))) %>%
  filter(!str_detect(periodo, regex("^\\-",    ignore_case = TRUE)))  # líneas que empiezan con "-"

path <- file.path(dir_empleo, "ENOE_empleo.xlsx")
ENOE_empleo <- read_excel(path, skip = 5)
ENOE_empleo <- ENOE_empleo %>%
  rename(periodo = 1) %>%
  filter(!is.na(periodo)) %>%
  filter(!str_detect(periodo, regex("^fuente", ignore_case = TRUE)))

path <- file.path(dir_empleo, "ENOE_desempleo.xlsx")
ENOE_desempleo <- read_excel(path, skip = 4)
ENOE_desempleo <- ENOE_desempleo %>%
  rename(periodo = 1) %>%
  filter(!is.na(periodo)) %>%
  filter(!str_detect(periodo, regex("^fuente", ignore_case = TRUE))) %>%
  filter(!str_detect(periodo, regex("^notas",  ignore_case = TRUE))) %>%
  filter(!str_detect(periodo, regex("^El 27",  ignore_case = TRUE))) %>%
  filter(!str_detect(periodo, regex("^\\-",    ignore_case = TRUE)))  # líneas que empiezan con "-"
  


#---- Cargar datos población INEGI 2 ----

cargar_censo_inegi <- function(ruta_base, anio, encoding = "Latin1") {
  
  # Carpeta raíz de censos
  carpeta_censos <- list.files(
    ruta_base,
    pattern = "^Censos_conteos_vivienda_INEGI$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(carpeta_censos) == 0) stop("No se encontró 'Censos_conteos_vivienda_INEGI' en ruta_base.")
  carpeta_censos <- carpeta_censos[1]
  
  # Ruta del año -> conjunto_de_datos
  ruta_anio <- file.path(carpeta_censos, as.character(anio), "conjunto_de_datos")
  if (!dir.exists(ruta_anio)) stop("No existe la ruta: ", ruta_anio)
  
  # Patrón por año (usa tus nombres)
  patrones <- list(
    "2000" = "^cgpv2000_iter_00.*\\.csv$",
    "2005" = "^cpv2005_iter_00.*\\.csv$",
    "2010" = "^iter_00_cpv2010.*\\.csv$",
    "2020" = "^conjunto_de_datos_iter_00CSV20.*\\.csv$"
  )
  
  patron <- patrones[[as.character(anio)]]
  if (is.null(patron)) stop("Año no soportado: ", anio, " (solo 2000, 2005, 2010, 2020)")
  
  # Buscar archivo
  archivos <- list.files(ruta_anio, pattern = patron, full.names = TRUE, ignore.case = TRUE)
  if (length(archivos) == 0) stop("No se encontró CSV para ", anio, " con patrón: ", patron)
  if (length(archivos) > 1) {
    message("Más de un match para ", anio, ". Tomando el primero:\n", paste(archivos, collapse = "\n"))
  }
  
  archivo <- archivos[1]
  
  # Leer CSV
  df <- readr::read_csv(
    archivo,
    show_col_types = FALSE,
    locale = readr::locale(encoding = encoding)
  ) %>%
    dplyr::rename_with(tolower)
  
  
  df
}

cargar_censos_inegi <- function(ruta_base, anios = c(2000, 2005, 2010, 2020), encoding = "Latin1") {
  out <- lapply(anios, function(a) cargar_censo_inegi(ruta_base, a, encoding = encoding))
  names(out) <- as.character(anios)
  out
}

# Cargarlos
censos <- cargar_censos_inegi(ruta_base)

# Crear clave municipal 
censos_municipales <- lapply(censos, function(df) {
  df %>%
    filter(loc == "0000") %>%
    mutate(
      entidad = str_pad(entidad, 2, pad = "0"),
      mun     = str_pad(mun, 3, pad = "0"),
      cvegeo  = paste0(entidad, mun)
    )
})


# ---- Cargar poblaciones INEGI
to_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("N/D", "ND", "NA", "N.A.", "S/D", "", ".", "NULL")] <- NA_character_
  readr::parse_number(x)
}

pick_first_existing <- function(df, candidates) {
  nm <- names(df)
  hit <- candidates[candidates %in% nm]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

fix_cvegeo5 <- function(x) str_pad(as.character(x), 5, pad = "0")

filter_municipios <- function(df) {
  df %>%
    rename_with(tolower) %>%
    mutate(cvegeo = fix_cvegeo5(cvegeo)) %>%
    filter(loc == "0000") %>%
    filter(str_detect(cvegeo, "^\\d{5}$")) %>%
    filter(cvegeo != "00000") %>%
    filter(substr(cvegeo, 3, 5) != "000")
}


# 1) Lengua indígena (HLI): pob_hli + den_hli + prop_hli
armonizar_pob_indigena <- function(df, anio) {
  
  df <- df %>% dplyr::rename_with(tolower)
  
  # ---- 2000: no tienes denominador 5+, usar total municipal (fallback) ----
  if (anio == 2000) {
    if (!"p5_hli" %in% names(df)) stop("No existe p5_hli en 2000")
    
    df <- df %>% dplyr::mutate(pob_hli = to_num(p5_hli))
    
    den_var_total <- pick_first_existing(df, c(
      "pobtot","p_total","ptotal","pob_total","poblacion","ptot"
    ))
    if (is.na(den_var_total)) {
      stop("2000: no encuentro población total municipal (pobtot/p_total/etc.).")
    }
    
    df <- df %>%
      dplyr::mutate(
        den_hli  = to_num(.data[[den_var_total]]),
        prop_hli = dplyr::if_else(den_hli > 0, pob_hli / den_hli, NA_real_)
      )
    
    return(df)
  }
  
  # ---- 2005: 5+ ----
  if (anio == 2005) {
    if (!"p5ymahli" %in% names(df)) stop("No existe p5ymahli en 2005")
    df <- df %>% dplyr::mutate(pob_hli = to_num(p5ymahli))
    
    den_var <- pick_first_existing(df, c("p5ym_tot","p5_tot","p5_total","p5ymas","p_5ymas","p5ym"))
    if (is.na(den_var)) stop("2005: no encontré denominador elegible 5+ (p5ym_tot o similar).")
    
    df <- df %>%
      dplyr::mutate(
        den_hli  = to_num(.data[[den_var]]),
        prop_hli = dplyr::if_else(den_hli > 0, pob_hli / den_hli, NA_real_)
      )
    return(df)
  }
  
  # ---- 2010: 3+ ----
  if (anio == 2010) {
    if (!"p3ym_hli" %in% names(df)) stop("No existe p3ym_hli en 2010")
    df <- df %>% dplyr::mutate(pob_hli = to_num(p3ym_hli))
    
    den_var <- pick_first_existing(df, c("p3ym_tot","p3_tot","p3_total","p_3ymas","p3ymas","p3ym"))
    if (is.na(den_var)) stop("2010: no encontré denominador elegible 3+ (p3ym_tot o similar).")
    
    df <- df %>%
      dplyr::mutate(
        den_hli  = to_num(.data[[den_var]]),
        prop_hli = dplyr::if_else(den_hli > 0, pob_hli / den_hli, NA_real_)
      )
    return(df)
  }
  
  # ---- 2020: 5+ ----
  if (anio == 2020) {
    if (!"p5_hli" %in% names(df)) stop("No existe p5_hli en 2020")
    df <- df %>% dplyr::mutate(pob_hli = to_num(p5_hli))
    
    den_var <- pick_first_existing(df, c("p5_tot","p5_total","p5ymas","p_5ymas","p5ym_tot","p5ym"))
    if (is.na(den_var)) stop("2020: no encontré denominador elegible 5+ (p5_tot o similar).")
    
    df <- df %>%
      dplyr::mutate(
        den_hli  = to_num(.data[[den_var]]),
        prop_hli = dplyr::if_else(den_hli > 0, pob_hli / den_hli, NA_real_)
      )
    return(df)
  }
  
  stop("Año no soportado: ", anio)
}

# ============================================================
# 2) Seguridad social: ss_si + ss_no + den_ss + prop_ss_no
# ============================================================
seguridad_social_censo <- function(df, anio) {
  
  df <- df %>% rename_with(tolower)
  
  if (anio == 2000) {
    req <- c("pderimss","pderiste","psderss")
    miss <- setdiff(req, names(df)); if (length(miss)) stop("Faltan vars 2000: ", paste(miss, collapse=", "))
    df <- df %>%
      mutate(
        pderimss = to_num(pderimss),
        pderiste = to_num(pderiste),
        psderss  = to_num(psderss),
        ss_si    = pderimss + pderiste,
        ss_no    = psderss
      )
    return(df)
  }
  
  if (anio == 2005) {
    req <- c("p_imss","p_issste","p_segpop","p_sinder")
    miss <- setdiff(req, names(df)); if (length(miss)) stop("Faltan vars 2005: ", paste(miss, collapse=", "))
    df <- df %>%
      mutate(
        p_imss    = to_num(p_imss),
        p_issste  = to_num(p_issste),
        p_segpop  = to_num(p_segpop),
        p_sinder  = to_num(p_sinder),
        ss_si     = p_imss + p_issste + p_segpop,
        ss_no     = p_sinder
      )
    return(df)
  }
  
  if (anio == 2010) {
    req <- c("pder_imss","pder_iste","pder_istee","pder_segp","psinder")
    miss <- setdiff(req, names(df)); if (length(miss)) stop("Faltan vars 2010: ", paste(miss, collapse=", "))
    df <- df %>%
      mutate(
        pder_imss  = to_num(pder_imss),
        pder_iste  = to_num(pder_iste),
        pder_istee = to_num(pder_istee),
        pder_segp  = to_num(pder_segp),
        psinder    = to_num(psinder),
        ss_si      = pder_imss + pder_iste + pder_istee + pder_segp,
        ss_no      = psinder
      )
    return(df)
  }
  
  if (anio == 2020) {
    req <- c("pder_ss","psinder","pafil_ipriv")
    miss <- setdiff(req, names(df)); if (length(miss)) stop("Faltan vars 2020: ", paste(miss, collapse=", "))
    df <- df %>%
      mutate(
        pder_ss     = to_num(pder_ss),
        psinder     = to_num(psinder),
        pafil_ipriv = to_num(pafil_ipriv),
        # Público = total con SS - afiliación privada (no permitir negativos)
        ss_si       = pmax(pder_ss - pafil_ipriv, 0),
        # No = sin SS + privada
        ss_no       = psinder + pafil_ipriv
      )
    return(df)
  }
  
  stop("Año no soportado: ", anio)
}

# ============================================================
# 3) Educación: edu_low/medium/high + den_edu + prop_edu_*
#    (denominador interno 15+)
# ============================================================
educacion_censo <- function(df, anio) {
  
  df <- df %>% rename_with(tolower)
  
  vars <- if (anio == 2000) {
    list(
      low    = c("p15_sinstr","p15_sprima","p15_cprima"),
      medium = c("p15_ssecu","p15_csecu","p15_sinsec","p15_consec"),
      high   = c("p15_cmedss")
    )
  } else if (anio == 2005) {
    list(
      low    = c("p15ymase","p15ym_ebin"),
      medium = c("p15ym_ebc"),
      high   = c("p15ymapb")
    )
  } else if (anio %in% c(2010, 2020)) {
    list(
      low    = c("p15ym_se","p15pri_in","p15pri_co"),
      medium = c("p15sec_in","p15sec_co"),
      high   = c("p18ym_pb")
    )
  } else {
    stop("Año no soportado: ", anio)
  }
  
  req  <- unlist(vars, use.names = FALSE)
  miss <- setdiff(req, names(df)); if (length(miss)) stop("Faltan vars edu ", anio, ": ", paste(miss, collapse=", "))
  
  df %>%
    mutate(
      edu_low    = rowSums(across(all_of(vars$low),    to_num), na.rm = TRUE),
      edu_medium = rowSums(across(all_of(vars$medium), to_num), na.rm = TRUE),
      edu_high   = rowSums(across(all_of(vars$high),   to_num), na.rm = TRUE)
    )
}

# ============================================================
# Pipeline completo
# censos_municipales: lista nombrada "2000","2005","2010","2020"
# ============================================================
stopifnot(is.list(censos_municipales))
stopifnot(all(grepl("^\\d{4}$", names(censos_municipales))))

censos_municipales_proc <- purrr::imap(censos_municipales, ~{
  anio <- as.integer(.y)
  .x %>%
    rename_with(tolower) %>%
    armonizar_pob_indigena(anio) %>%
    seguridad_social_censo(anio) %>%
    educacion_censo(anio)
})

# quedarte con municipios
censos_mun <- purrr::imap(censos_municipales_proc, ~ filter_municipios(.x))

# ============================================================
# Tabla final de denominadores y proporciones (por municipio y año censal)
# ============================================================
denoms_censo <- purrr::imap_dfr(censos_mun, ~{
  anio_censo <- as.integer(.y)
  .x %>%
    transmute(
      cvegeo = fix_cvegeo5(cvegeo),
      anio_censo = anio_censo,
      
      # Educación (15+)
      den_edu = edu_low + edu_medium + edu_high,
      prop_edu_low    = if_else(den_edu > 0, edu_low / den_edu, NA_real_),
      prop_edu_medium = if_else(den_edu > 0, edu_medium / den_edu, NA_real_),
      prop_edu_high   = if_else(den_edu > 0, edu_high / den_edu, NA_real_),
      
      # Seguridad social (interno)
      den_ss = ss_si + ss_no,
      prop_ss_no = if_else(den_ss > 0, ss_no / den_ss, NA_real_),
      prop_ss_si = if_else(den_ss > 0, ss_si / den_ss, NA_real_),
      
      # Lengua indígena (elegible 3+/5+ según año, si se encontró denominador)
      den_hli = den_hli,
      prop_hli = prop_hli
    )
})

# ============================================================
# Chequeos rápidos
# ============================================================

# Educación: sumas ~1 (cuando den_edu>0)
denoms_censo %>%
  mutate(prop_sum_edu = prop_edu_low + prop_edu_medium + prop_edu_high) %>%
  summarise(
    edu_p01 = quantile(prop_sum_edu, 0.01, na.rm = TRUE),
    edu_med = median(prop_sum_edu, na.rm = TRUE),
    edu_p99 = quantile(prop_sum_edu, 0.99, na.rm = TRUE)
  )

# SS: sumas ~1
denoms_censo %>%
  mutate(prop_sum_ss = prop_ss_no + prop_ss_si) %>%
  summarise(
    ss_p01 = quantile(prop_sum_ss, 0.01, na.rm = TRUE),
    ss_med = median(prop_sum_ss, na.rm = TRUE),
    ss_p99 = quantile(prop_sum_ss, 0.99, na.rm = TRUE)
  )

# HLI: cuántos NA por falta de denominador elegible
denoms_censo %>%
  summarise(
    n_total = n(),
    n_hli_na = sum(is.na(prop_hli)),
    pct_hli_na = mean(is.na(prop_hli)) * 100
  )

# ============================================================
# (Opcional) Integración a mortalidad municipal-año usando anio_censo/anio_rezago
# Supone que en tu mortalidad tienes:
#   mortalidad_mun: cvegeo, ANIO y una columna anio_censo (o anio_rezago) ya asignada
# ============================================================

# Ejemplo si tu columna se llama anio_rezago (mapeo año->censo más cercano o previo)
# mortalidad_mun <- mortalidad_mun %>%
#   left_join(denoms_censo, by = c("cvegeo" = "cvegeo", "anio_rezago" = "anio_censo"))

diag_hli <- purrr::imap_dfr(censos_mun, ~{
  anio <- as.integer(.y)
  tibble(
    anio_censo = anio,
    n = nrow(.x),
    n_na_prop_hli = sum(is.na(.x$prop_hli)),
    n_na_den_hli  = sum(is.na(.x$den_hli)),
    den_hli_min = suppressWarnings(min(.x$den_hli, na.rm = TRUE)),
    den_hli_max = suppressWarnings(max(.x$den_hli, na.rm = TRUE))
  )
})

diag_hli

#---- Calculadr denominador CONAPO/INEGI
head(poblacion_conapo_quinquenal)

pob_long <- poblacion_conapo_quinquenal %>%
  mutate(
    cvegeo = str_pad(as.character(CLAVE), 5, pad = "0"),
    ANIO   = AÑO
  ) %>%
  pivot_longer(
    cols = starts_with("POB_"),
    names_to = "grupo_edad",
    values_to = "poblacion"
  ) %>%
  select(cvegeo, SEXO, ANIO, grupo_edad, poblacion)

pob_long <- pob_long %>%
  group_by(cvegeo, ANIO, grupo_edad) %>%
  summarise(poblacion = sum(poblacion), .groups = "drop")

pob_long <- pob_long %>%
  mutate(
    anio_censo = case_when(
      ANIO >= 1990 & ANIO < 2005 ~ 2000,
      ANIO >= 2005 & ANIO < 2010 ~ 2005,
      ANIO >= 2010 & ANIO < 2020 ~ 2010,
      ANIO >= 2020               ~ 2020
    )
  )


pob_long <- pob_long %>%
  left_join(denoms_censo,
            by = c("cvegeo", "anio_censo"))

pob_long <- pob_long %>%
  mutate(
    pob_edu_low    = round(poblacion * prop_edu_low),
    pob_edu_medium = round(poblacion * prop_edu_medium),
    pob_edu_high   = poblacion - pob_edu_low - pob_edu_medium,
    pob_ss_si      = round (poblacion * prop_ss_si),
    pob_ss_no      = round (poblacion * prop_ss_no),
    pob_indi       = round (poblacion * prop_hli)
  )

denominadores_tasas <- pob_long %>%
  dplyr::group_by(ANIO, grupo_edad) %>%
  dplyr::summarise(
    poblacion = sum(poblacion, na.rm = TRUE),
    pob_edu_low    = sum(pob_edu_low, na.rm = TRUE),
    pob_edu_medium = sum(pob_edu_medium, na.rm = TRUE),
    pob_edu_high   = sum(pob_edu_high, na.rm = TRUE),
    pob_ss_si   = sum(pob_ss_si, na.rm = TRUE),
    pob_ss_no   = sum(pob_ss_no, na.rm = TRUE),
    pob_indi   = sum(pob_indi, na.rm = TRUE),
    .groups = "drop"
  )


# Crear tasa estratificada por educacion
# diccionario edad
deaths_nat_age <- defunciones_grupo_edad %>%
  dplyr::mutate(ANIO = as.integer(ANIO)) %>%
  dplyr::rename(DEATHS = def_cvd) %>%
  dplyr::select(ANIO, grupo_edad, DEATHS)

mort_nat_age <- deaths_nat_age %>%
  dplyr::left_join(
    denominadores_tasas %>% dplyr::select(ANIO, grupo_edad, poblacion),
    by = c("ANIO","grupo_edad")
  ) %>%
  dplyr::mutate(rate_age = DEATHS / poblacion)

w_std <- denominadores_tasas %>%
  dplyr::filter(ANIO == 2010) %>%
  dplyr::transmute(
    grupo_edad,
    w = poblacion / sum(poblacion)
  )

tasas_ajustadas_2 <- mort_nat_age %>%
  dplyr::left_join(w_std, by = "grupo_edad") %>%
  dplyr::group_by(ANIO) %>%
  dplyr::summarise(
    tasa_adj = sum(rate_age * w, na.rm = TRUE),
    tasa_adj_100k = tasa_adj * 100000,
    .groups = "drop"
  )

map_edad <- dplyr::tibble(
  EDAD_AGRU  = sprintf("%02d", 6:23),
  grupo_edad = c(
    "POB_00_04","POB_05_09","POB_10_14","POB_15_19","POB_20_24","POB_25_29",
    "POB_30_34","POB_35_39","POB_40_44","POB_45_49","POB_50_54","POB_55_59",
    "POB_60_64","POB_65_69","POB_70_74","POB_75_79","POB_80_84","POB_85_mm"
  )
)

def_cat_age <- defunciones_edad_categorias %>%
  dplyr::mutate(
    ANIO = as.integer(ANIO),
    EDAD_AGRU = sprintf("%02d", as.integer(EDAD_AGRU))
  ) %>%
  dplyr::left_join(map_edad, by = "EDAD_AGRU") %>%
  dplyr::filter(!is.na(grupo_edad)) %>%
  dplyr::select(ANIO, grupo_edad, categoria, nivel, DEATHS)

#Educacion
den_edu_long <- denominadores_tasas %>%
  dplyr::select(ANIO, grupo_edad, pob_edu_low, pob_edu_medium, pob_edu_high) %>%
  tidyr::pivot_longer(starts_with("pob_edu_"),
                      names_to = "var", values_to = "POP") %>%
  dplyr::mutate(
    categoria = "educacion",
    nivel = dplyr::case_when(
      var == "pob_edu_low"    ~ "low_education",
      var == "pob_edu_medium" ~ "medium_education",
      var == "pob_edu_high"   ~ "high_education"
    )
  ) %>%
  dplyr::select(ANIO, grupo_edad, categoria, nivel, POP)

tasa_adj_edu <- def_cat_age %>%
  dplyr::filter(categoria == "educacion") %>%
  dplyr::left_join(den_edu_long, by = c("ANIO","grupo_edad","categoria","nivel")) %>%
  dplyr::mutate(rate_age = DEATHS / POP) %>%
  dplyr::left_join(w_std, by = "grupo_edad") %>%
  dplyr::group_by(ANIO, nivel) %>%
  dplyr::summarise(
    tasa_adj_100k = sum(rate_age * w, na.rm = TRUE) * 100000,
    .groups = "drop"
  )

#Lengua indigena
den_len_long <- denominadores_tasas %>%
  dplyr::select(ANIO, grupo_edad, poblacion, pob_indi) %>%
  dplyr::mutate(pob_no_indi = pmax(poblacion - pob_indi, 0)) %>%
  tidyr::pivot_longer(c(pob_indi, pob_no_indi),
                      names_to = "var", values_to = "POP") %>%
  dplyr::mutate(
    categoria = "lengua",
    nivel = dplyr::if_else(var == "pob_indi", "indigena", "no_indigena")
  ) %>%
  dplyr::select(ANIO, grupo_edad, categoria, nivel, POP)

tasa_adj_lengua <- def_cat_age %>%
  dplyr::filter(categoria == "lengua") %>%
  dplyr::left_join(den_len_long, by = c("ANIO","grupo_edad","categoria","nivel")) %>%
  dplyr::mutate(rate_age = DEATHS / POP) %>%
  dplyr::left_join(w_std, by = "grupo_edad") %>%
  dplyr::group_by(ANIO, nivel) %>%
  dplyr::summarise(
    tasa_adj_100k = sum(rate_age * w, na.rm = TRUE) * 100000,
    .groups = "drop"
  )

#Seguridad social
den_ss_long <- denominadores_tasas %>%
  dplyr::select(ANIO, grupo_edad, pob_ss_no, pob_ss_si) %>%
  tidyr::pivot_longer(c(pob_ss_no, pob_ss_si),
                      names_to = "var", values_to = "POP") %>%
  dplyr::mutate(
    categoria = "derechohabiencia",
    nivel = dplyr::if_else(var == "pob_ss_si", "derechohabiencia", "no_derechohabiencia")
  ) %>%
  dplyr::select(ANIO, grupo_edad, categoria, nivel, POP)

tasa_adj_ss <- def_cat_age %>%
  dplyr::filter(categoria == "derechohabiencia") %>%
  dplyr::left_join(den_ss_long, by = c("ANIO","grupo_edad","categoria","nivel")) %>%
  dplyr::mutate(rate_age = DEATHS / POP) %>%
  dplyr::left_join(w_std, by = "grupo_edad") %>%
  dplyr::group_by(ANIO, nivel) %>%
  dplyr::summarise(
    tasa_adj_100k = sum(rate_age * w, na.rm = TRUE) * 100000,
    .groups = "drop"
  )

summary(tasa_adj_edu$tasa_adj_100k)
summary(tasa_adj_lengua$tasa_adj_100k)
summary(tasa_adj_ss$tasa_adj_100k)


plot_df <- bind_rows(
  tasa_adj_edu   %>% mutate(variable = "Educación"),
  tasa_adj_lengua %>% mutate(variable = "Lengua indígena"),
  tasa_adj_ss    %>% mutate(variable = "Seguridad social")
) %>%
  mutate(
    variable = factor(variable, levels = c("Educación","Lengua indígena","Seguridad social")),
    # etiquetas más legibles
    nivel_lbl = case_when(
      variable == "Educación" & nivel == "low_education"    ~ "Baja",
      variable == "Educación" & nivel == "medium_education" ~ "Media",
      variable == "Educación" & nivel == "high_education"   ~ "Alta",
      variable == "Lengua indígena" & nivel == "indigena"      ~ "Indígena",
      variable == "Lengua indígena" & nivel == "no_indigena"   ~ "No indígena",
      variable == "Seguridad social" & nivel == "no_derechohabiencia" ~ "Sin SS",
      variable == "Seguridad social" & nivel == "derechohabiencia"    ~ "Con SS",
      TRUE ~ as.character(nivel)
    )
  )

ggplot(filter(plot_df, variable == "Educación"),
       aes(ANIO, tasa_adj_100k, group = nivel_lbl, linetype = nivel_lbl)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Año", y = "Tasa ajustada por edad (por 100,000)", linetype = "") +
  theme_minimal(base_size = 12)


ggplot(filter(plot_df, variable == "Lengua indígena"),
       aes(ANIO, tasa_adj_100k, group = nivel_lbl, linetype = nivel_lbl)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Año", y = "Tasa ajustada por edad (por 100,000)", linetype = "") +
  theme_minimal(base_size = 12)

ggplot(filter(plot_df, variable == "Seguridad social"),
       aes(ANIO, tasa_adj_100k, group = nivel_lbl, linetype = nivel_lbl)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Año", y = "Tasa ajustada por edad (por 100,000)", linetype = "") +
  theme_minimal(base_size = 12)

#---- Calcular tasas estratificadas ----
registros_edad_quinq <- clasificacion_edad_quinq(classified_deaths_id, anios)

def_year <- purrr::imap_dfr(
  registros_edad_quinq,
  ~{
    anio <- as.integer(.y)
    
    .x %>%
      dplyr::mutate(
        cvegeo    = sprintf("%05d", suppressWarnings(as.integer(as.character(cvegeo)))),
        ENT_RESID = substr(cvegeo, 1, 2),
        ANIO      = anio,
        SEXO      = dplyr::recode(as.character(SEXO),
                                  `1` = "Male",
                                  `2` = "Female",
                                  .default = NA_character_)
      ) %>%
      dplyr::filter(!is.na(grupo_edad), !is.na(SEXO)) %>%
      dplyr::group_by(cvegeo, ENT_RESID, ANIO, SEXO, grupo_edad) %>%
      dplyr::summarise(
        def_total = dplyr::n(),
        def_cvd   = sum(muerte_cvd, na.rm = TRUE),
        .groups = "drop"
      )
  }
)

#Por sexo
niveles_edad <- c(
  "POB_00_04","POB_05_09","POB_10_14","POB_15_19","POB_20_24",
  "POB_25_29","POB_30_34","POB_35_39","POB_40_44","POB_45_49",
  "POB_50_54","POB_55_59","POB_60_64","POB_65_69","POB_70_74",
  "POB_75_79","POB_80_84","POB_85_mm"
)

pob_long_sexo <- poblacion_conapo_quinquenal %>%
  dplyr::mutate(
    cvegeo = sprintf("%05d", as.integer(CLAVE)),
    ENT_RESID = sprintf("%02d", as.integer(CLAVE_ENT)),
    ANIO = as.integer(AÑO),
    SEXO = ifelse(SEXO == "HOMBRES", "Male", "Female")
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(niveles_edad),
    names_to = "grupo_edad",
    values_to = "poblacion"
  ) %>%
  dplyr::mutate(
    grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE)
  ) %>%
  dplyr::group_by(cvegeo, ENT_RESID, ANIO, SEXO, grupo_edad) %>%
  dplyr::summarise(poblacion = sum(poblacion), .groups = "drop")

df_age_mun_sexo <- pob_long_sexo %>%
  dplyr::left_join(def_year,
                   by = c("cvegeo","ENT_RESID","ANIO","SEXO","grupo_edad")) %>%
  dplyr::mutate(
    def_cvd = dplyr::coalesce(def_cvd, 0L)
  ) %>%
  dplyr::filter(poblacion > 0)

pesos_2010 <- pob_long_sexo %>%
  dplyr::filter(ANIO == 2010) %>%
  dplyr::group_by(grupo_edad) %>%
  dplyr::summarise(pob_std = sum(poblacion), .groups = "drop") %>%
  dplyr::mutate(peso = pob_std / sum(pob_std)) %>%
  dplyr::select(grupo_edad, peso)

tasa_pais_cvd_sexo <- df_age_mun_sexo %>%
  dplyr::group_by(ANIO, SEXO, grupo_edad) %>%
  dplyr::summarise(
    defunciones = sum(def_cvd),
    poblacion = sum(poblacion),
    .groups = "drop"
  ) %>%
  dplyr::left_join(pesos_2010, by = "grupo_edad") %>%
  dplyr::mutate(
    tasa_esp = defunciones / poblacion,
    tasa_ponderada = tasa_esp * peso
  ) %>%
  dplyr::group_by(ANIO, SEXO) %>%
  dplyr::summarise(
    tasa_ajustada = sum(tasa_ponderada, na.rm = TRUE) * 1e5,
    .groups = "drop"
  )


# ---- 1) Mapeo EDAD_AGRU (INEGI) -> grupos quinquenales tipo POB_*
map_edadagru_a_pobgrupo <- function(edad_agru){
  x <- suppressWarnings(as.integer(as.character(edad_agru)))
  out <- dplyr::case_when(
    x %in% 1:5   ~ "POB_00_04",
    x == 6       ~ "POB_05_09",
    x == 7       ~ "POB_10_14",
    x == 8       ~ "POB_15_19",
    x == 9       ~ "POB_20_24",
    x == 10      ~ "POB_25_29",
    x == 11      ~ "POB_30_34",
    x == 12      ~ "POB_35_39",
    x == 13      ~ "POB_40_44",
    x == 14      ~ "POB_45_49",
    x == 15      ~ "POB_50_54",
    x == 16      ~ "POB_55_59",
    x == 17      ~ "POB_60_64",
    x == 18      ~ "POB_65_69",
    x == 19      ~ "POB_70_74",
    x == 20      ~ "POB_75_79",
    x == 21      ~ "POB_80_84",
    x %in% 22:29 ~ "POB_85_mm",   # 85+ (85-89, 90-94,...,120)
    TRUE         ~ NA_character_  # 30 = No especificada y otros
  )
  out
}

# ---- 2) Niveles de edad esperados en CONAPO quinquenal
niveles_edad <- c(
  "POB_00_04","POB_05_09","POB_10_14","POB_15_19","POB_20_24",
  "POB_25_29","POB_30_34","POB_35_39","POB_40_44","POB_45_49",
  "POB_50_54","POB_55_59","POB_60_64","POB_65_69","POB_70_74",
  "POB_75_79","POB_80_84","POB_85_mm"
)

# ---- 3) Pesos estándar 2010 (ambos sexos combinados)
pesos_2010 <- poblacion_conapo_quinquenal %>%
  dplyr::filter(AÑO == 2010) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(niveles_edad),
    names_to = "grupo_edad",
    values_to = "poblacion"
  ) %>%
  dplyr::group_by(grupo_edad) %>%
  dplyr::summarise(pob_std = sum(as.numeric(poblacion), na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE),
    peso = pob_std / sum(pob_std, na.rm = TRUE)
  ) %>%
  dplyr::select(grupo_edad, peso)

# ---- 4) Función: calcula tasas (cruda y ajustada) para un año, sexo y tipo (Total/CVD)
calc_tasas_anio_sexo <- function(anio, sexo_label, tipo = c("Total","CVD")){
  tipo <- match.arg(tipo)
  
  # (a) Defunciones (numerador) por grupo_edad
  df_def <- classified_deaths_id[[as.character(anio)]]
  if (is.null(df_def)) return(NULL)
  
  df_def <- df_def %>%
    dplyr::transmute(
      ANIO = as.integer(anio),
      sexo = dplyr::case_when(
        suppressWarnings(as.integer(as.character(SEXO))) == 1 ~ "HOMBRES",
        suppressWarnings(as.integer(as.character(SEXO))) == 2 ~ "MUJERES",
        TRUE ~ NA_character_
      ),
      grupo_edad = map_edadagru_a_pobgrupo(EDAD_AGRU),
      muerte_cvd = suppressWarnings(as.integer(as.character(muerte_cvd)))
    ) %>%
    dplyr::filter(!is.na(sexo), !is.na(grupo_edad))
  
  if (tipo == "CVD") {
    df_def <- df_def %>% dplyr::filter(muerte_cvd == 1)
  }
  # si Total, no filtramos
  
  def_age <- df_def %>%
    dplyr::filter(sexo == sexo_label) %>%
    dplyr::count(grupo_edad, name = "defunciones") %>%
    dplyr::mutate(grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE))
  
  # (b) Población (denominador) por grupo_edad
  pob_age <- poblacion_conapo_quinquenal %>%
    dplyr::filter(AÑO == anio) %>%
    dplyr::mutate(
      sexo = dplyr::case_when(
        toupper(as.character(SEXO)) %in% c("HOMBRES","MUJERES") ~ toupper(as.character(SEXO)),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(sexo == sexo_label) %>%
    dplyr::summarise(
      POB_TOTAL = sum(as.numeric(POB_TOTAL), na.rm = TRUE),
      dplyr::across(dplyr::all_of(niveles_edad), ~ sum(as.numeric(.x), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(niveles_edad),
      names_to = "grupo_edad",
      values_to = "poblacion"
    ) %>%
    dplyr::mutate(grupo_edad = factor(grupo_edad, levels = niveles_edad, ordered = TRUE))
  
  # (c) Unir y calcular tasas específicas y ajustada
  df_age <- pob_age %>%
    dplyr::left_join(def_age, by = "grupo_edad") %>%
    dplyr::mutate(defunciones = dplyr::coalesce(defunciones, 0L)) %>%
    dplyr::left_join(pesos_2010, by = "grupo_edad") %>%
    dplyr::mutate(
      tasa_esp = dplyr::if_else(poblacion > 0, defunciones / poblacion, NA_real_)
    )
  
  tasa_ajustada <- sum(df_age$peso * df_age$tasa_esp, na.rm = TRUE) * 1e5
  def_total <- sum(df_age$defunciones, na.rm = TRUE)
  pob_total <- sum(df_age$poblacion, na.rm = TRUE)
  tasa_cruda <- ifelse(pob_total > 0, def_total / pob_total * 1e5, NA_real_)
  
  dplyr::tibble(
    ANIO = as.integer(anio),
    sexo = sexo_label,
    tipo = tipo,
    defunciones = as.integer(def_total),
    poblacion = as.numeric(pob_total),
    tasa_cruda = as.numeric(tasa_cruda),
    tasa_ajustada = as.numeric(tasa_ajustada)
  )
}

# ---- 5) Correr para todos los años y ambos sexos, para Total y CVD
tasas_long <- purrr::map_dfr(anios, function(anio){
  bind_rows(
    calc_tasas_anio_sexo(anio, "HOMBRES", "CVD"),
    calc_tasas_anio_sexo(anio, "MUJERES", "CVD")
  )
})

# ---- 6) Pasar a formato final: una fila por año, columnas por sexo y tipo
tasas_estratificacion <- tasas_long %>%
  tidyr::pivot_wider(
    names_from = c(tipo, sexo),
    values_from = c(defunciones, poblacion, tasa_cruda, tasa_ajustada),
    names_glue = "{.value}_{tolower(tipo)}_{tolower(sexo)}"
  ) %>%
  dplyr::arrange(ANIO)

View(tasas_estratificacion)

ggplot(filter(tasas_estratificacion, variable == "Sexo"),
       aes(ANIO, tasa_adj_100k, group = nivel_lbl, linetype = nivel_lbl)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Año", y = "Tasa ajustada por edad (por 100,000)", linetype = "") +
  theme_minimal(base_size = 12)


#Asignar SLI
# Asignar SLI/Rezago social a mortalidad municipal (por municipio-año)
library(dplyr)
library(stringr)

# 1) Normalizar clave municipal a 5 dígitos (ENT 2 + MUN 3)
norm_cvegeo5 <- function(x){
  x <- as.character(x)
  x <- str_extract(x, "\\d+")        # quedarte con dígitos
  x <- str_pad(x, 5, pad = "0")      # mínimo 5 (rellena a la izquierda)
  
  # Si viene más largo (6,7,...): asumir "ENT..." y "MUN" al final -> ENT = primeros 2, MUN = últimos 3
  too_long <- !is.na(x) & nchar(x) > 5
  x[too_long] <- paste0(
    substr(x[too_long], 1, 2),
    substr(x[too_long], nchar(x[too_long]) - 2, nchar(x[too_long]))
  )
  
  x
}

# 2) Preparar rezago_df (clave + año) y asegurar unicidad para evitar duplicar filas en el join
rezago_key <- rezago_df %>%
  mutate(
    cvegeo5 = norm_cvegeo5(cvegeo),
    anio_rezago = as.integer(anio_rezago)
  ) %>%
  select(cvegeo5, anio_rezago, indice_rezago, grado_rezago) %>%
  distinct(cvegeo5, anio_rezago, .keep_all = TRUE)

# 3) Preparar mortalidad (clave + año) y asignar el año de rezago más cercano
anios_rezago <- sort(unique(rezago_key$anio_rezago))

anio_mas_cercano <- function(y, candidatos){
  dif <- abs(candidatos - y)
  idx <- which(dif == min(dif))
  max(candidatos[idx])   # en empate, usa el más reciente (si prefieres el más viejo, cambia a min())
}

mortalidad_mun2 <- mortalidad_mun %>%
  mutate(
    cvegeo5 = norm_cvegeo5(cvegeo),
    ANIO = as.integer(ANIO),
    anio_rezago = vapply(ANIO, anio_mas_cercano, integer(1), candidatos = anios_rezago)
  )

# 4) Join final (municipio-año) para asignar índice/grado de rezago
mortalidad_mun_final <- mortalidad_mun2 %>%
  left_join(rezago_key, by = c("cvegeo5","anio_rezago"))

# 5) Diagnósticos rápidos
# Longitud de claves (evita NA con na.rm)
rezago_key %>%
  summarise(min=min(nchar(cvegeo5), na.rm=TRUE),
            max=max(nchar(cvegeo5), na.rm=TRUE),
            n_na=sum(is.na(cvegeo5)))

mortalidad_mun2 %>%
  summarise(min=min(nchar(cvegeo5), na.rm=TRUE),
            max=max(nchar(cvegeo5), na.rm=TRUE),
            n_na=sum(is.na(cvegeo5)))

# Cobertura global
mortalidad_mun_final %>%
  summarise(n_total=n(),
            n_sin_rezago=sum(is.na(indice_rezago)),
            prop_sin_rezago=mean(is.na(indice_rezago)))

# Cobertura por año
mortalidad_mun_final %>%
  group_by(ANIO) %>%
  summarise(n=n(),
            n_sin_rezago=sum(is.na(indice_rezago)),
            prop_sin_rezago=mean(is.na(indice_rezago))) %>%
  arrange(desc(prop_sin_rezago))

# Claves municipales únicas en mortalidad
mortalidad_mun2 %>% summarise(n_cvegeo5 = n_distinct(cvegeo5))

# Ejemplos de claves que NO existen en rezago
mortalidad_mun2 %>%
  distinct(cvegeo5) %>%
  anti_join(rezago_key %>% distinct(cvegeo5), by="cvegeo5") %>%
  slice(1:20)


# 11) (Opcional) Si quieres crear una categoría estándar SLI (ejemplo)
# Ajusta los cortes/labels a tu criterio
mortalidad_mun_final <- mortalidad_mun_final %>%
  mutate(
    sli_std = case_when(
      grado_rezago %in% c("Muy bajo", "Bajo") ~ "SLI bajo",
      grado_rezago == "Medio"                ~ "SLI medio",
      grado_rezago %in% c("Alto", "Muy alto")~ "SLI alto",
      TRUE                                   ~ NA_character_
    )
  )

# Tablas rápidas de distribución
mortalidad_mun_final %>%
  count(ANIO, sli_std) %>%
  group_by(ANIO) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(ANIO, sli_std)

names(mortalidad_mun_final)

library(dplyr)

mortalidad_mun_final <- mortalidad_mun_final %>%
  mutate(
    indice_rezago_ok = `indice_rezago.y.y.y`,
    grado_rezago_ok  = `grado_rezago.y.y.y`
  ) %>%
  select(
    -matches("^indice_rezago(\\.|$)"),
    -matches("^grado_rezago(\\.|$)")
  ) %>%
  rename(
    indice_rezago = indice_rezago_ok,
    grado_rezago  = grado_rezago_ok
  )

grep("^indice_rezago|^grado_rezago", names(mortalidad_mun_final), value = TRUE)

names(mortalidad_mun_final)

#Excess mortality 
mortalidad_mun_final <- mortalidad_mun_final %>%
  mutate(
    periodo = case_when(
      ANIO %in% 2015:2019 ~ "pre",
      ANIO %in% 2020:2023 ~ "pandemia",
      TRUE ~ NA_character_
    )
  )

ref_pre <- mortalidad_mun_final %>%
  filter(periodo == "pre") %>%
  group_by(cvegeo) %>%
  summarise(
    tasa_ref = mean(tasa_ajustada_edad_mun, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(tasa_ref), tasa_ref > 0)

MORTALITY <- mortalidad_mun_final %>%
  filter(ANIO %in% 2015:2023) %>%   # 🔑 filtrar aquí
  mutate(
    periodo = case_when(
      ANIO %in% 2015:2019 ~ "pre",
      ANIO %in% 2020:2023 ~ "pandemia"
    )
  ) %>%
  left_join(ref_pre, by = "cvegeo") %>%
  mutate(
    EXCESS_MORTALITY = if_else(periodo == "pandemia",
                               tasa_ajustada_edad_mun - tasa_ref, NA_real_),
    PERCENT_EXCESS_MORTALITY = if_else(periodo == "pandemia",
                                       100 * (tasa_ajustada_edad_mun - tasa_ref) / tasa_ref, NA_real_),
    RR_EXCESS = if_else(periodo == "pandemia",
                        tasa_ajustada_edad_mun / tasa_ref, NA_real_)
  )

table(MORTALITY$periodo, useNA = "ifany")
MORTALITY_PAN <- MORTALITY %>%
  filter(periodo == "pandemia")

#Crear tabla
MORTALITY_PAN %>%
  group_by(ANIO) %>%
  summarise(
    med = median(EXCESS_MORTALITY, na.rm = TRUE),
    q1  = quantile(EXCESS_MORTALITY, 0.25, na.rm = TRUE),
    q3  = quantile(EXCESS_MORTALITY, 0.75, na.rm = TRUE)
  )

library(ggplot2)

#Figura 1
ggplot(
  MORTALITY_PAN %>%
    group_by(ANIO) %>%
    summarise(
      med = median(EXCESS_MORTALITY, na.rm = TRUE),
      q1  = quantile(EXCESS_MORTALITY, 0.25, na.rm = TRUE),
      q3  = quantile(EXCESS_MORTALITY, 0.75, na.rm = TRUE)
    ),
  aes(x = ANIO, y = med)
) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    y = "Excess mortality (per 100,000)",
    x = "Year"
  ) +
  theme_minimal()


#Figura 2 
nat_obs <- mortalidad_mun_final %>%
  filter(ANIO %in% 2015:2023) %>%
  group_by(ANIO) %>%
  summarise(
    tasa_obs = weighted.mean(tasa_ajustada_edad_mun, w = poblacion_total, na.rm = TRUE),
    pob_nat  = sum(poblacion_total, na.rm = TRUE),
    .groups = "drop"
  )

## Colapsar a periodos 
nat_periodos <- nat_obs %>%
  mutate(
    periodo = case_when(
      ANIO %in% 2015:2019 ~ "2015–2019",
      ANIO %in% 2020:2023 ~ "2020–2023",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(periodo)) %>%
  group_by(periodo) %>%
  summarise(
    tasa_obs = weighted.mean(tasa_obs, w = pob_nat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(periodo = factor(periodo, levels = c("2015–2019","2020–2023")))

tasa_ref <- nat_periodos %>%
filter(periodo == "2015–2019") %>%
pull(tasa_obs) %>%
first()

df_plot <- nat_periodos %>%
mutate(
tasa_esperada = tasa_ref,
exceso = pmax(tasa_obs - tasa_esperada, 0)
) %>%
select(periodo, tasa_esperada, exceso) %>%
pivot_longer(c(tasa_esperada, exceso),
names_to = "componente",
values_to = "tasa")

ggplot(df_plot, aes(periodo, tasa, fill = componente)) +
  geom_col(width = 0.6) +
  scale_fill_manual(
    values = c(tasa_esperada = "#9ecae1", exceso = "#de2d26"),
    labels = c(tasa_esperada = "Tasa esperada (2015–2019)",
               exceso = "Exceso de mortalidad")
  ) +
  labs(x = NULL, y = "Tasa ajustada por edad", fill = NULL) +
  theme_minimal(base_size = 13)


# Modelo binominal negativo
Dataset.glm <- mortalidad_mun %>%
  filter(ANIO <= 2019)

Dataset.glm <- mortalidad_mun %>%
  mutate(
    # ---- Crear variables base ----
    EDU_NONLOW = medium_education_cv + high_education_cv,
    EMP_MH = ocu_high_cv + ocu_medium_cv,
    EMP_UL = ocu_unemployed_cv + ocu_low_cv,
    
    # ---- Crear ratios ----
    RATIO_SEX = (defunciones_cv_hombres) / (defunciones_cv_mujeres),
    RATIO_INDIGENOUS = (indigena_cv) / (no_indigena_cv),
    RATIO_ASISTENCIA = (asistencia_no_cv) / (asistencia_si_cv),
    RATIO_AMBULATORIO_HOST = (lugar_ambulatory_cv) / (lugar_hospital_cv),
    RATIO_URBANO_RURAL = (urbana_cv) / (rural_cv),
    RATIO_EDU_INF = (low_education_cv) / (EDU_NONLOW),
    RATIO_UL_JOB = (EMP_UL) / (EMP_MH),
    RATIO_LOW_JOB = (ocu_low_cv) / (EMP_MH),
    RATIO_SS = (no_derechohabiencia_cv) / (derechohabiencia_cv)
  )


# Reemplazar Inf por NA 
Dataset.glm$RATIO_SEX[Dataset.glm$RATIO_SEX == Inf] <- NA
Dataset.glm$RATIO_INDIGENOUS[Dataset.glm$RATIO_INDIGENOUS == Inf] <- NA
Dataset.glm$RATIO_ASISTENCIA[Dataset.glm$RATIO_ASISTENCIA == Inf] <- NA
Dataset.glm$RATIO_AMBULATORIO_HOST[Dataset.glm$RATIO_AMBULATORIO_HOST == Inf] <- NA
Dataset.glm$RATIO_URBANO_RURAL[Dataset.glm$RATIO_URBANO_RURAL == Inf] <- NA
Dataset.glm$RATIO_UL_JOB[Dataset.glm$RATIO_UL_JOB == Inf] <- NA
Dataset.glm$RATIO_LOW_JOB[Dataset.glm$RATIO_LOW_JOB == Inf] <- NA
Dataset.glm$RATIO_EDU_INF[Dataset.glm$RATIO_EDU_INF == Inf] <- NA
Dataset.glm$RATIO_SS[Dataset.glm$RATIO_SS == Inf] <- NA

mod_1 <- MASS::glm.nb(
  tasa_adj_mun_cvd ~
    scale(RATIO_SEX) + scale (RATIO_ASISTENCIA) + 
    scale (RATIO_AMBULATORIO_HOST) + 
    scale (RATIO_UL_JOB) + scale (RATIO_EDU_INF) + scale (RATIO_SS),
  data = Dataset.glm
)

jtools::summ(mod_1, exp = TRUE)

mod_1.5 <- MASS::glm.nb(
  tasa_adj_mun_cvd ~
    scale(RATIO_SEX) +
    scale(RATIO_ASISTENCIA) +
    scale(RATIO_INDIGENOUS) +
    scale(RATIO_AMBULATORIO_HOST) +
    scale(RATIO_URBANO_RURAL) +
    scale(RATIO_UL_JOB) +
    scale(RATIO_EDU_INF) +
    scale(RATIO_SS),
  data = Dataset.glm
)
jtools::summ(mod_1.5, exp = TRUE)


#Modelo 2
Dataset.glm_2 <- mortalidad_mun %>%
  dplyr::filter(ANIO >= 2012 & ANIO <= 2019)

Dataset.glm_2 <- Dataset.glm_2 %>%
  mutate(
    # ---- Crear variables base ----
    EDU_NONLOW = medium_education_cv + high_education_cv,
    EMP_MH = ocu_high_cv + ocu_medium_cv,
    EMP_UL = ocu_unemployed_cv + ocu_low_cv,
    
    # ---- Crear ratios ----
    RATIO_SEX = (defunciones_cv_hombres) / (defunciones_cv_mujeres),
    RATIO_INDIGENOUS = (indigena_cv) / (no_indigena_cv),
    RATIO_ASISTENCIA = (asistencia_no_cv) / (asistencia_si_cv),
    RATIO_AMBULATORIO_HOST = (lugar_ambulatory_cv) / (lugar_hospital_cv),
    RATIO_URBANO_RURAL = (urbana_cv) / (rural_cv),
    RATIO_EDU_INF = (low_education_cv) / (EDU_NONLOW),
    RATIO_UL_JOB = (EMP_UL) / (EMP_MH),
    RATIO_LOW_JOB = (ocu_low_cv) / (EMP_MH),
    RATIO_SS = (no_derechohabiencia_cv) / (derechohabiencia_cv)
  )

# Reemplazar Inf por NA 
Dataset.glm_2$RATIO_SEX[Dataset.glm_2$RATIO_SEX == Inf] <- NA
Dataset.glm_2$RATIO_INDIGENOUS[Dataset.glm_2$RATIO_INDIGENOUS == Inf] <- NA
Dataset.glm_2$RATIO_ASISTENCIA[Dataset.glm_2$RATIO_ASISTENCIA == Inf] <- NA
Dataset.glm_2$RATIO_AMBULATORIO_HOST[Dataset.glm_2$RATIO_AMBULATORIO_HOST == Inf] <- NA
Dataset.glm_2$RATIO_URBANO_RURAL[Dataset.glm_2$RATIO_URBANO_RURAL == Inf] <- NA
Dataset.glm_2$RATIO_UL_JOB[Dataset.glm_2$RATIO_UL_JOB == Inf] <- NA
Dataset.glm_2$RATIO_LOW_JOB[Dataset.glm_2$RATIO_LOW_JOB == Inf] <- NA
Dataset.glm_2$RATIO_EDU_INF[Dataset.glm_2$RATIO_EDU_INF == Inf] <- NA
Dataset.glm_2$RATIO_SS[Dataset.glm_2$RATIO_SS == Inf] <- NA



mod_2 <- MASS::glm.nb(
  tasa_adj_mun_cvd ~
    scale(RATIO_SEX) +
    scale(RATIO_ASISTENCIA) +
    scale(RATIO_INDIGENOUS) +
    scale(RATIO_AMBULATORIO_HOST) +
    scale(RATIO_URBANO_RURAL) +
    scale(RATIO_UL_JOB) +
    scale(RATIO_EDU_INF) +
    scale(RATIO_SS),
  data = Dataset.glm_2
)

summary(mod_2)
jtools::summ(mod_2, exp = TRUE)

#Modelo 3
Dataset.glm_3 <- mortalidad_mun %>%
  mutate(
    TOTAL_EDU = low_education_cv + medium_education_cv + high_education_cv,
    TOTAL_JOB = ocu_high_cv + ocu_medium_cv + ocu_low_cv + ocu_unemployed_cv,
    TOTAL_ASIST = asistencia_no_cv + asistencia_si_cv,
    TOTAL_LOC = lugar_ambulatory_cv + lugar_hospital_cv,
    TOTAL_RES = urbana_cv + rural_cv,
    TOTAL_SS = no_derechohabiencia_cv + derechohabiencia_cv,
    
    PROP_LOW_EDU = low_education_cv / TOTAL_EDU,
    PROP_UNEMP_LOW = (ocu_unemployed_cv + ocu_low_cv) / TOTAL_JOB,
    PROP_NO_ASIST = asistencia_no_cv / TOTAL_ASIST,
    PROP_AMBULATORY = lugar_ambulatory_cv / TOTAL_LOC,
    PROP_URBAN = urbana_cv / TOTAL_RES,
    PROP_INDIGENOUS = indigena_cv / (indigena_cv + no_indigena_cv),
    PROP_NO_SS = no_derechohabiencia_cv / TOTAL_SS
  )

mod_3 <- MASS::glm.nb(
  tasa_adj_mun_cvd ~
    scale(PROP_LOW_EDU) +
    scale(PROP_UNEMP_LOW) +
    scale(PROP_NO_ASIST) +
    scale(PROP_AMBULATORY) +
    scale(PROP_URBAN) +
    scale(PROP_INDIGENOUS) +
    scale(PROP_NO_SS),
  data = Dataset.glm_3
)

summary(mod_3)
jtools::summ(mod_3, exp = TRUE)

#---- Modelo de regresión logística 1 ----
library(dplyr)
library(stringr)
library(rlang)

# 1) Agregar por año
year_sum <- mortalidad_mun %>%
  group_by(ANIO) %>%
  summarise(
    def_cvd    = sum(def_cvd, na.rm = TRUE),
    def_noncvd = sum(def_noncvd, na.rm = TRUE),
    across(ends_with("_all"), ~sum(.x, na.rm = TRUE)),
    across(ends_with("_cv"),  ~sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# 2) Detectar pares *_all y *_cv
all_vars <- names(year_sum)[str_ends(names(year_sum), "_all")]
cv_vars  <- names(year_sum)[str_ends(names(year_sum), "_cv")]

roots_all <- str_remove(all_vars, "_all$")
roots_cv  <- str_remove(cv_vars,  "_cv$")

common_roots <- intersect(roots_all, roots_cv)
if (length(common_roots) == 0) stop("No encontré pares *_all y *_cv en year_sum.")

cv_cols  <- paste0(common_roots, "_cv")
all_cols <- paste0(common_roots, "_all")

# 3) Construir expresiones para NO CVD: (root_all - root_cv)
noncv_exprs <- setNames(
  lapply(common_roots, function(r) {
    expr( !!sym(paste0(r, "_all")) - !!sym(paste0(r, "_cv")) )
  }),
  common_roots
)

# 4) Base con 2 filas por año (CVD / no CVD)
base_cvd <- year_sum %>%
  select(ANIO, def_cvd, all_of(cv_cols)) %>%
  mutate(
    Def_CVD = 1L,
    deaths  = def_cvd
  ) %>%
  select(ANIO, Def_CVD, deaths, all_of(cv_cols))

names(base_cvd) <- c("ANIO", "Def_CVD", "deaths", common_roots)

base_noncvd <- year_sum %>%
  select(ANIO, def_noncvd, all_of(all_cols), all_of(cv_cols)) %>%
  mutate(
    Def_CVD = 0L,
    deaths  = def_noncvd
  ) %>%
  transmute(
    ANIO, Def_CVD, deaths,
    !!!noncv_exprs
  )

base_2rows <- bind_rows(base_cvd, base_noncvd) %>%
  arrange(ANIO, desc(Def_CVD))

#Education
base_agrup_edu <- df_individual %>%
  group_by(muerte_cvd) %>%
  summarise(
    deaths  = n(),
    edu_low = sum(NIVEL_ESCOLARIDAD == "low_education", na.rm = TRUE),
    .groups = "drop"
  )

#Education
glm(
  cbind(low_education, deaths - low_education) ~ Def_CVD,
  family = binomial,
  data = base_2rows
)


#SS
glm(cbind(no_derechohabiencia, deaths - no_derechohabiencia) ~ Def_CVD,
    family = binomial,
    data = base_2rows)

#Ocupacion
glm(cbind(ocu_unemployed, deaths - ocu_unemployed) ~ Def_CVD,
    family = binomial,
    data = base_2rows)

#Asistencia
glm(cbind(asistencia_no, deaths - asistencia_no) ~ Def_CVD,
    family = binomial,
    data = base_2rows)

#Urbana/rural
glm(cbind(urbana, deaths - urbana) ~ Def_CVD,
    family = binomial,
    data = base_2rows)

#Lugar/contexto
glm(cbind(lugar_ambulatory, deaths - lugar_ambulatory) ~ Def_CVD,
    family = binomial,
    data = base_2rows)

#Indigenous
glm(cbind(indigena, deaths - indigena) ~ Def_CVD,
    family = binomial,
    data = base_2rows)

#Modelo multivariable
glm(cbind(low_education, deaths - low_education) ~ 
      Def_CVD + lugar_ambulatory + ocu_unemployed + urbana + indigena + asistencia_no,
    family = binomial,
    data = base_2rows)


####
install.packages("speedglm")
library(speedglm)

df_individual_30 <- df_individual %>%
  slice_sample(prop = .50)

m1 <- speedglm(muerte_cvd ~ NIVEL_OCUPACION, data=df_individual_30, family=binomial())
summary(m1)

table(df_individual_1$NIVEL_OCUPACION)

m4 <- speedglm(muerte_cvd ~ NIVEL_OCUPACION + (1 | ANIO), data = df_individual_1, family = binomial)
jtools::summ(m4,exp=T)


#---- Modelo de regresión logística 2 ----
df_individual <- imap_dfr(classified_deaths_id, function(df_base, year) {
  
  year <- as.character(year)
  
  # Validar que todas las listas tengan mismo nrow
  stopifnot(
    nrow(df_base) == nrow(regis_edu_all[[year]]),
    nrow(df_base) == nrow(regis_dh_all[[year]]),
    nrow(df_base) == nrow(regis_lengua_all[[year]]),
    nrow(df_base) == nrow(regis_ocup_all[[year]]),
    nrow(df_base) == nrow(regis_asist_all[[year]]),
    nrow(df_base) == nrow(regis_area_all[[year]]),
    nrow(df_base) == nrow(contexto_classified_deaths_id[[year]])
  )
  
  tibble(
    cvegeo = sprintf("%05d", as.integer(df_base$cvegeo)),
    ANIO = as.integer(year),
    SEXO = suppressWarnings(as.integer(as.character(df_base$SEXO))),
    EDAD = suppressWarnings(as.numeric(df_base$EDAD)),
    muerte_cvd = suppressWarnings(as.integer(as.character(df_base$muerte_cvd))),
    
    NIVEL_ESCOLARIDAD = regis_edu_all[[year]]$NIVEL_ESCOLARIDAD,
    seguro = regis_dh_all[[year]]$seguro,
    LENGUA_INDIGENA = regis_lengua_all[[year]]$LENGUA_INDIGENA,
    NIVEL_OCUPACION = regis_ocup_all[[year]]$NIVEL_OCUPACION,
    ASISTENCIA = regis_asist_all[[year]]$ASISTENCIA,
    area = regis_area_all[[year]]$area,
    LUGAR_MUERTE = contexto_classified_deaths_id[[year]]$LUGAR_MUERTE
  )
})

df_individual <- df_individual %>%
  mutate(
    muerte_cvd = as.integer(muerte_cvd) 
    )

df_individual <- df_individual %>%
  mutate(
    SEXO = factor(SEXO),
    NIVEL_ESCOLARIDAD = factor(NIVEL_ESCOLARIDAD),
    NIVEL_OCUPACION = factor(NIVEL_OCUPACION),
    area = factor(area),
    LUGAR_MUERTE = factor(LUGAR_MUERTE),
    ANIO = factor(ANIO)
  )

dim(df_individual)
colSums(is.na(df_individual))
table(df_individual$muerte_cvd)

#Modelo por periodos
df_individual <- df_individual %>%
  mutate(
    ANIO_num = as.integer(as.character(ANIO)),
    quinquenio = case_when(
      ANIO_num >= 1990 & ANIO_num <= 1994 ~ "1990_1994",
      ANIO_num >= 1995 & ANIO_num <= 1999 ~ "1995_1999",
      ANIO_num >= 2000 & ANIO_num <= 2004 ~ "2000_2004",
      ANIO_num >= 2005 & ANIO_num <= 2009 ~ "2005_2009",
      ANIO_num >= 2010 & ANIO_num <= 2014 ~ "2010_2014",
      ANIO_num >= 2015 & ANIO_num <= 2019 ~ "2015_2019",
      ANIO_num >= 2020 & ANIO_num <= 2024 ~ "2020_2024"
    )
  )

#Modelo
install.packages("fixest")
library(fixest)

m_cluster <- feglm(
  muerte_cvd ~ 
    NIVEL_ESCOLARIDAD +
    seguro +
    LENGUA_INDIGENA +
    NIVEL_OCUPACION +
    ASISTENCIA +
    area +
    LUGAR_MUERTE +
    SEXO +
    scale(EDAD) +
    quinquenio,
  data = df_individual,
  family = binomial(),
  cluster = ~ cvegeo
)

summary(m_cluster)

