##########################
### LOAD LIBRARIES ######
##########################
library(plumber)
library(logger)
library(tictoc)
library(sf)
library(jsonlite)
library(dplyr)
library(rgbif)
library(geojsonsf)


sf::sf_use_s2(TRUE)

##########################
### WORKING DIRECTORY ###
##########################
wd <- "/media/docker/sRedList/sredlist-server"
if (!dir.exists(wd)) stop("Working directory not found: ", wd)
setwd(wd)

##########################
### LOGGING SETUP #######
##########################
logfile <- file.path("logs", "rina_api.log")
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
log_appender(appender_tee(logfile))

##########################
### SOURCE HELPERS ######
##########################
safe_source <- function(file, lines = NULL) {
  if (!file.exists(file)) stop("File not found: ", file)
  
  if (!is.null(lines)) {
    txt <- readLines(file)[lines]
    con <- textConnection(txt)
    source(con)
    close(con)
  } else {
    source(file, local = .GlobalEnv)
  }
}

safe_source("sRLfun_ShinyEditPoints.R")
safe_source("sRLfun_ShinyEditPolyg.R")
safe_source("sRLfun_ShinyDD.R")
safe_source("server.R", lines = 1:73)


############################
### SF → GEOJSON HELPER ###
############################
sf_to_geojson_clean <- function(sf_obj) {
  if (!inherits(sf_obj, "sf")) stop("Object must be of class 'sf'")
  
  # Transform to WGS84 and make valid geometries
  sf_obj <- sf::st_transform(sf::st_make_valid(sf_obj), 4326)
  
  # Convert to GeoJSON
  geojson_text <- geojsonsf::sf_geojson(sf_obj)
  geojson_list <- jsonlite::fromJSON(geojson_text, simplifyVector = FALSE)
  
  # Feature Cleaning
  geojson_list$features <- lapply(geojson_list$features, function(feat) {
    feat$type <- "Feature"
    
    # Property cleanup: scalars remain scalars, empty lists become NULL
    feat$properties <- lapply(feat$properties, function(x) {
      if (is.null(x) || length(x) == 0) return(NULL)
      if (is.list(x) && length(x) == 1) return(x[[1]])
      if (length(x) == 1) return(x[[1]])
      x
    })
    
    # Geometry type as a string
    feat$geometry$type <- as.character(feat$geometry$type)
    
    # Recursive flatten coordinates
    flatten_coords <- function(coords) {
      if (is.list(coords)) {
        if (all(sapply(coords, is.numeric))) {
          return(as.numeric(coords))
        } else {
          return(lapply(coords, flatten_coords))
        }
      }
      coords
    }
    feat$geometry$coordinates <- flatten_coords(feat$geometry$coordinates)
    
    feat
  })
  
  geojson_list
}

###############################
### UNIFIED JSON NORMALIZER ###
###############################
.normalize_datetime <- function(x) {
  if (inherits(x, "POSIXt")) return(format(as.POSIXct(x, tz="UTC"), "%Y-%m-%dT%H:%M:%SZ"))
  if (inherits(x, "Date")) return(format(as.Date(x), "%Y-%m-%d"))
  x
}

.coerce_column <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  .normalize_datetime(x)
}

serialize_table_unified <- function(obj) {
  if (inherits(obj, "sf")) {
    df <- sf::st_drop_geometry(obj)
    geom <- tryCatch(sf::st_geometry_type(obj), error = function(e) NULL)
    if (!is.null(geom) && all(grepl("POINT", geom))) {
      coords <- sf::st_coordinates(obj)
      df$lon <- coords[,1]
      df$lat <- coords[,2]
    } else {
      df$geometry_wkt <- sf::st_as_text(sf::st_geometry(obj))
    }
  } else if (is.data.frame(obj)) {
    df <- obj
  } else {
    return(normalize_json_value(obj))
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  for (n in names(df)) df[[n]] <- .coerce_column(df[[n]])

  lapply(seq_len(nrow(df)), function(i) {
    lapply(df[i, , drop=FALSE], normalize_json_value)
  })
}

normalize_json_value <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  if (inherits(x, "sf") || is.data.frame(x)) return(serialize_table_unified(x))
  if (is.list(x)) return(lapply(x, normalize_json_value))
  if (length(x) == 1) return(jsonlite::unbox(.coerce_column(x)))
  x
}

normalize_json_table <- function(x) normalize_json_value(x)

##########################
### LOAD DDPRIO DATA ####
##########################
DD <- readRDS("resources/resources_Shiny_DD/DD_prepared_for_ShinyREALMS.rds")

##########################
### LOAD WOS (lazy) #####
##########################
load_WOS <- function() {
  # If WOS_extract is already in memory, do nothing
  if (exists("WOS_extract", envir = .GlobalEnv)) return(invisible(TRUE))
  
  # Check that the DDfun_RefWos function exists
  if (!exists("DDfun_RefWos", mode = "function")) {
    stop("The DDfun_RefWos function was not loaded. Please ensure that sRLfun_ShinyDD.R has been sourced correctly.")
  }
  
  # Check that the RDS file exists
  rds_file <- "resources/resources_Shiny_DD/WOS_extract_moreinfo.rds"
  if (!file.exists(rds_file)) stop("RDS file not found: ", rds_file)
  
  # Load and process data
  WOS_extract <<- tryCatch(
    {
      readRDS(rds_file) |>
        subset(!is.na(Title)) |>
        DDfun_RefWos()
    },
    error = function(e) {
      stop("Error reading or processing WOS_extract: ", e$message)
    }
  )
  
  invisible(TRUE)
}


################################### ENDPOINTS ###################################


##########################
### DDPRIO ENDPOINTS ####
##########################

#* Get DD priority by group
#* @param group Taxonomic group
#* @post /ddprio/by_group
function(req, res, group = req$argsQuery$group) {
  if (is.null(group) || group == "") {
    res$status <- 400
    return(list(error = "Missing group parameter"))
  }

  load_WOS()
  sub <- subset(DD, Group == group & DD)

  sub$inc.GBIF <- as.integer(round(100*(sub$nb_GBIFgeo - sub$nb_GBIFgeoASS)/(1+sub$nb_GBIFgeoASS)))
  sub$inc.WOS  <- as.integer(round(100*(sub$WOS - sub$WOSASS)/(1+sub$WOSASS)))

  sub$WOSList <- lapply(sub$scientific_name, function(sp) {
    refs <- subset(WOS_extract, Species == sp)$Ref
    if (length(refs) == 0) list() else as.list(refs)
  })

  list(group = group, count = nrow(sub), data = normalize_json_table(sub))
}

#* Get DD priority for all species
#* @post /ddprio/all
function() {
  load_WOS()
  sub <- subset(DD, DD)
  sub$WOSList <- lapply(sub$scientific_name, function(sp) {
    refs <- subset(WOS_extract, Species == sp)$Ref
    if (length(refs) == 0) list() else as.list(refs)
  })
  list(count = nrow(sub), data = normalize_json_table(sub))
}

##########################
### SPECIES ENDPOINTS ###
##########################

#* Get DD priority points for species
#* @param species Scientific name
#* @post /species/<species>/DDprio_points
function(species, res) {
  load_WOS()
  sub <- subset(DD, scientific_name == species & DD)
  if (nrow(sub) == 0) {
    res$status <- 404
    return(list(error = "Species not found"))
  }
  sub$WOSList <- subset(WOS_extract, Species == species)$Ref
  list(species = species, data = normalize_json_table(sub))
}

#* Get WOS references for species
#* @param species Scientific name
#* @post /species/<species>/WOS
function(species) {
  load_WOS()
  refs <- subset(WOS_extract, Species == species)$Ref
  list(species = species, WOS = as.list(refs))
}

##########################
### STORAGE / FLAGS #####
##########################

#* Read storage
#* @post /storage
function(req, res) {
  q <- req$argsQuery
  if (is.null(q$sci_name) || is.null(q$user)) {
    res$status <- 400
    return(list(error = "Missing sci_name or user"))
  }
  normalize_json_table(sRL_StoreRead(q$sci_name, q$user, MANDAT = 1))
}

#* Get flags
#* @post /flags
function(req, res) {
  q <- req$argsQuery
  sp <- sRL_StoreRead(q$sci_name, q$user, MANDAT = 1)
  normalize_json_table(sp$flags)
}

#* Get points
#* @post /points
function(req) {
  q <- req$argsQuery
  sp <- sRL_StoreRead(q$sci_name, q$user, MANDAT = 1)
  pts <- sp$flags
  if (inherits(pts, "sf")) pts <- sf::st_transform(pts, 4326)
  normalize_json_table(pts)
}

##########################
### NON-GEO GBIF ########
##########################

#* Get non-georeferenced GBIF records
#* @post /nongeo
function(req, res) {
  query    <- req$argsQuery
  sci_name <- query$sci_name

  if (is.null(sci_name)) {
    res$status <- 400
    return(list(error = "Please provide 'sci_name' parameter"))
  }

  gbif_data <- tryCatch(
    rgbif::occ_data(scientificName = sci_name, hasCoordinate = FALSE, limit = 1000)$data,
    error = function(e){ res$status <- 500; return(list(error = paste("GBIF request failed:", e$message))) }
  )

  if (is.list(gbif_data) && !is.data.frame(gbif_data)) return(gbif_data)
  if (is.null(nrow(gbif_data))) return(list(warning="No non-georeferenced records found", count=0, data=list()))

  if (!"country" %in% names(gbif_data))  gbif_data$country  <- NA
  if (!"locality" %in% names(gbif_data)) gbif_data$locality <- NA

  Tab <- gbif_data %>%
    dplyr::mutate(Link = paste0("https://gbif.org/occurrence/", gbifID)) %>%
    dplyr::select(
      scientificName, basisOfRecord, eventDate, higherGeography,
      continent, country, locality, institutionCode, collectionCode,
      occurrenceRemarks, identifiedBy, Link
    )

  Tab_summ <- gbif_data %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(Localities = paste(unique(na.omit(locality)), collapse="<br>"),
                     Number_records = dplyr::n(), .groups="drop")

  warn_limit <- nrow(gbif_data) == 1000

  list(
    warning_limit1000 = jsonlite::unbox(warn_limit),
    count   = jsonlite::unbox(nrow(Tab)),
    Raw     = normalize_json_table(Tab),
    Summary = normalize_json_table(Tab_summ)
  )
}

##############################
### COUNTRIES - LOAD DATA ###
##############################
#* Get species countries + distribution polygon
#* @get /species/<sci_name>/countries
function(req, res, sci_name, user = req$argsQuery$user) {

  # User parameter control
  if (missing(user) || user == "") {
    res$status <- 400
    return(list(error = "Missing user parameter"))
  }
  
  # Read data from storage
  Storage_SP <- tryCatch(
    sRL_StoreRead(sci_name, user, MANDAT = 1),
    error = function(e) NULL
  )
  
  if (is.null(Storage_SP)) {
    res$status <- 404
    return(list(error = "Species storage not found"))
  }
  
  COO <- Storage_SP$coo
  coo_occ <- Storage_SP$coo_occ
  
  # Align presence, origin, seasonality
  if (!is.null(coo_occ)) {
    COO$presence <- coo_occ$presence[match(COO$lookup, coo_occ$lookup)]
    COO$origin   <- coo_occ$origin[match(COO$lookup, coo_occ$lookup)]
    COO$seasonal <- coo_occ$seasonal[match(COO$lookup, coo_occ$lookup)]
  }
  
  # Level0_occupied
  COO$Level0_occupied <- COO$SIS_name0 %in% subset(COO, Level1_occupied == TRUE)$SIS_name0
  
  # GeoJSON of countries
  countries_geojson <- sf_to_geojson_clean(COO)
  
  # Table without geometry
  COO_table <- COO %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct(lookup, .keep_all = TRUE)
  
  # ----------------------------
  # Pre-calculated distribution polygon
  # ----------------------------
  dist_poly <- Storage_SP$distSP_saved
  if (!inherits(dist_poly, "sf")) {
    dist_poly <- NULL
    area_km2 <- 0
    dist_geojson <- NULL
  } else {
    # Let's make sure WGS84 and valid geometries
    dist_poly <- sf::st_transform(sf::st_make_valid(dist_poly), 4326)
    
    # Convert to JSON-ready R list (not string)
    dist_geojson <- jsonlite::fromJSON(geojsonsf::sf_geojson(dist_poly), simplifyVector = FALSE)
    
    # Calculate area km²
    area_km2 <- as.numeric(sf::st_area(dist_poly)) / 1e6
  }

  list(
    species = sci_name,
    countries_geojson = countries_geojson,
    distribution_polygon = dist_geojson,
    area_km2 = jsonlite::unbox(area_km2),
    table = normalize_json_table(COO_table)
  )
}

###################################
### COUNTRIES - SAVE #############
###################################

#* Save countries (update + save)
#* @post /species/<sci_name>/countries/save
function(req, res, sci_name, user) {
  
  body <- tryCatch(jsonlite::fromJSON(req$postBody, simplifyVector = FALSE),
                   error = function(e) NULL)
  
  if (is.null(body$changes) || length(body$changes) == 0) {
    res$status <- 400
    return(list(error = "Missing changes"))
  }

  Storage_SP <- tryCatch(sRL_StoreRead(sci_name, user, MANDAT = 1),
                         error = function(e) NULL)
  
  if (is.null(Storage_SP)) {
    res$status <- 404
    return(list(error = "Species storage not found"))
  }

  COO <- Storage_SP$coo

  # Logical columns to normalize
  logical_cols <- c("Level1_occupied", "Level0_occupied")

  for (chg in body$changes) {
    if (!is.list(chg) || is.null(chg$lookup) || is.null(chg$field) || is.null(chg$value)) next
    idx <- which(COO$lookup == chg$lookup)
    if (length(idx) == 0) next

    # Logical conversion if necessary
    if (chg$field %in% logical_cols) {
      COO[idx, chg$field] <- as.logical(as.numeric(chg$value))
    } else {
      COO[idx, chg$field] <- chg$value
    }
  }

  # Recalculate Level0_occupied consistently
  COO$Level0_occupied <- COO$SIS_name0 %in% subset(COO, Level1_occupied == TRUE)$SIS_name0

  # Update coo_occ
  Storage_SP$coo <- COO
  Storage_SP$coo_occ <- COO[, c("lookup", "presence", "origin", "seasonal")]

  # Save to disk
  tryCatch({
    sRL_StoreSave(sci_name, user, Storage_SP)
  }, error = function(e) {
    res$status <- 500
    return(list(error = paste("Save failed:", e$message)))
  })

  list(
    message = "Changes applied and saved",
    data = normalize_json_table(sf::st_drop_geometry(COO))
  )
}








##########################
### CORS FILTER #########
##########################
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "*")
  if (req$REQUEST_METHOD == "OPTIONS") return(list(status = "OK"))
  plumber::forward()
}

##########################
### LOGGING HOOKS #######
##########################
#* @hook preroute
function(req) { tictoc::tic() }

#* @hook postroute
function(req, res) {
  tm <- tictoc::toc(quiet = TRUE)
  log_info("{req$REQUEST_METHOD} {req$PATH_INFO} {res$status} {round(tm$toc - tm$tic,3)}s")
}
