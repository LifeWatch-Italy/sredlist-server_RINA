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
library(concaveman)
library(MASS)
library(smoothr)
library(rmapshaper)
library(lwgeom)
# shiny is needed because sRLPolyg_PrepareHydro -> sRLPolyg_CreatePopup builds
# popup HTML with actionButton/tooltip/icon. The popup itself is rebuilt on the
# Angular client, so we drop the R-generated Popup column before serializing.
suppressMessages(library(shiny))
suppressMessages(library(shinyWidgets))
suppressMessages(library(bslib))

sf::sf_use_s2(TRUE)

#* Return the real R error message instead of plumber's generic
#* "500 - Internal server error", so the client can show it to the user
#* @plumber
function(pr) {
  pr$setErrorHandler(function(req, res, err) {
    log_error("{req$PATH_INFO}: {conditionMessage(err)}")
    res$status <- 500
    list(status = "error", message = conditionMessage(err))
  })
}

# Global null-coalescing helper (some endpoints also define it locally)
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
safe_source <- function(file, upto_marker = NULL) {
  if (!file.exists(file)) stop("File not found: ", file)

  if (!is.null(upto_marker)) {
    all_lines <- readLines(file)
    end <- which(grepl(upto_marker, all_lines, fixed = TRUE))[1]
    if (is.na(end)) stop("Marker '", upto_marker, "' not found in ", file)
    con <- textConnection(all_lines[1:(end - 1)])
    source(con, local = .GlobalEnv)
    close(con)
  } else {
    source(file, local = .GlobalEnv)
  }
}

safe_source("sRLfun_ShinyEditPoints.R")
safe_source("sRLfun_ShinyEditPolyg.R")
safe_source("sRLfun_ShinyDD.R")
# Only load server.R's library()/source() setup (helper functions, config), not its
# plumber routes — stop right before the marker line so this stays correct even if
# server.R's line count changes, instead of a brittle hardcoded line range.
safe_source("server.R", upto_marker = "### RINA_API_SETUP_END")


############################
### HYDROBASINS LAYERS ####
############################
# Lazy-load the hydrobasins layers (level 8 is large, so we only read them on the
# first hydrobasins request rather than at API startup). Mirrors load_WOS().
# - hydro_raw : level 8 with attributes (hybas_id, next_down, next_sink), CRS = Mollweide
# - hydro3    : level 3 overlay (blue "Hydrobasins-3" layer)
# - distCountries_mapping : Red List countries (server.R:104), used by the load flow
load_hydro <- function() {
  if (exists("hydro_raw", envir = .GlobalEnv) &&
      exists("hydro3", envir = .GlobalEnv)) return(invisible(TRUE))

  log_info("Loading hydrobasins layers (first request)...")

  hydro_raw <<- tryCatch({
    h <- sf::st_read(config$hydrobasins_path, quiet = TRUE)
    sf::st_crs(h) <- CRSMOLL
    h
  }, error = function(e) stop("Failed to load hydro_raw: ", e$message))

  # Precompute the per-feature bounding boxes ONCE. st_filter over the full 234k-basin
  # layer rebuilds a spatial index on every call (~6.6s); a cheap numeric bbox pre-filter
  # against this cached matrix reduces the candidate set to a few hundred first (~0.05s).
  hydro_bbox <<- tryCatch(
    do.call(rbind, lapply(sf::st_geometry(hydro_raw), function(g) sf::st_bbox(g))),
    error = function(e) { log_warn("Failed to precompute hydro_bbox: {e$message}"); NULL }
  )

  hydro3 <<- tryCatch(
    sf::st_read(sub("level8_", "level3_", config$hydrobasins_path), quiet = TRUE),
    error = function(e) stop("Failed to load hydro3: ", e$message)
  )

  if (!exists("distCountries_mapping", envir = .GlobalEnv)) {
    distCountries_mapping <<- tryCatch({
      d <- sf::st_read("Species/Map countries/Red_List_countries_msSimplif0.05_MOLL.shp", quiet = TRUE)
      sf::st_crs(d) <- CRSMOLL
      d
    }, error = function(e) { log_warn("Failed to load distCountries_mapping: {e$message}"); NULL })
  }

  invisible(TRUE)
}


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

##############################
### POLYGON HELPERS #########
##############################

create_mcp <- function(points_sf) {
  st_convex_hull(st_union(points_sf))
}

create_alpha <- function(points_sf, concavity = 2) {
  concaveman(points_sf, concavity = concavity)
}

create_kernel <- function(df, level = 0.5, n = 100) {
  kde <- MASS::kde2d(df$lon, df$lat, n = n)
  
  lvl <- quantile(kde$z, level)
  cls <- contourLines(kde$x, kde$y, kde$z, levels = lvl)
  
  if (length(cls) == 0) return(NULL)
  
  polys <- lapply(cls, function(cl) {
    coords <- cbind(cl$x, cl$y)
    if (nrow(coords) < 4) return(NULL)
    st_polygon(list(coords))
  })
  
  polys <- polys[!sapply(polys, is.null)]
  if (length(polys) == 0) return(NULL)
  
  st_sfc(polys, crs = 4326)
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
  if (is.null(q$sci_name) || is.null(q$username)) {
    res$status <- 400
    return(list(error = "Missing sci_name or username"))
  }
  normalize_json_table(sRL_StoreRead(q$sci_name, q$username, MANDAT = 1))
}

#* Get flags
#* @post /flags
function(req, res) {
  q <- req$argsQuery
  sp <- sRL_StoreRead(q$sci_name, q$username, MANDAT = 1)
  normalize_json_table(sp$flags)
}

#* Get points
#* @post /points
function(req, res) {
  q <- req$argsQuery
  sci_name <- gsub("_", " ", utils::URLdecode(q$sci_name %||% ""))
  sp <- tryCatch(sRL_StoreRead(sci_name, q$username, MANDAT = 1),
                 error = function(e) { log_warn("points StoreRead failed: {e$message}"); NULL })
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
function(req, res, sci_name, username = req$argsQuery$username) {

  # sci_name arrives URL-encoded from the path (e.g. "Saara%20loricata"); normalize it
  sci_name <- gsub("_", " ", utils::URLdecode(sci_name))

  # User parameter control
  if (missing(username) || username == "") {
    res$status <- 400
    return(list(error = "Missing username parameter"))
  }
  
  # Read data from storage
  Storage_SP <- tryCatch(
    sRL_StoreRead(sci_name, username, MANDAT = 1),
    error = function(e) NULL
  )
  
  if (is.null(Storage_SP)) {
    res$status <- 404
    return(list(error = "Species storage not found"))
  }
  
  COO <- Storage_SP$coo
  coo_occ <- Storage_SP$coo_occ

  # The countries-of-occurrence table may not exist yet (COO step not run for this species)
  if (is.null(COO) || !("Level1_occupied" %in% names(COO))) {
    res$status <- 404
    return(list(error = "No countries of occurrence available for this species. Run the Countries (COO) step first."))
  }

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
    
    # Calculate area km² (sum across features so unbox() gets a scalar)
    area_km2 <- sum(as.numeric(sf::st_area(dist_poly)), na.rm = TRUE) / 1e6
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
function(req, res, sci_name, username) {

  # sci_name arrives URL-encoded from the path (e.g. "Loxodonta%20africana"); normalize it
  sci_name <- gsub("_", " ", utils::URLdecode(sci_name))

  body <- tryCatch(jsonlite::fromJSON(req$postBody, simplifyVector = FALSE),
                   error = function(e) NULL)
  
  if (is.null(body$changes) || length(body$changes) == 0) {
    res$status <- 400
    return(list(error = "Missing changes"))
  }

  Storage_SP <- tryCatch(sRL_StoreRead(sci_name, username, MANDAT = 1),
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
      # Empty string from the client means "cleared" → store NA, not "" (matches original)
      val <- chg$value
      if (is.character(val) && val == "") val <- NA
      COO[idx, chg$field] <- val
    }
  }

  # Vectorised "is TRUE" that treats NA/other as FALSE
  occupied <- COO$Level1_occupied %in% TRUE

  # When a country is not occupied, clear presence/origin/seasonal to NA (port of
  # Shiny_Countries.R:201) so downstream code and colours treat it as empty.
  empty_rows <- which(!occupied)
  if (length(empty_rows) > 0) {
    for (col in c("presence", "origin", "seasonal")) {
      if (col %in% names(COO)) COO[empty_rows, col] <- NA
    }
  }

  # Recalculate Level0_occupied consistently
  COO$Level0_occupied <- COO$SIS_name0 %in% subset(COO, Level1_occupied == TRUE)$SIS_name0

  # coo_occ = only the occupied countries, collapsing marine/terrestrial duplicates
  # (port of Shiny_Countries.R:223-225)
  coo_occ <- COO[COO$Level1_occupied %in% TRUE, ]
  cols_keep <- intersect(c("SIS_name0", "SIS_name1", "lookup", "lookup_SIS0", "presence", "origin", "seasonal"), names(coo_occ))
  coo_occ <- sf::st_drop_geometry(coo_occ)[, cols_keep, drop = FALSE]
  coo_occ <- dplyr::distinct(coo_occ, lookup, .keep_all = TRUE)

  Storage_SP$coo <- COO
  Storage_SP$coo_occ <- coo_occ

  # Regenerate the stored leaflet for the RMD report (best-effort)
  tryCatch({
    if (exists("sRLCountry_CreateLeaflet", mode = "function")) {
      Storage_SP$Leaflet_COO <- sRLCountry_CreateLeaflet(COO, Storage_SP, FALSE)
    }
  }, error = function(e) log_warn("Could not rebuild Leaflet_COO: {e$message}"))

  # Record usage counter COO_EditShiny (yes + count++)
  tryCatch({
    op <- Storage_SP$Output
    if (!is.null(op) && "COO_EditShiny" %in% op$Parameter) {
      op$Value[op$Parameter == "COO_EditShiny"] <- "yes"
      op$Count[op$Parameter == "COO_EditShiny"] <- as.numeric(op$Count[op$Parameter == "COO_EditShiny"]) + 1
      Storage_SP$Output <- op
    }
  }, error = function(e) log_warn("Could not update COO_EditShiny counter: {e$message}"))

  # Save to disk
  save_ok <- tryCatch({
    sRL_StoreSave(sci_name, username, Storage_SP)
    TRUE
  }, error = function(e) { res$status <- 500; log_error("COO save failed: {e$message}"); FALSE })
  if (!isTRUE(save_ok)) return(list(error = "Save failed"))

  list(
    message = "Changes applied and saved",
    data = normalize_json_table(sf::st_drop_geometry(COO))
  )
}







##########################
### POLYGONS ###########
##########################
#* Get points + polygons together (FINAL CLEAN VERSION)
#* @post /mapdata
function(req, res) {
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  q <- req$argsQuery
  sci_name <- q$sci_name
  username     <- q$username
  method   <- tolower(q$method %||% "mcp")
  
  if (is.null(sci_name) || is.null(username)) {
    res$status <- 400
    return(list(error = "Missing sci_name or username"))
  }
  
  log_info("MAPDATA START species={sci_name} username={username} method={method}")
  
  # ----------------------------
  # LOAD STORAGE
  # ----------------------------
  Storage_SP <- tryCatch(
    sRL_StoreRead(sci_name, username, MANDAT = 1),
    error = function(e) {
      log_error("Storage read failed: {e$message}")
      NULL
    }
  )
  
  if (is.null(Storage_SP)) {
    res$status <- 404
    return(list(error = "Storage not found"))
  }
  
  pts_raw <- Storage_SP$flags
  log_info("Flags class: {paste(class(pts_raw), collapse=', ')}")
  
  # ----------------------------
  # FORCE → SF
  # ----------------------------
  pts_sf <- NULL
  
  if (inherits(pts_raw, "sf")) {
    
    pts_sf <- sf::st_transform(pts_raw, 4326)
    
  } else if (is.data.frame(pts_raw)) {
    
    if (!all(c("decimalLongitude", "decimalLatitude") %in% names(pts_raw))) {
      log_warn("Missing coordinates in flags")
    } else {
      
      pts_clean_df <- pts_raw %>%
        dplyr::mutate(
          decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)),
          decimalLatitude  = suppressWarnings(as.numeric(decimalLatitude))
        ) %>%
        dplyr::filter(
          !is.na(decimalLongitude),
          !is.na(decimalLatitude)
        )
      
      if (nrow(pts_clean_df) > 0) {
        pts_sf <- sf::st_as_sf(
          pts_clean_df,
          coords = c("decimalLongitude", "decimalLatitude"),
          crs = 4326,
          remove = FALSE
        )
      }
    }
  } else {
    log_warn("Flags is neither sf nor data.frame")
  }
  
  # ----------------------------
  # NO VALID POINTS → keep going, polygons may still be available
  # ----------------------------
  no_points <- is.null(pts_sf) || nrow(pts_sf) == 0
  points_warning <- NULL

  if (no_points) {

    log_warn("No valid points after conversion")
    points_warning <- "No valid points"
    points_out <- list()
    pts_unique <- if (is.null(pts_sf)) data.frame() else pts_sf[0, ]

  } else {

    # ----------------------------
    # EXTRACT COORDS
    # ----------------------------
    coords <- sf::st_coordinates(pts_sf)
    pts_sf$lon <- coords[,1]
    pts_sf$lat <- coords[,2]

    # ----------------------------
    # FILTER VALID POINTS (SHINY-COMPATIBLE)
    # ----------------------------
    valid_cols <- c(".val", ".equ", ".zer", ".cap", ".cen", ".gbf")
    existing_cols <- valid_cols[valid_cols %in% names(pts_sf)]

    pts_filtered <- pts_sf

    if (length(existing_cols) > 0) {

      pts_filtered <- pts_filtered %>%
        dplyr::mutate(
          dplyr::across(
            all_of(existing_cols),
            ~ as.logical(as.character(.))
          )
        ) %>%
        dplyr::filter(
          dplyr::if_all(all_of(existing_cols), ~ . == TRUE)
        )

      log_info("Filtered valid points: {nrow(pts_filtered)}")

    } else {
      log_warn("No validation columns found → skipping filter")
    }

    # ----------------------------
    # REMOVE DUPLICATES
    # ----------------------------
    pts_unique <- pts_filtered %>%
      dplyr::distinct(lon, lat, .keep_all = TRUE)

    log_info("Points: total={nrow(pts_sf)} valid={nrow(pts_filtered)} unique={nrow(pts_unique)}")

    # ----------------------------
    # SERIALIZE POINTS
    # ----------------------------
    points_out <- normalize_json_table(pts_unique)
  }
  
  # ----------------------------
  # POLYGON GENERATION
  # ----------------------------
  poly_geojson <- list()
  warning_msg <- NULL
  
  dist_poly <- Storage_SP$distSP_saved
  
  # 🔹 USE PRECOMPUTED (Shiny logic)
  if (!is.null(dist_poly) && inherits(dist_poly, "sf")) {
    
    log_info("Using precomputed polygon from storage")
    
    # Same as sRLPolyg_InitDistri: split multipolygons into individual
    # polygons so each one gets its own popup/attributes on the client.
    dist_poly <- sf_polygons_only(sf::st_transform(dist_poly, 4326))
    poly_geojson <- sf_to_geojson_clean(dist_poly)
    
  } else {
    
    log_warn("No precomputed polygon → fallback")
    
    # 🔹 MCP fallback
    if (method == "mcp" && nrow(pts_unique) >= 3) {
      
      log_info("Fallback MCP computation")
      
      hull_sf <- tryCatch({
        
        pts_proj <- sf::st_transform(pts_unique, 3857)
        geom <- sf::st_convex_hull(sf::st_union(pts_proj))
        
        sf::st_as_sf(
          data.frame(id = 1),
          geometry = sf::st_transform(geom, 4326)
        )
        
      }, error = function(e) {
        log_error("Hull error: {e$message}")
        NULL
      })
      
      if (!is.null(hull_sf) && !sf::st_is_empty(hull_sf)) {
        poly_geojson <- sf_to_geojson_clean(hull_sf)
      } else {
        warning_msg <- "Polygon generation failed"
      }
      
    } else if (method == "individual") {
      
      log_info("Method = individual → no polygon")
      poly_geojson <- list()
      
    } else {
      
      warning_msg <- paste("Polygon not available for method:", method)
      log_warn(warning_msg)
    }
  }
  
  # ----------------------------
  # TEXT DEFAULTS for the attribute panel (port of Shiny_EditPoly.R:325-346)
  # ----------------------------
  text_defaults <- tryCatch({
    la <- if (!is.null(dist_poly) && inherits(dist_poly, "sf")) sf::st_drop_geometry(dist_poly) else data.frame()
    names(la) <- tolower(names(la))
    firstNonNull <- function(col, default) {
      if (col %in% names(la)) { v <- la[[col]][!is.na(la[[col]])]; if (length(v) > 0) return(as.character(v[1])) }
      default
    }
    list(
      source     = firstNonNull("source", "sRedList platform"),
      yrcompiled = firstNonNull("yrcompiled", format(Sys.time(), "%Y")),
      citation   = firstNonNull("citation", "IUCN (International Union for Conservation of Nature)"),
      compiler   = firstNonNull("compiler", tryCatch(sRL_userformatted(username), error = function(e) username)),
      island     = firstNonNull("island", ""),
      data_sens  = 0,
      sens_comm  = firstNonNull("sens_comm", ""),
      dist_comm  = firstNonNull("dist_comm", "")
    )
  }, error = function(e) list())

  # ----------------------------
  # FINAL RESPONSE
  # ----------------------------
  combined_warning <- paste(c(points_warning, warning_msg), collapse = "; ")
  if (combined_warning == "") combined_warning <- NULL

  list(
    species = sci_name,
    method = method,
    n_points = nrow(pts_unique),
    points = points_out,
    polygons = poly_geojson,
    warning = combined_warning,
    text_defaults = text_defaults
  )
}







###############################
# Helper: safe GeoJSON → sf → GeoJSON
###############################
# st_make_valid / smooth / split may return GEOMETRYCOLLECTIONs or degenerate
# lines/points: keep only polygonal parts and explode them into single POLYGON
# features (a direct st_cast to POLYGON fails or drops multipolygon parts)
sf_polygons_only <- function(sf_obj) {
  sf_obj <- sf::st_make_valid(sf_obj)
  if (any(sf::st_geometry_type(sf_obj) == "GEOMETRYCOLLECTION")) {
    sf_obj <- sf::st_collection_extract(sf_obj, "POLYGON", warn = FALSE)
  }
  sf_obj <- sf_obj[grepl("POLYGON", sf::st_geometry_type(sf_obj)) & !sf::st_is_empty(sf_obj), ]
  if (nrow(sf_obj) == 0) stop("No polygon left after the operation")
  sf::st_cast(sf::st_cast(sf_obj, "MULTIPOLYGON", warn = FALSE), "POLYGON", warn = FALSE)
}

sf_to_geojson_safe <- function(sf_obj) {
  if (!inherits(sf_obj, "sf")) stop("Object must be of class 'sf'")
  sf_obj <- sf::st_make_valid(sf_obj)
  sf_obj <- sf::st_transform(sf_obj, 4326)
  geojson_text <- geojsonsf::sf_geojson(sf_obj)
  geojson_list <- jsonlite::fromJSON(geojson_text, simplifyVector = FALSE)
  geojson_list
}

###############################
# POST /simplify
#* Simplify polygon geometries
#* @post /simplify
function(req, res) {
  body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)

  poly_geojson <- body$polygon
  tol <- as.numeric(body$tolerance %||% 0)

  if (is.null(poly_geojson)) {
    res$status <- 400
    return(list(status="error", message="Missing polygon"))
  }

  sf_obj <- sf::st_read(jsonlite::toJSON(poly_geojson, auto_unbox=TRUE), quiet=TRUE)
  sf_obj <- sf_polygons_only(sf_obj)

  simplified_sf <- rmapshaper::ms_simplify(sf_obj, keep = min(1, 3000 / nrow(sf::st_coordinates(sf_obj))), keep_shapes = TRUE)
  simplified_sf <- sf::st_make_valid(simplified_sf)
  simplified_sf <- simplified_sf[grepl("POLYGON", sf::st_geometry_type(simplified_sf)), ]

  list(status="ok", polygon=sf_to_geojson_safe(simplified_sf))
}

###############################
# POST /smooth
#* Smooth polygon geometries
#* @post /smooth
function(req, res) {
  body <- tryCatch(
    jsonlite::fromJSON(req$postBody, simplifyVector = FALSE),
    error = function(e) {
      res$status <- 400
      return(list(status="error", message=paste("Invalid JSON body:", e$message)))
    }
  )
  
  poly_geojson <- body$polygon
  smooth_val <- as.numeric(body$smoothness %||% 0)

  if (is.null(poly_geojson)) {
    res$status <- 400
    return(list(status = "error", message = "Missing polygon"))
  }
  # Convert GeoJSON to sf
  sf_obj <- sf::st_read(jsonlite::toJSON(poly_geojson, auto_unbox = TRUE), quiet = TRUE)
  # Ensure valid geometries
  sf_obj <- sf_polygons_only(sf_obj)
  # PROJECT TO METERS (EPSG:3857) before smooth
  sf_obj_m <- sf::st_transform(sf_obj, 3857)
  # Compute smooth parameter as in Shiny
  smooth_par <- exp(smooth_val / 20) - 1
  # Apply smooth (smoothr::smooth requires smoothness > 0; a value of 0 means "no smoothing")
  if (smooth_par > 0) {
    smoothed_sf_m <- smoothr::smooth(sf_obj_m, method = "ksmooth", smoothness = smooth_par, max_distance = 10000)
  } else {
    smoothed_sf_m <- sf_obj_m
  }
  # Keep polygons only & validate
  smoothed_sf_m <- sf_polygons_only(smoothed_sf_m)
  # Transform back to WGS84
  smoothed_sf <- sf::st_transform(smoothed_sf_m, 4326)
  # Convert to GeoJSON safely
  geojson_out <- sf_to_geojson_safe(smoothed_sf)
  list(status = "ok", polygon = geojson_out)
}

#* Split polygon by line
#* @post /split
function(req, res) {
  
  library(sf)
  library(lwgeom)
  library(geojsonsf)
  library(jsonlite)
  
  # 🔥 usa direttamente il body raw
  body_txt <- req$postBody
  
  # parse SOLO per estrarre i pezzi
  body <- jsonlite::fromJSON(body_txt, simplifyVector = FALSE)
  
  if (is.null(body$polygon) || is.null(body$line)) {
    res$status <- 400
    return(list(status = "error", message = "Missing polygon or line"))
  }
  
  # 🔥 IMPORTANTISSIMO: NON riconvertire con toJSON()
  poly_txt <- jsonlite::toJSON(body$polygon, auto_unbox = TRUE)
  line_txt <- jsonlite::toJSON(body$line, auto_unbox = TRUE)
  
  # 🔥 Converti direttamente
  poly_sf <- geojsonsf::geojson_sf(poly_txt)
  line_sf <- geojsonsf::geojson_sf(line_txt)
  
  # Fix geometrie
  poly_sf <- sf_polygons_only(poly_sf)
  line_sf <- st_cast(line_sf, "LINESTRING")
  
  # 🔥 SPLIT (only polygons intersecting the line, as in Shiny_EditPoly.R)
  hits <- lengths(st_intersects(poly_sf, line_sf)) > 0
  if (!any(hits)) {
    res$status <- 400
    return(list(status = "error", message = "No split produced: the line does not intersect any polygon"))
  }
  split <- lwgeom::st_split(poly_sf[hits, ], st_union(line_sf))
  split_sf <- sf_polygons_only(split)
  split_sf <- rbind(poly_sf[!hits, ], split_sf[, names(poly_sf)])
  
  split_sf <- st_transform(split_sf, 4326)
  
  return(list(
    polygons = sf_to_geojson_clean(split_sf)
  ))
}

###############################
# POST /discard
#* Discard polygon changes and restore original
#* @post /discard
function(req, res) {
  body <- tryCatch(
    jsonlite::fromJSON(req$postBody, simplifyVector = FALSE),
    error = function(e) {
      res$status <- 400
      return(list(status="error", message=paste("Invalid JSON body:", e$message)))
    }
  )
  original_poly <- body$original_polygon
  if (is.null(original_poly)) {
    res$status <- 400
    return(list(status="error", message="Missing original polygon"))
  }
  list(status="ok", polygon=original_poly)
}

###################################
# POST /savePolygon
###################################
#* Save edited polygon
#* @post /savePolygon
function(req, res) {

  is_error <- function(x) {
    isTRUE(is.list(x) && !is.null(x$status) && x$status == "error")
  }

  # ----------------------------
  # Parse JSON
  # ----------------------------
  body <- tryCatch(
    jsonlite::fromJSON(req$postBody, simplifyVector = FALSE),
    error = function(e) {
      res$status <- 400
      return(list(status="error", message=paste("Invalid JSON:", e$message)))
    }
  )

  if (is_error(body)) return(body)

  poly_geojson <- body$polygon
  sci_name     <- body$sci_name
  username         <- body$username

  if (is.null(poly_geojson) || is.null(sci_name) || is.null(username)) {
    res$status <- 400
    return(list(status="error", message="Missing polygon, sci_name or username"))
  }

  # ----------------------------
  # GeoJSON → sf
  # ----------------------------
  sf_obj <- tryCatch({
    sf::st_read(jsonlite::toJSON(poly_geojson, auto_unbox = TRUE), quiet = TRUE)
  }, error = function(e) {
    res$status <- 400
    return(list(status="error", message=paste("Invalid GeoJSON:", e$message)))
  })

  if (is_error(sf_obj)) return(sf_obj)

  # ----------------------------
  # Validazione geometria
  # ----------------------------
  sf_obj <- tryCatch({
    sf_polygons_only(sf_obj)
  }, error = function(e) {
    res$status <- 400
    return(list(status="error", message=paste("Invalid geometry:", e$message)))
  })

  if (is_error(sf_obj)) return(sf_obj)

  sf_obj <- sf::st_transform(sf_obj, 4326)

  # ----------------------------
  # Default attributes for newly drawn polygons
  # (leaflet-geoman features arrive with empty properties; without these
  # defaults they get filtered out of distribution-attributes by presence/
  # seasonal/origin, mirroring Victor's Shiny_EditPoly.R behaviour)
  # ----------------------------
  if (!"presence" %in% names(sf_obj)) sf_obj$presence <- NA
  if (!"origin"   %in% names(sf_obj)) sf_obj$origin   <- NA
  if (!"seasonal" %in% names(sf_obj)) sf_obj$seasonal <- NA
  if (!"binomial" %in% names(sf_obj)) sf_obj$binomial <- NA

  sf_obj$presence[is.na(sf_obj$presence)] <- 1
  sf_obj$origin[is.na(sf_obj$origin)]     <- 1
  sf_obj$seasonal[is.na(sf_obj$seasonal)] <- 1
  sf_obj$binomial[is.na(sf_obj$binomial)] <- sci_name

  # ----------------------------
  # Text attributes (Fill text attributes panel)
  # ----------------------------
  text_attr <- body$text_attributes
  if (!is.null(text_attr)) {
    if (!is.null(text_attr$source))     sf_obj$source     <- text_attr$source
    if (!is.null(text_attr$yrcompiled)) sf_obj$yrcompiled <- text_attr$yrcompiled
    if (!is.null(text_attr$citation))   sf_obj$citation   <- text_attr$citation
    if (!is.null(text_attr$compiler))   sf_obj$compiler   <- text_attr$compiler
    if (!is.null(text_attr$island))     sf_obj$island     <- text_attr$island
    if (!is.null(text_attr$dist_comm))  sf_obj$dist_comm  <- text_attr$dist_comm
    if (!is.null(text_attr$data_sens))  sf_obj$data_sens  <- as.numeric(text_attr$data_sens)
  }

  # ----------------------------
  # LOAD STORAGE
  # ----------------------------
  Storage_SP <- tryCatch({
    sRL_StoreRead(sci_name, username, MANDAT = 1)
  }, error = function(e) {
    res$status <- 500
    return(list(status="error", message=paste("Storage read failed:", e$message)))
  })

  if (is_error(Storage_SP)) return(Storage_SP)

  if (is.null(Storage_SP)) {
    res$status <- 404
    return(list(status="error", message="Storage not found"))
  }

  # ----------------------------
  # UPDATE POLYGON
  # ----------------------------
  Storage_SP$distSP_saved <- sf_obj

  # ----------------------------
  # SAVE (QUI STA LA FIX)
  # ----------------------------
  save_result <- tryCatch({

    sRL_StoreSave(sci_name, username, Storage_SP)

    list(status="ok", message="Polygon saved successfully")

  }, error = function(e) {
    res$status <- 500
    list(status="error", message=paste("Save failed:", e$message))
  })

  return(save_result)
}


################################
### HYDROBASINS HELPERS #######
################################
# hybas_id / next_down / next_sink are 10-digit integers that lose precision when
# serialized as JS doubles. Force them to character on an sf object before geojson.
sRL_hydroIdsToChar <- function(sf_obj) {
  for (col in c("hybas_id", "next_down", "next_sink")) {
    if (col %in% names(sf_obj)) {
      sf_obj[[col]] <- format(sf_obj[[col]], scientific = FALSE, trim = TRUE)
    }
  }
  sf_obj
}

# Fast candidate pre-filter: numeric bbox-overlap test against the cached hydro_bbox
# matrix (no GEOS/s2), returning a small subset of hydro_raw whose bounding boxes
# intersect the query distribution's buffered bbox. PrepareHydro then runs its exact
# st_filter on this reduced set instead of all 234k basins.
sRL_hydroPrefilter <- function(distSP, buffer_m = 50000) {
  if (is.null(hydro_bbox)) return(hydro_raw)  # fall back to full layer if cache missing
  bb <- sf::st_bbox(sf::st_buffer(sf::st_as_sfc(sf::st_bbox(distSP)), buffer_m))
  keep <- hydro_bbox[, "xmin"] <= bb["xmax"] & hydro_bbox[, "xmax"] >= bb["xmin"] &
          hydro_bbox[, "ymin"] <= bb["ymax"] & hydro_bbox[, "ymax"] >= bb["ymin"]
  hydro_raw[keep, ]
}


######################################
### HYDROBASINS - LOAD DATA (A) #####
######################################
#* Load hydrobasins for a species (native Edit-Poly hydrobasins variant)
#* @post /species/<sci_name>/hydrobasins/load
function(req, res, sci_name, username = req$argsQuery$username) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # sci_name arrives URL-encoded from the path (e.g. "Saara%20loricata"); normalize it
  sci_name <- gsub("_", " ", utils::URLdecode(sci_name))

  if (missing(username) || is.null(username) || username == "") {
    res$status <- 400
    return(list(error = "Missing username parameter"))
  }

  # Parse body { method, buffer_km }
  body <- tryCatch(
    jsonlite::fromJSON(req$postBody, simplifyVector = TRUE),
    error = function(e) list()
  )
  method    <- tolower(body$method %||% "hydro8")
  buffer_km <- as.numeric(body$buffer_km %||% 50)

  # hydroMCP produces a level-8 hydrobasins distribution at creation time (the MCP is only
  # used to select which level-8 basins to include), so in the editor it behaves like hydro8.
  if (method == "hydromcp") method <- "hydro8"

  if (!method %in% c("hydro8", "hydro10", "hydro12")) {
    res$status <- 400
    return(list(error = paste("Unsupported method:", method)))
  }

  log_info("HYDRO LOAD species={sci_name} username={username} method={method} buffer={buffer_km}")


  # Ensure hydro layers are available
  load_hydro()

  # Load storage
  Storage_SP <- tryCatch(sRL_StoreRead(sci_name, username, MANDAT = 1),
                         error = function(e) { log_error("Storage read failed: {e$message}"); NULL })
  if (is.null(Storage_SP)) { res$status <- 404; return(list(error = "Species storage not found")) }

  Output <- Storage_SP$Output
  getOut <- function(par) { v <- Output$Value[Output$Parameter == par]; if (length(v) == 0) NA else v[1] }

  # Decide source distribution + SRC_created + HydroLev
  dist_source_is_created <- (isTRUE(getOut("Distribution_Source") == "Created")) &&
    (TRUE %in% grepl("hydro", Output$Value))

  if (buffer_km >= 100 && !is.null(Storage_SP$distSP_saved_tempoHydro)) {
    # Expand: reload from the last saved hydro state, using a wider buffer
    dist_src    <- Storage_SP$distSP_saved_tempoHydro
    SRC_created  <- "yes"
    HydroLev     <- if (!is.na(getOut("Mapping_Start"))) getOut("Mapping_Start") else method
  } else if (dist_source_is_created) {
    dist_src    <- Storage_SP$distSP3_BeforeCrop
    SRC_created  <- "yes"
    HydroLev     <- if (!is.na(getOut("Mapping_Start"))) getOut("Mapping_Start") else method
  } else {
    dist_src    <- Storage_SP$distSP_saved
    SRC_created  <- "no"
    HydroLev     <- method
  }

  # Guard BEFORE any transform: the source distribution may be missing
  if (is.null(dist_src) || !inherits(dist_src, "sf") || nrow(dist_src) == 0) {
    res$status <- 404
    return(list(error = "No source distribution available to build hydrobasins. Create the distribution (Step 3) before editing it as hydrobasins."))
  }

  dist_loaded0 <- sf::st_transform(dist_src, sf::st_crs(hydro_raw))


  # Fast numeric bbox pre-filter (60km margin ≥ the 50km buffer used inside PrepareHydro),
  # so PrepareHydro's exact st_filter runs on a few hundred basins instead of 234k.
  hydro_sub <- sRL_hydroPrefilter(dist_loaded0, buffer_m = 60000)

  # Prepare hydrobasins (heavy geo, reuses the Shiny function; skip the Shiny popup — rebuilt on the client)
  hydro_ready <- tryCatch(
    sRLPolyg_PrepareHydro(dist_loaded0, hydro_sub, HydroLev, SRC_created = SRC_created, make_popup = FALSE),
    error = function(e) { log_error("PrepareHydro failed: {e$message}"); NULL }
  )
  if (is.null(hydro_ready) || length(hydro_ready) == 0 ||
      is.null(hydro_ready$hydroSP) || nrow(hydro_ready$hydroSP) == 0) {
    res$status <- 422
    return(list(error = "The distribution does not overlap with hydrobasins"))
  }

  hydroSP    <- hydro_ready$hydroSP      # simplified, WGS84, for display
  hydroSP_HQ <- hydro_ready$hydroSP_HQ   # high quality, for save

  # Stash HQ geometry in storage so /save can swap to it without shipping it to the client
  Storage_SP$hydroSP_HQ <- hydroSP_HQ
  tryCatch(sRL_StoreSave(sci_name, username, Storage_SP),
           error = function(e) log_warn("Could not stash hydroSP_HQ: {e$message}"))

  # Level-3 overlay within 10km of the hydrobasins (decorative background layer, hidden by
  # default) — simplify aggressively: it dominates the payload otherwise (~750KB → ~40KB).
  hydro3_out <- tryCatch({
    h3 <- sf::st_filter(hydro3, sf::st_buffer(sf::st_transform(hydroSP, sf::st_crs(hydro3)), 10))
    if (nrow(h3) > 0) {
      h3 <- tryCatch(rmapshaper::ms_simplify(h3, keep = 0.02, keep_shapes = TRUE), error = function(e) h3)
      sf_to_geojson_clean(sf::st_transform(h3, 4326))
    } else NULL
  }, error = function(e) { log_warn("hydro3 overlay failed: {e$message}"); NULL })

  # Occurrence points (reuse the same sources as /mapdata)
  points_out <- list()
  tryCatch({
    if (!is.null(Storage_SP$dat_proj_saved)) {
      pts <- sf::st_transform(Storage_SP$dat_proj_saved, 4326)
      cc <- sf::st_coordinates(pts); pts$decimalLongitude <- cc[,1]; pts$decimalLatitude <- cc[,2]
      points_out <- normalize_json_table(pts)
    } else if (!is.null(Storage_SP$data_comparison)) {
      pts <- sf::st_as_sf(Storage_SP$data_comparison,
                          coords = c("decimalLongitude", "decimalLatitude"),
                          crs = 4326, remove = FALSE)
      points_out <- normalize_json_table(pts)
    }
  }, error = function(e) log_warn("Points extraction failed: {e$message}"))

  # Text defaults for the attribute panel (port of Shiny_EditPoly.R:325-346)
  text_defaults <- tryCatch({
    la <- sf::st_drop_geometry(hydroSP)
    names(la) <- tolower(names(la))
    firstNonNull <- function(col, default) {
      if (col %in% names(la)) { v <- la[[col]][!is.na(la[[col]])]; if (length(v) > 0) return(as.character(v[1])) }
      default
    }
    # Distribution comment autofill (port of Shiny_EditPoly.R:338-345)
    dist_source <- getOut("Distribution_Source")
    existing_comm <- if ("dist_comm" %in% names(la)) { v <- la[["dist_comm"]][!is.na(la[["dist_comm"]])]; if (length(v) > 0) as.character(v[1]) else NA } else NA
    if (isTRUE(dist_source == "Created")) {
      dist_comm_val <- if (is.na(existing_comm)) "" else existing_comm
    } else {
      prefix <- if (is.na(existing_comm)) {
        paste0("The distribution was taken from the sRedList platform (", dist_source, ") and was ")
      } else {
        paste0(existing_comm, " It was ")
      }
      dist_comm_val <- paste0(prefix, "manually edited on the sRedList platform on the ", Sys.Date(), ".")
    }
    dist_comm_val <- substr(dist_comm_val, 1, 254)

    list(
      source     = firstNonNull("source", "sRedList platform"),
      yrcompiled = format(Sys.time(), "%Y"),
      citation   = firstNonNull("citation", "IUCN (International Union for Conservation of Nature)"),
      compiler   = firstNonNull("compiler", tryCatch(sRL_userformatted(username), error = function(e) username)),
      island     = firstNonNull("island", ""),
      data_sens  = 0,
      sens_comm  = firstNonNull("sens_comm", ""),
      dist_comm  = dist_comm_val
    )
  }, error = function(e) list())

  # Drop the R-generated Shiny popup (rebuilt on the client) + internal helper cols
  drop_cols <- intersect(c("Popup", "cols", "ar", "N", "Grid_cells", "dist_comm", "binomial"), names(hydroSP))
  if (length(drop_cols) > 0) hydroSP <- hydroSP[, setdiff(names(hydroSP), drop_cols)]


  # Simplify geometries for DISPLAY only (the HQ geometry is stashed separately for /save),
  # shrinking the payload substantially. Attributes/hybas_id are preserved.
  hydroSP <- tryCatch(rmapshaper::ms_simplify(hydroSP, keep = 0.1, keep_shapes = TRUE),
                      error = function(e) hydroSP)

  # Serialize hydrobasins FeatureCollection (ids as strings)
  hydro_geojson <- sf_to_geojson_clean(sRL_hydroIdsToChar(sf::st_transform(hydroSP, 4326)))

  # Ensure empty (NA) presence/origin/seasonal serialize as JSON null, not {}
  hydro_geojson$features <- lapply(hydro_geojson$features, function(feat) {
    for (k in c("presence", "origin", "seasonal")) {
      v <- feat$properties[[k]]
      if (is.null(v) || length(v) == 0 || (is.list(v) && length(v) == 0)) {
        feat$properties[[k]] <- NULL
      }
    }
    feat
  })

  list(
    species     = sci_name,
    method      = method,
    buffer_km   = buffer_km,
    src_created = SRC_created,
    hydrobasins = hydro_geojson,
    hydro3      = hydro3_out,
    points      = points_out,
    text_defaults = text_defaults
  )
}


######################################
### HYDROBASINS - SAVE (C) ##########
######################################
#* Save the edited hydrobasins distribution
#* @post /species/<sci_name>/hydrobasins/save
function(req, res, sci_name, username = req$argsQuery$username) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # sci_name arrives URL-encoded from the path (e.g. "Saara%20loricata"); normalize it
  sci_name <- gsub("_", " ", utils::URLdecode(sci_name))

  if (missing(username) || is.null(username) || username == "") {
    res$status <- 400
    return(list(status = "error", message = "Missing username parameter"))
  }

  body <- tryCatch(jsonlite::fromJSON(req$postBody, simplifyVector = FALSE),
                   error = function(e) NULL)
  if (is.null(body)) { res$status <- 400; return(list(status = "error", message = "Invalid JSON body")) }

  attributes_all <- body$attributes_all %||% list()
  text_fields    <- body$text_fields %||% list()

  if (length(attributes_all) == 0) {
    res$status <- 400
    return(list(status = "error", message = "No hydrobasins attributes provided"))
  }

  dist_comm <- as.character(text_fields$dist_comm %||% "")
  if (nchar(dist_comm) > 254) {
    res$status <- 400
    return(list(status = "error", message = paste0("Distribution comment too long (", nchar(dist_comm), " chars, max 254)")))
  }

  log_info("HYDRO SAVE species={sci_name} username={username} n_attr={length(attributes_all)}")

  load_hydro()

  Storage_SP <- tryCatch(sRL_StoreRead(sci_name, username, MANDAT = 1),
                         error = function(e) NULL)
  if (is.null(Storage_SP)) { res$status <- 404; return(list(status = "error", message = "Species storage not found")) }

  hydroSP_HQ <- Storage_SP$hydroSP_HQ
  if (is.null(hydroSP_HQ) || !inherits(hydroSP_HQ, "sf")) {
    res$status <- 409
    return(list(status = "error", message = "No HQ hydrobasins in storage; reload the hydrobasins before saving"))
  }

  # Build attribute frame keyed by hybas_id (as character)
  attr_df <- do.call(rbind, lapply(attributes_all, function(a) {
    data.frame(
      hybas_id = as.character(a$hybas_id),
      presence = if (is.null(a$presence)) NA_real_ else as.numeric(a$presence),
      origin   = if (is.null(a$origin))   NA_real_ else as.numeric(a$origin),
      seasonal = if (is.null(a$seasonal)) NA_real_ else as.numeric(a$seasonal),
      stringsAsFactors = FALSE
    )
  }))

  # Reconstruct dist_tosave from HQ geometry (HQ swap), join attributes by hybas_id
  dist_tosave <- hydroSP_HQ
  key <- as.character(format(dist_tosave$hybas_id, scientific = FALSE, trim = TRUE))
  dist_tosave$presence <- attr_df$presence[match(key, attr_df$hybas_id)]
  dist_tosave$origin   <- attr_df$origin[match(key, attr_df$hybas_id)]
  dist_tosave$seasonal <- attr_df$seasonal[match(key, attr_df$hybas_id)]

  # Text + identity attributes (port of Shiny_EditPoly.R:1068-1078)
  dist_tosave$binomial   <- sci_name
  dist_tosave$id_no      <- tryCatch(sRL_CalcIdno(sci_name), error = function(e) NA)
  dist_tosave$source     <- as.character(text_fields$source %||% "sRedList platform")
  dist_tosave$yrcompiled <- as.character(text_fields$yrcompiled %||% format(Sys.time(), "%Y"))
  dist_tosave$citation   <- as.character(text_fields$citation %||% "")
  dist_tosave$compiler   <- as.character(text_fields$compiler %||% "")
  dist_tosave$data_sens  <- ifelse(isTRUE(text_fields$data_sens) || identical(text_fields$data_sens, 1), 1, 0)
  dist_tosave$sens_comm  <- as.character(text_fields$sens_comm %||% "")
  dist_tosave$island     <- as.character(text_fields$island %||% "")
  dist_tosave$dist_comm  <- dist_comm

  dist_tosave <- sf::st_transform(sf::st_make_valid(dist_tosave), 4326)

  # Keep the full (incl. NA) set for discard
  Storage_SP$distSP_saved_tempoHydro <- dist_tosave

  # Drop empty basins (presence NA) for the real save
  dist_final <- dist_tosave[!is.na(dist_tosave$presence), ]
  if (nrow(dist_final) == 0) {
    res$status <- 400
    return(list(status = "error", message = "Distribution is empty (no occupied hydrobasins)"))
  }

  # Assign attributes to occurrence points by intersection
  tryCatch({
    if (!is.null(Storage_SP$dat_proj_saved)) {
      pts <- sf::st_transform(Storage_SP$dat_proj_saved, 4326)
      pts_inter <- sf::st_join(pts, dist_final, join = sf::st_intersects)
      if ("gbifID" %in% names(Storage_SP$dat_proj_saved)) {
        m <- match(Storage_SP$dat_proj_saved$gbifID, pts_inter$gbifID)
        Storage_SP$dat_proj_saved$presence <- pts_inter$presence[m]
        Storage_SP$dat_proj_saved$origin   <- pts_inter$origin[m]
        Storage_SP$dat_proj_saved$seasonal <- pts_inter$seasonal[m]
      }
    }
  }, error = function(e) log_warn("Point intersection failed: {e$message}"))

  # Regenerate report plot (port of Shiny_EditPoly.R:1112-1118)
  tryCatch({
    dist_final$cols <- sRL_ColourDistrib(dist_final)$cols
    plot_dist <- ggplot2::ggplot()
    if (!is.null(Storage_SP$CountrySP_saved)) {
      plot_dist <- plot_dist + ggplot2::geom_sf(data = Storage_SP$CountrySP_saved, fill = "gray96", col = "gray50")
    }
    plot_dist <- plot_dist +
      ggplot2::geom_sf(data = dist_final, fill = dist_final$cols) +
      ggplot2::theme_void() + ggplot2::ggtitle("")
    plot_path <- paste0("resources/AOH_stored/", sub(" ", "_", sci_name), "_", sRL_userdecode(username),
                        "/Plots/plot_manually_edited.png")
    dir.create(dirname(plot_path), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(plot_path, plot_dist, width = 10, height = 8)
  }, error = function(e) log_warn("Report plot generation failed: {e$message}"))

  # Persist
  Storage_SP$distSP_saved <- dist_final
  Storage_SP$Output$Value[Storage_SP$Output$Parameter == "Gbif_EditPoly"] <- "yes"
  save_ok <- tryCatch({
    sRL_StoreSave(sci_name, username, Storage_SP)
    tryCatch(sRL_saveMapDistribution(sci_name, Storage_SP), error = function(e) log_warn("saveMapDistribution failed: {e$message}"))
    TRUE
  }, error = function(e) { log_error("Hydro save failed: {e$message}"); res$status <- 500; FALSE })

  if (!isTRUE(save_ok)) return(list(status = "error", message = "Save failed"))

  area_km2 <- tryCatch(as.numeric(sum(sf::st_area(sf::st_transform(dist_final, CRSMOLL)))) / 1e6,
                       error = function(e) NA)

  list(
    status      = "ok",
    message     = "Hydrobasins distribution saved successfully",
    area_km2    = area_km2,
    hydrobasins = sf_to_geojson_clean(sRL_hydroIdsToChar(dist_final))
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
