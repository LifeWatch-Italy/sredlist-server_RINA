#!/usr/bin/env Rscript

library(terra)
library(fs)
library(DBI)
library(RSQLite)
library(png)

rasters_path <- "resources/resources_Shiny_DD/1.GFC_final/"
tiles_path   <- "resources/resources_Shiny_DD/1.GFC_tiles/"
dir_create(tiles_path)

rasters <- dir_ls(rasters_path, glob="*_RANGEsmall.tif")

pal_fun <- colorRampPalette(c("#d01c8b","#a1d76a"))
pal <- pal_fun(100)

for(ras_path in rasters){

  sp_name <- gsub("_RANGEsmall.tif","", path_file(ras_path))
  mbtiles_file <- file.path(tiles_path, paste0(sp_name,".mbtiles"))

  if(file_exists(mbtiles_file)){
    cat("MBTiles already created:", sp_name,"\n")
    next
  }

  cat("\nProcessing:", sp_name,"\n")

  r <- rast(ras_path)

  vals <- values(r)
  vals[vals==0] <- NA
  values(r) <- vals

  vals_scaled <- round(vals)
  vals_scaled[vals_scaled<1] <- 1
  vals_scaled[vals_scaled>100] <- 100

  cols <- rep(NA,length(vals_scaled))
  idx <- !is.na(vals_scaled)
  cols[idx] <- pal[vals_scaled[idx]]

  rgb_vals <- matrix(NA,nrow=length(cols),ncol=3)
  rgb_vals[idx,] <- t(col2rgb(cols[idx]))

  r_rgb <- rast(r,nlyrs=3)
  values(r_rgb) <- rgb_vals

  alpha <- rast(r)
  values(alpha) <- ifelse(is.na(vals),0,255)

  r_rgba <- c(r_rgb,alpha)

  tmp_rgba <- tempfile(fileext=".tif")
  writeRaster(r_rgba,tmp_rgba,datatype="INT1U",overwrite=TRUE)

  tmp_3857 <- tempfile(fileext=".tif")

  system2("gdalwarp",c(
    "-t_srs","EPSG:3857",
    "-r","near",
    "-dstalpha",
    tmp_rgba,
    tmp_3857
  ))

  tiles_dir <- tempfile()
  dir_create(tiles_dir)

  system2("gdal2tiles.py", c(
    "--xyz",
    "--zoom=5-18",
    "--processes=4",
    tmp_3857,
    tiles_dir
  ))

  # -------------------
  # create MBTiles
  # -------------------

  con <- dbConnect(SQLite(),mbtiles_file)

  dbExecute(con,"
  CREATE TABLE tiles (
    zoom_level INTEGER,
    tile_column INTEGER,
    tile_row INTEGER,
    tile_data BLOB
  )")

  dbExecute(con,"
  CREATE TABLE metadata (
    name TEXT,
    value TEXT
  )")

  png_files <- dir_ls(tiles_dir, recurse=TRUE, glob="*.png")

  stmt <- dbSendStatement(con,
  "INSERT INTO tiles (zoom_level,tile_column,tile_row,tile_data)
  VALUES (?,?,?,?)"
  )

  for(png in png_files){

    parts <- strsplit(png,"/")[[1]]

    z <- as.integer(parts[length(parts)-2])
    x <- as.integer(parts[length(parts)-1])
    y <- as.integer(gsub(".png","",parts[length(parts)]))

    # XYZ -> TMS
    y <- (2^z - 1) - y

    raw_tile <- readBin(png,"raw",file.info(png)$size)

    dbBind(stmt, list(z, x, y, I(list(raw_tile))))
  }

  dbClearResult(stmt)

  dbExecute(con,"
  INSERT INTO metadata VALUES
  ('name', ?)",params=list(sp_name))

  dbDisconnect(con)

  cat("*** MBTiles created:",mbtiles_file,"\n")

  # cleanup 
  try(unlink(tiles_dir, recursive = TRUE), silent = TRUE)
  try(unlink(tmp_rgba), silent = TRUE)
  try(unlink(tmp_3857), silent = TRUE)

  gc()
}