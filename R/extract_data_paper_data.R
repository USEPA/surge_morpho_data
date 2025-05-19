library(readr)
library(dplyr)
library(sf)
library(here)
library(terra)
library(tidyr)

################################################
# Get Bathy files to fill in missing site depths
################################################

  bathy_files <- list.files(here("data/bathymetry_raster"), full.names = TRUE)
  surge_poly <- st_read(here("data/all_lakes.gpkg"), layer = "all_lakes") |>
    st_transform(5072)
  surge_pts <- st_read(here("data/all_lakes.gpkg"), layer = "points") |>
    st_transform(5072)

  get_and_crop <- function(bathy_file, res, res_pts, overwrite = TRUE){
    if(overwrite){
      bathy <- rast(bathy_file) |>
        project(terra::crs(res))
      bathy_ext <- ext(bathy)
      bathy_ext_vect <- as.polygons(bathy_ext, crs = crs(bathy))
      bathy_ext_polygons <- st_as_sf(bathy_ext_vect)
      reservoir <- st_intersection(bathy_ext_polygons, res) |>
        st_cast("MULTIPOLYGON")
      points <- st_intersection(bathy_ext_polygons, res_pts)
      bathy_crop <- crop(bathy, reservoir) * -1
      if(reservoir$lake_id == "1006"){
        # It appears the bathymetry for  Buckhorn Reservoir in Kentuky is in Feet
        # Converting to meters.
        bathy_crop <- bathy_crop * 0.3048
      }
      writeRaster(bathy_crop, paste0(here("data/bathymetry_clipped"), "/",
                                     reservoir$lake_id, "_",
                                     basename(bathy_file)),
                  overwrite = TRUE)
      write_sf(reservoir, here("data/bathymetry_reservoirs.gpkg"),
               layer = "reservoirs", append = TRUE)
      write_sf(points, here("data/bathymetry_reservoirs.gpkg"),
               layer = "points", append = TRUE)
    }
  }

  lapply(bathy_files, get_and_crop, res = surge_poly, res_pts = surge_pts,
         overwrite = TRUE)

################################################
# Fill in missing site depths from bathy
################################################
surge_res_morpho <- read_csv(here("data/all_lakes_lakemorpho.csv")) |>
  mutate(lake_id = as.numeric(lake_id))
surge_res_points <- st_read(here("data/all_lakes.gpkg"), "points") |>
  st_transform(5072)
surge_res_points_missing <- surge_res_points |>
  filter(.by = lake_id, all(is.na(site_depth)))
missing_ids <- surge_res_points_missing$lake_id |>
  unique()
bathy_files <- list.files(here("data/bathymetry_clipped"), full.names = TRUE)

fill_in_missing_depth <- function(i){
  surge_res_points_missing_i <- surge_res_points_missing |>
    filter(lake_id == i)
  bathy_file <- bathy_files[grep(i, bathy_files)]
  bathy_missing <- terra::rast(bathy_file) |>
    project(terra::crs(surge_res_points_missing_i))
  # site_id
  bathy_pts <- terra::extract(bathy_missing, surge_res_points_missing_i,
                              bind = TRUE) |>
    as.data.frame() |>
    mutate(lake_id = i)
  names(bathy_pts)[length(names(bathy_pts))] <- "bathy_depth"
  bathy_pts <- mutate(bathy_pts, site_depth = bathy_depth) |>
    select(-bathy_depth)
  surge_res_points_missing_i |>
    select(-site_depth) |>
    left_join(bathy_pts)
}

missing <- lapply(missing_ids, fill_in_missing_depth)
surge_res_points_filled <- bind_rows(surge_res_points, missing) |>
  filter(!is.na(site_depth))

################################################
# Voronoi polygons volumes
################################################

# Removing pts that are outside of the reservoir area that doesn't have bathymetry
# For example Lake Loramie, bathy tif stops short of actual reservori poly but
# some samples occured in the area outside fo the bathy
surge_poly <- st_read(here("data/all_lakes.gpkg"), layer = "all_lakes") |>
  st_transform(5072)
surge_res_points_filled_int <- st_intersection(surge_res_points_filled, surge_poly)
get_volumes <- function(id, reservoir, pts){
  #if(id == "1019"){browser()}
  reservoir <- filter(reservoir, lake_id == id)
  pts <- filter(pts, lake_id == id)
  # Should we use the spsurvey weights? I am gonna say no, because some
  # reservoirs will have other points as well that wouldn't have weights.
  mean_depth = mean(pts$site_depth, na.rm = TRUE)
  area = as.numeric(st_area(reservoir))
  vol_formula <- mean_depth * area

  # Voronoi Poloygons
  # From https://gis.stackexchange.com/questions/362134/i-want-to-create-a-voronoi-diagram-while-retaining-the-data-in-the-data-frame
  st_voronoi_point <- function(points){
    ## points must be POINT geometry
    # check for point geometry and execute if true
    if(!all(st_geometry_type(points) == "POINT")){
      stop("Input not  POINT geometries")
    }
    g = st_combine(st_geometry(points)) # make multipoint
    v = st_voronoi(g)
    v = st_collection_extract(v)
    return(v[unlist(st_intersects(points, v))])
  }

  vor_poly <- st_voronoi_point(pts)
  vor_poly <- st_intersection(vor_poly, st_as_sfc(reservoir))
  #vor_poly <- st_cast(vor_poly, "POLYGON") #Not working Multi-polygons messing things up.
  vor_poly <- st_set_geometry(pts, vor_poly)
  vol_vor <- sum(as.numeric(st_area(vor_poly)*vor_poly$site_depth), na.rm = TRUE)

  data.frame(lake_id  = as.character(id), volume_formula = vol_formula,
             volume_voronoi = vol_vor)

}

sampled_volumes <- lapply(surge_poly$lake_id, get_volumes,
                          reservoir = surge_poly,  pts = surge_res_points_filled_int) |>
  bind_rows() |>
  mutate(lake_id = as.numeric(lake_id)) |>
  pivot_longer(volume_formula:volume_voronoi, names_to = "variables",
               values_to = "values") |>
  mutate(source = case_when(variables == "volume_formula" ~
                              "sampled pts formula",
                            variables == "volume_voronoi" ~
                              "sampled pts voronoi",
                            TRUE ~ NA_character_),
         variables = "volume")

################################################
# max, mean, volume from points
################################################
surge_res_points_filled_nogeo <- surge_res_points_filled
st_geometry(surge_res_points_filled_nogeo) <- NULL
surge_sites_depth <- surge_res_points_filled_nogeo |>
  summarize(
    .by = lake_id,
    mean_depth = mean(site_depth, na.rm = TRUE),
    max_depth = max(site_depth, na.rm = TRUE),
  )

sampled_volume <- filter(sampled_volumes, source == "sampled pts voronoi") |>
  select(lake_id, volume = values)

surge_sites <- left_join(surge_sites_depth, sampled_volume) |>
  mutate(littoral_fraction =
           (1 - (1 - (2.5 / max_depth)^((max_depth / mean_depth) - 1))) * 100) |>
  mutate(source = "surge_sites") |>
  left_join(select(surge_res_morpho, lake_id, lake_name)) |>
  pivot_longer(cols = mean_depth:littoral_fraction, names_to = "variables",
               values_to = "values")

################################################
# Other morpho metrics
################################################
surge_morpho <- surge_res_morpho |>
  select(lake_id, lake_name, surfaceArea:shorelineDevelopment,
         maxWidth:maxLength, -fetch) |>
  mutate(source = "surge_morpho") |>
  pivot_longer(cols = surfaceArea:maxLength, names_to = "variables",
               values_to = "values") |>
  bind_rows(surge_sites) |>
  arrange(lake_id) |>
  mutate(values = round(values, 2))

write_csv(surge_morpho, here("data/surge_morpho.csv"))
sessioninfo::session_info(to_file = "data_paper_session.txt")
