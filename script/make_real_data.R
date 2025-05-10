# prep

install.packages(c(
  "curl",      # reliable downloader
  "readr", "dplyr", "tidyr",
  "terra",     # raster handling (replaces raster/sp)
  "sf",        # simple‑features for shapefiles
  "ape"        # phylogeny
))
# For MODIS NDVI you also need MODIStsp (CRAN) & GDAL on your system.
#   install.packages("MODIStsp")


library(curl)
library(readr); library(dplyr); library(tidyr)
library(terra); library(sf)
library(ape)

#---------------------------------------------------------------------------#
# 1.  Download AVONET species‑average trait table (∼9 MB CSV, CC‑BY‑4.0)    #
#---------------------------------------------------------------------------#
avonet_url <- "https://figshare.com/ndownloader/files/34106754"  # Suppl. Data 1
csv_path   <- tempfile(fileext = ".csv")
curl_download(avonet_url, csv_path, mode = "wb")
traits_all <- read_csv(csv_path, show_col_types = FALSE)

#  keep only columns we need & rename to simple names
traits_all <- traits_all %>%
  transmute(
    Species      = BirdTree_species,           # already in BirdTree format
    Diet         = Diet_Primary,               # categorical
    Habitat      = Habitat_broad,              # categorical
    BodyMass     = Body_mass_g,
    WingLen      = Wing_length_mm,
    BeakLength   = Culmen_length_mm
  ) %>%
  filter(!is.na(Species))

#---------------------------------------------------------------------------#
# 2.  Choose ~300 species across the tree for a compact example             #
#---------------------------------------------------------------------------#
set.seed(2025)
species_300 <- traits_all %>%
  group_by(substr(Species, 1, 1)) %>%          # stratify by first letter
  slice_sample(prop = 300 / nrow(traits_all)) %>%
  ungroup() %>% pull(Species)

traits_300 <- traits_all %>% filter(Species %in% species_300)

# Introduce 10 % missing values at random for demonstration
for (col in c("BodyMass", "WingLen", "BeakLength"))
  traits_300[sample.int(nrow(traits_300), 0.1*nrow(traits_300)), col] <- NA

#---------------------------------------------------------------------------#
# 3.  Environmental covariates                                              #
#---------------------------------------------------------------------------#
# 3.1  Get species centroids (AVONET gives range centroids)
centroids <- read_csv(csv_path, show_col_types = FALSE) %>%
  select(BirdTree_species, Centroid_lat = Centroid_Lat, Centroid_lon = Centroid_Long) %>%
  filter(BirdTree_species %in% species_300)

# 3.2  Download WorldClim v2.1 Bio1 (mean annual temp) & Bio12 (precip)
wc_tmp   <- tempfile(); dir.create(wc_tmp)
wc_base  <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_"
curl_download(paste0(wc_base, "bio_1.tif"),  file.path(wc_tmp, "bio1.tif"))
curl_download(paste0(wc_base, "bio_12.tif"), file.path(wc_tmp, "bio12.tif"))
temp_rast <- rast(file.path(wc_tmp, "bio1.tif")) / 10        # °C
prec_rast <- rast(file.path(wc_tmp, "bio12.tif"))            # mm

# 3.3  Elevation (WorldClim elevation layer)
curl_download("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_elev.tif",
              file.path(wc_tmp, "elev.tif"))
elev_rast <- rast(file.path(wc_tmp, "elev.tif"))

# 3.4  NDVI  — use MODIS MOD13C2 annual mean (already aggregated, 0.05°)
# NOTE: this is ∼ 1 MB; we download 2020 composite
ndvi_url <- "https://zenodo.org/record/6496624/files/MODIS_NDVI_2020_mean_0.05.tif?download=1"
curl_download(ndvi_url, file.path(wc_tmp, "ndvi.tif"))
ndvi_rast <- rast(file.path(wc_tmp, "ndvi.tif"))

# 3.5  Extract values for each species centroid
pts <- vect(centroids, geom = c("Centroid_lon", "Centroid_lat"), crs = "EPSG:4326")
env_tbl <- data.frame(
  Species = centroids$BirdTree_species,
  Temp    = terra::extract(temp_rast, pts)[,2],
  Prec    = terra::extract(prec_rast, pts)[,2],
  Elev    = terra::extract(elev_rast, pts)[,2],
  NDVI    = terra::extract(ndvi_rast,  pts)[,2]
)

#---------------------------------------------------------------------------#
# 4.  BirdTree phylogeny subset (Jetz et al. 2012)                           #
#---------------------------------------------------------------------------#
#  BirdTree offers a static tar.gz of 10000 posterior trees; we grab one tree.
tree_url <- "https://birdtree.org/downloaddata/HackettStage2.tgz"
tgz <- tempfile(fileext = ".tgz")
curl_download(tree_url, tgz, mode = "wb")
untar(tgz, files = "HackettStage2/AllBirdsHackett1.nex", exdir = wc_tmp)
full_tree <- read.nexus(file.path(wc_tmp, "HackettStage2/AllBirdsHackett1.nex"))

# Prune to 300 species
tree300 <- drop.tip(full_tree, setdiff(full_tree$tip.label, species_300))

#---------------------------------------------------------------------------#
# 5.  Save compressed .rda files to package data/                            #
#---------------------------------------------------------------------------#
dir.create("../data", showWarnings = FALSE)
save(traits_300, file = "../data/avonet300.rda", compress = "xz")
save(env_tbl,    file = "../data/env300.rda",    compress = "xz")
save(tree300,    file = "../data/tree300.rda",   compress = "xz")
cat("Saved data to ../data/  (sizes):\n")
print(file.info(list.files("../data", full.names = TRUE))[,"size"] / 1024, digits = 2)