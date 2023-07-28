#' synanthropy: A Research Compendium
#' 
#' @description 
#' Adaptation of: https://github.com/lomorel/SynAnthrop
#' 
#' @author Nicolas Casajus \email{nicolas.casajus@fondationbiodiversite.fr}
#' 
#' @date 2023/07/28



## Install Dependencies (listed in DESCRIPTION) ----

remotes::install_deps(upgrade = "never")


## Load Project Addins (R Functions and Packages) ----

pkgload::load_all(here::here())


## Import France raster (CartNat) ----

fra <- raster::raster(here::here("data", 
                                 "Layer4_FINAL.tif"))


## Import species abundances (STOC) ----

spp <- read.csv2(here::here("data", 
                            "data_FrenchBBS_square_Barnagaud_allSp_2001_2023.csv"))


## Select columns ----

spp <- spp[ , c("scientific_name", "annee", "abondance_brut", 
                "longitude_grid_wgs84", "latitude_grid_wgs84")]

colnames(spp) <- c("Species", "Year", "Abundance", "X", "Y")


## Subset species ----

species <- c("Columba livia", "Phoenicurus ochruros", "Apus apus", 
             "Prunella modularis", "Turdus merula", "Parus major", 
             "Carduelis carduelis", "Dryocopus martius", "Lullula arborea", 
             "Phylloscopus bonelli")

new_spp <- spp[which(spp$"Species" %in% species), ]


## Remove missing coords ----

new_spp <- new_spp[which(!is.na(new_spp$"X")), ]


## Transform CRS ----

dat <- sf::st_as_sf(new_spp, coords = c("X", "Y"), crs = "epsg:4326")
dat <- sf::st_transform(dat, raster::crs(fra))

dat_data  <- sf::st_drop_geometry(dat)
dat_coord <- sf::st_coordinates(dat)
dat <- cbind(dat_data, dat_coord)


## Compute SSI ----

ssi_results <- ssi(r = fra, x = dat, resolution = c(100, 250, 500), sim = 200, 
                   threshold = 30)


## Score distribution within the studied taxa ----

sp_ssi <- ssi_results[[1]]
sp_ssi$Index <- as.integer(sp_ssi$"Index")

mean_ssi_by_resolution <- data.frame(sp_ssi %>%
                                       dplyr::group_by(Species, Resolution) %>%
                                       dplyr::summarise(Index = mean(Index), 
                                                        n = dplyr::n()))

effsize_res <- ssi_results[[2]]


# then scale by scale
sub_effsize_res <- subset(effsize_res, Resolution == "250")
mean_ssi_by_resolution <- subset(mean_ssi_by_resolution, Resolution == "250")

sub_effsize_res <- merge(sub_effsize_res, mean_ssi_by_resolution,
                         by = "Species")
sub_effsize_res$"Index" <- as.factor(sub_effsize_res$"Index")

ggplot(sub_effsize_res, aes(x = reorder(Species, -effsize), 
                            y = -effsize, fill = Index)) +
  geom_hline(yintercept = 0.0, color = "darkgrey", size=0.8, 
             linetype="dashed") +
  geom_boxplot() + 
  coord_flip() +
  scale_fill_brewer(name="Score", palette = "RdYlGn") +
  ylab("Effect size") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = c(0.9, 0.2)) +
  theme_bw()


## Impacts of resolution ----

mean_index_by_reso <- data.frame(
  sp_ssi %>%
    dplyr::group_by(Species, Resolution) %>%
    dplyr::summarise(Index = mean(Index), n = dplyr::n()))

mean_index_by_reso$"Resolution" <- as.factor(mean_index_by_reso$"Resolution")

CGPfunctions::newggslopegraph(
  dataframe = mean_index_by_reso,
  Resolution,
  Index,
  Grouping = Species,
  Title = "Synanthropy scores for birds species in France",
  SubTitle = NULL,
  Caption = NULL)


## Maps ----

sub_distri <- subset(ssi_results[[3]], Resolution == "250")
sub_distri_obs <- subset(sub_distri, variable  == "Obs")
distri <- rbind(sub_distri_obs, 
                dplyr::sample_n(subset(sub_distri, variable  == "Null" ), 
                                nrow(sub_distri_obs)))
distri <- sf_transform_xy(distri, raster::crs(fra), 2154)

ggplot() + 
  geom_point(data = subset(distri, variable == "Null"), aes(x = x, y = y, 
                                                            color = variable), 
             size = 1.5, colour = "#ff6600",alpha = 1/2) +
  geom_point(data = subset(distri, variable == "Obs"), aes(x = x, y = y, 
                                                           color = variable), 
             size = 1.5, colour = "#660099",alpha = 1/2) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
  facet_wrap(. ~ Species) +
  theme_void() +
  theme(legend.position="none")
