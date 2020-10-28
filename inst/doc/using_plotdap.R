## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = FALSE---------------------------------------------------------
### needed libraries
library(mapdata)
library(plotdap)
library(rerddap)


## ---- eval = FALSE------------------------------------------------------------
#  install.packages("plotdap")

## ---- eval = FALSE------------------------------------------------------------
#  devtools::install_github('ropensci/plotdap')

## ----world, fig.align = 'center', fig.height = 4, fig.width = 5---------------
library(plotdap)
plotdap()
plotdap("base")

## ----southPole, fig.align = 'center', fig.height = 4, fig.width = 5-----------
plotdap("base",
  mapTitle = "MODIS South Pole Stereographic", 
  mapFill = "transparent", 
  mapColor = "steelblue",
  crs = "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=-63 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
)

## ----alaska, fig.align = 'center', fig.height = 4, fig.width = 5--------------
alaska <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
plotdap("base", crs = alaska)

## ---- echo = FALSE------------------------------------------------------------
alaska <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

## ----usmap, fig.align = 'center', fig.height = 4, fig.width = 5---------------
library(sf)
library(mapdata)
w <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
us <- st_transform(subset(w, ID == "USA"), alaska)
plotdap(mapData = us)

## ----get_viirsSST, eval = FALSE-----------------------------------------------
#  sstInfo <- rerddap::info('erdVHsstaWS3day')
#  # get latest 3-day composite sst
#  viirsSST <- rerddap::griddap(sstInfo,
#                               latitude = c(41., 31.),
#                               longitude = c(-128., -115),
#                               time = c('last','last'),
#                               fields = 'sst')
#  

## ----viirs_hires,echo = TRUE, eval = FALSE------------------------------------
#  w <- map("worldHires", xlim = c(-140., -114), ylim = c(30., 42.),
#           fill = TRUE, plot = FALSE)
#  # map using that outline,  temperature color from cmocean
#  add_griddap(plotdap(mapData = w), viirsSST, ~sst, fill = "thermal" )
#  

## ----world2Hires, fig.align = 'center', fig.height = 3, fig.width = 5,  message = FALSE, warning = FALSE----
xpos <- c(135.25, 240.25)
ypos <- c(20.25, 60.25)
zpos <- c(70.02, 70.02)
remove <- c("UK:Great Britain", "France", "Spain", "Algeria", "Mali", "
            Burkina Faso", "Ghana", "Togo")
#subset world2Hires with those countries removed
w <- map("world2Hires", plot = FALSE, fill = TRUE, ylim = ypos, xlim = xpos)
w <- map("world2Hires", regions = w$names[!(w$names %in% remove)], 
         plot = FALSE, fill = TRUE, ylim = ypos, xlim = xpos)
# plot result
plotdap(mapData = w)


## ----cairo, eval = FALSE,  echo = TRUE, warning = FALSE, message = FALSE------
#  # write plot to disk using the Cairo package
#  library(Cairo)
#  # (latitude limits) / (longitude limits)
#  r <- 85 / 120
#  CairoPNG("myPlot.png", height = 400 * r, width = 400, res = 96)
#  # alter default margins for base plotting
#  # (leaving just enough space for a title)
#  par(mar = c(0, 0, 1, 0))
#  plotdap("base", mapData = us, mapTitle = "Albers projection of Alaska")
#  dev.off()

## ----get_sardines, echo = TRUE, eval = FALSE----------------------------------
#  my_url <- 'https://coastwatch.pfeg.noaa.gov/erddap/'
#  sardines <- tabledap(
#    'FRDCPSTrawlLHHaulCatch',
#    fields = c('latitude',  'longitude', 'time', 'scientific_name',
#               'subsample_count'),
#    'time>=2010-01-01', 'time<=2012-01-01', 'scientific_name="Sardinops sagax"',
#     url = my_url)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  p1 <- add_tabledap(
#    plotdap(crs = "+proj=robin",  mapTitle = "subsample count"),
#    sardines,
#    ~subsample_count
#  )
#  p2 <- add_tabledap(
#    plotdap(crs = "+proj=robin", mapTitle = "Log subsample count"),
#    sardines,
#    ~log2(subsample_count)
#  )
#  
#  p1
#  p2
#  

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  p1 <- add_tabledap(
#    plotdap(crs = "+proj=robin", mapTitle = "Sardines - change color"),
#    sardines,
#    ~subsample_count,
#    color = "dense",
#  )
#  p2  <- add_tabledap(
#    plotdap(crs = "+proj=robin", mapTitle = "Sardines - change shape and size"),
#    sardines,
#    ~subsample_count,
#    shape = 4,
#    size = 1.
#  )
#  p1
#  p2

## ----get_mur, eval = FALSE----------------------------------------------------
#  murSST_west <- griddap(
#    'jplMURSST41',
#    latitude = c(22, 51),
#    longitude = c(-140, -105),
#    time = c('last', 'last'),
#    fields = 'analysed_sst'
#    )

## ----plot_mur, echo = TRUE, eval = FALSE--------------------------------------
#  add_griddap(
#    plotdap(crs = "+proj=robin"),
#    murSST_west,
#    ~analysed_sst,
#    maxpixels = 50000
#  )

## ----get_wind, eval = FALSE---------------------------------------------------
#  wind <- griddap(
#    'erdQMwindmday',
#    time = c('2016-04-16', '2016-06-16'),
#    latitude = c(30, 50),
#    longitude = c(210, 240),
#    fields = 'y_wind'
#  )

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  p1 <- add_griddap(
#    plotdap(mapTitle = "Mean Meridional Wind"),
#    wind,
#    ~y_wind,
#    fill = "delta",
#    time = mean
#  )
#  my_func <- function(x) var(x, na.rm = TRUE)
#  p2 <- add_griddap(
#    plotdap(mapTitle = "Variance of Meridional Wind"),
#    wind,
#    ~y_wind,
#    fill = "delta",
#    time = my_func
#  )
#  p1
#  p2
#  

## ----viirsSST_gridland, echo = TRUE, eval = FALSE-----------------------------
#  plotdap(mapTitle = "Grid over Land") %>%
#        add_griddap(
#          viirsSST,
#          ~sst,
#          fill = "thermal"
#          )
#  

## ----viirsSST_landgrid, echo = TRUE, eval = FALSE-----------------------------
#  plotdap(mapTitle = "Land Over Grid") %>%
#      add_griddap(
#        viirsSST,
#        ~sst,
#        fill = "thermal"
#        ) %>%
#      print(landmask = TRUE)

## ----maxpix_10, echo = TRUE, eval = FALSE-------------------------------------
#  plotdap(mapTitle = "maxpixels = 10,000") %>%
#      add_griddap(
#        viirsSST,
#        ~sst,
#        fill = "thermal",
#        maxpixels = 10000
#        )

## ----maxpix_50, echo = TRUE, eval = FALSE-------------------------------------
#  plotdap(mapTitle = "maxpixels = 50,000") %>%
#       add_griddap(
#         viirsSST,
#         ~sst,
#         fill = "thermal",
#         maxpixels = 50000
#         )

## ----maxpix_100, echo = TRUE, eval = FALSE------------------------------------
#  plotdap(mapTitle = "maxpixels = 100,000") %>%
#      add_griddap(
#        viirsSST,
#        ~sst,
#        fill = "thermal",
#        maxpixels = 100000
#        )

## ----maxpix_100_mask, echo = TRUE, eval = FALSE-------------------------------
#  plotdap(mapTitle = "maxpixels = 100,000, landmask") %>%
#      add_griddap(
#        viirsSST,
#        ~sst,
#        fill = "thermal",
#        maxpixels = 100000
#        ) %>%
#      print(landmask = TRUE)

## ----hires--------------------------------------------------------------------
w <- map("worldHires", xlim = c(-130., -114), ylim = c(30., 42.), 
         fill = TRUE, plot = FALSE)

## ----viirs_hires_mask, echo = TRUE, eval = FALSE------------------------------
#  plotdap(mapData = w) %>%
#      add_griddap(
#        viirsSST,
#        ~sst,
#        fill = 'thermal',
#        maxpixels = 50000
#        ) %>%
#      print(landmask = TRUE)

## ----soda70_get, echo = TRUE, eval = FALSE------------------------------------
#  soda70Info <- rerddap::info('erdSoda331oceanmday')
#  xpos <- c(135.25, 240.25)
#  ypos <- c(20.25, 60.25)
#  zpos <- c(70.02, 70.02)
#  tpos <- c('2010-12-15', '2010-12-15')
#  soda70 <- rerddap::griddap(soda70Info,
#                    longitude = xpos,
#                    latitude = ypos,
#                    time = tpos,
#                    depth = zpos,
#                    fields = 'temp' )
#  

## ----soda70, echo= TRUE, eval = FALSE-----------------------------------------
#  remove <- c("UK:Great Britain", "France", "Spain", "Algeria", "Mali",
#              "Burkina Faso", "Ghana", "Togo")
#  #subset world2Hires with those countries removed
#  w <- map("mapdata::world2Hires", plot = FALSE, fill = TRUE,
#           ylim = ypos, xlim = xpos)
#  w <- map("mapdata::world2Hires", regions = w$names[!(w$names %in% remove)],
#           plot = FALSE, fill = TRUE, ylim = ypos, xlim = xpos)
#  # plot result
#  plotdap(mapData = w) %>%
#          add_griddap(
#            soda70,
#            ~temp,
#            fill = "thermal"
#            ) %>%
#          print(landmask = TRUE)
#  

## ----overlay, echo = TRUE, eval = FALSE---------------------------------------
#  p <-   add_griddap(plotdap(),
#      murSST_west,
#      ~analysed_sst
#      )
#  p1 <-   add_tabledap(p,
#      sardines,
#      ~subsample_count
#      )
#  print(p1, landmask = TRUE)

## ----modify, fig.align = 'center', fig.height = 4, fig.width = 5--------------
library(ggplot2)
plotdap(crs = "+proj=robin") %>%
  add_tabledap(
    sardines, 
    ~subsample_count, 
    size = 1
    ) %>%
  add_ggplot(
    labs(
      subtitle = "Sardinops sagax samples",
      caption = "Sardines are yummy"
    ), 
    theme_minimal(),
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  )

## ----modify2, echo = TRUE, eval = FALSE, timeit = TRUE------------------------
#  temp_color <- rerddap::colors$temperature
#  plotdap(mapTitle = "Reset colorscale limits") %>%
#      add_griddap(
#        viirsSST,
#        ~sst,
#        fill = "thermal"
#        ) %>%
#      add_ggplot(
#         scale_fill_gradientn(colours = temp_color, na.value = NA, limits = c(10, 20)),
#          scale_colour_gradientn(colors = temp_color, na.value = NA, limits = c(10, 20)),
#          guides(colour = FALSE)
#       ) %>%
#       print(landmask = TRUE)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  add_tabledap(
#    plotdap(crs = "+proj=robin"),
#    sardines,
#    ~subsample_count,
#    color = "dense",
#    shape = 4,
#    animate = TRUE
#  )

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  add_griddap(
#    plotdap(crs = "+proj=robin"),
#    wind, ~y_wind,
#    time = identity,
#    fill = 'delta',
#    animate = TRUE
#  )

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  p <- add_griddap(
#    plotdap(crs = "+proj=robin"),
#    wind, ~y_wind,
#    time = identity,
#    fill = 'delta',
#    animate = TRUE
#  )
#    ylim <- c(30, 50)
#    xlim <- c(-150, -120)
#    p <- bbox_set(p, xlim, ylim)

