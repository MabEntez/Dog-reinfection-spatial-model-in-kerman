rm(list=ls())  #Erase workspace
library(sf)
library(sp)
library(rstudioapi)
library(dplyr)
library(stringr)

setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

# Loading and processing data ##### 

kerman_locations <- read.csv("./Kerman Dogs.csv")
d <- data.frame(locations = c(kerman_locations[,"Sampled.location"], kerman_locations[1:44,"Locations.NOT.sampled"]))
d[1:83,2] <- NA
d[1:83,3] <- kerman_locations[1:83,2]
d[1:83,4] <- kerman_locations[1:83,3]
d[is.na(d)] <- 0
d[,5] <- d[,4] / d[,3]

d[, 1:2] <- str_split_fixed(d$locations, ", ", 2)
d <- stringr::str_split_fixed(locations, c("Long", "lat"))

d <- d %>% 
  rename(
    Latitude = locations,
    Longitude = V2,
    N_samples = V3,
    N_positives = V4,
    Prevelance = V5,
  )

d$Longitude <- as.numeric(d$Longitude)
d$Latitude <- as.numeric(d$Latitude)


d <- d[!d$N_samples==0,]

sum(d$N_positives) / sum(d$N_samples)

# Plotting the sample locations#####
library(leaflet)
library(raster)

r <- raster(ncol=150, nrow=150, xmn=56.8, xmx=57.4, ymn=30, ymx=30.6)
r[] <- 1:length(r)
plot(r)
pal <- colorNumeric("viridis", c(0, 1), na.color = "grey")
leaflet(d) %>%  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~Longitude, lat = ~Latitude, color = ~pal(Prevelance)) %>%
  addLegend("bottomright", pal = pal, values = ~Prevelance, title = "% of tests positive") %>%
  addScaleBar(position = c("bottomleft"))

# Creating a mesh around the sampled region #####
library(INLA)
coo <- cbind(d$Longitude, d$Latitude)
mesh <- inla.mesh.2d(loc = coo, max.edge = c(0.02, 3), cutoff = 0.005)

mesh$n

plot(mesh)
points(coo, col = "red")

# SPDE INLA Model #####
spde <- inla.spde2.matern(mesh = mesh, alpha = 2)

indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

A <- inla.spde.make.A(mesh = mesh, loc = coo)

dp <- rasterToPoints(r)
dim(dp)

coop <- dp[, c("x", "y")]


Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

#stack for estimation stk.e
stk.e <- inla.stack(tag = "est",
                    data = list(y = d$N_positive, numtrials = d$N_samples),
                    A = list(1, A),
                    effects = list(data.frame(b0 = rep(1, times = nrow(d))), s = indexs))

#stack for prediction stk.p
stk.p <- inla.stack(tag = "pred",
                    data = list(y = NA, numtrials = NA),
                    A = list(1, Ap),
                    effects = list(data.frame(b0 = rep(1, times = nrow(dp))), s = indexs))

#stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ 0 + b0 + f(s, model = spde)

res <- inla(formula, family = "binomial", Ntrials = numtrials,
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))

summary(res)
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]

save.image("INLAData_all.RData")


#Visualising model outputs #####
load("INLAData_all.RData")

r_prev_mean <- rasterize(x = coop, y = r, field = prev_mean, fun = mean)

pal <- colorNumeric("viridis", c(0.37,0.51), na.color = "transparent")

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_mean), title = "Relative Risk") %>%
  addScaleBar(position = c("bottomleft"))

r_prev_ll <- rasterize(x = coop, y = r, field = prev_ll, fun = mean)

pal <- colorNumeric("viridis", c(0.11,0.35), na.color = "transparent")

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_ll, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_ll), title = "LL") %>%
  addScaleBar(position = c("bottomleft"))

r_prev_ul <- rasterize(x = coop, y = r, field = prev_ul, fun = mean)

pal <- colorNumeric("viridis", c(0.50,0.81), na.color = "transparent")

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_ul, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_ul), title = "UL") %>%
  addScaleBar(position = c("bottomleft"))


prev_rr <- prev_mean / 0.4154947
round(prev_rr, digits = 5)
prev_rr <- replace(prev_rr, prev_rr > 0.99 & prev_rr < 1.01, NA)

r_prev_rr <- rasterize(x = coop, y = r, field = prev_rr, fun = mean)


pal <- colorNumeric("viridis", c(0.9,1.21), na.color = "transparent")

m <- leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_rr, colors = pal, opacity = 0.6) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_rr), title = "Relative Risk") %>%
  addScaleBar(position = c("bottomleft"))%>% 
  fitBounds(56.9, 30.15, 57.3, 30.55)




library(mapview)


mapshot(m, file = "./Rplot.png")




rang <- apply(mesh$loc[, c(1, 2)], 2, range)
proj <- inla.mesh.projector(mesh,
                            xlim = rang[, 1], ylim = rang[, 2],
                            dims = c(300, 300)
)

mean_s <- inla.mesh.project(proj, res$summary.random$s$mean)
sd_s <- inla.mesh.project(proj, res$summary.random$s$sd)

df <- expand.grid(x = proj$x, y = proj$y)
df$mean_s <- as.vector(mean_s)
df$sd_s <- as.vector(sd_s)

library(viridis)
library(cowplot)
library(ggplot2)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

plot_grid(gmean, gsd)
