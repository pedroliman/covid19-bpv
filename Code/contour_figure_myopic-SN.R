
# COuntour Figure
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(fields)
library(viridis)
library(rootSolve)
# Reading data from the Tableau folder:

# Reading all data:


fR0 <- function(x,R0){exp(-R0*(1-x))-x}


Sinf1.5 <- multiroot.1D(f=fR0,R0=1.5,start = .5,nspec= 1)[["root"]]
Sinf3 <- multiroot.1D(f=fR0,R0=3,start = .5,nspec= 1)[["root"]]
Sinf4.5 <- multiroot.1D(f=fR0,R0=4.5,start = .5,nspec= 1)[["root"]]
Sinf6 <- multiroot.1D(f=fR0,R0=6,start = .5,nspec= 1)[["root"]]


S1.5max <- 1/1.5
S3max <- 1/3
S4.5max <- 1/4.5
S6max <- 1/6

epi.size <- cbind(R0 = c(1.5,3,4.5,6),EpiInf = c(1-Sinf1.5,1-Sinf3,1-Sinf4.5,1-Sinf6),EpiMin = c(1-S1.5max,1-S3max,1-S4.5max,1-S6max))

all_cases = read.csv("INTERMEDIATE/sol.explore.all.csv")
all_cases2 <- read.csv("INTERMEDIATE/parallel.output2.csv")
all_cases3 <- read.csv("INTERMEDIATE/parallel.output3.csv")
all_casesR03 <- read.csv("INTERMEDIATE/output-R03.csv")
all_casesR045 <- read.csv("INTERMEDIATE/output-R045.csv")

all_cases4 <- rbind(all_cases,all_cases2,all_cases3)
all_cases4 <-  subset(all_cases4,solved==TRUE&R0!=3&R0!=4.5)
all_cases4 <- rbind(all_cases4,all_casesR03[,1:dim(all_cases4)[2]],all_casesR045[,1:dim(all_cases4)[2]])

#all_cases = all_cases %>% 
#  dplyr::filter(R0 %in% c(1.5,3,4.5))


#data = read.csv("./Code/INTERMEDIATE/sol.explore-step11.csv") %>%
#  dplyr::select(tau, s, i, R_tau, Rt, case.id)

#exp_design = read.csv("./Code/INTERMEDIATE/exper.design-step11.csv") %>%
#  dplyr::select(case.id, c, R0, i0, tFinal)

# Selecting data fror the contour plot:




## Myopic

myopic_data = read.csv("INTERMEDIATE/sol.explore-myopic20210730.csv") %>%
  dplyr::select(tau, s, i, R_tau, Rt, case.id)

myopic_exp_design = read.csv("INTERMEDIATE/exper.design-round13.csv") %>%
  dplyr::select(case.id, c, R0, i0, tFinal)

myopic_data <- left_join(myopic_data,myopic_exp_design)

# Selecting data fror the contour plot:
myopic_final_data = myopic_data %>%
  dplyr::filter(tFinal == tau&case.id>417) %>%
  mutate(EpiSize = 1-s) 


g <- myopic_final_data %>%
  ggplot(data = ., mapping = aes(x = tFinal, y = c, fill = EpiSize)) +
  geom_tile() + scale_fill_viridis_c()+
  facet_wrap(facets = ~R0)

g

ggsave(filename = "FIGS/myopic_tile_plot.png",device = "png", plot = g, width = 3, height = 2.5, units = "in", scale = 2.5)

g2 <-  all_cases4 %>%
  dplyr::filter(tFinal == tau, case.id > 417) %>%
  mutate(EpiSize = 1-s) %>%
  ggplot(data = ., mapping = aes(x = tFinal, y = c, fill = EpiSize)) +
  geom_tile() + scale_fill_viridis_c()+
  facet_grid(i0 ~R0)

ggsave(filename = "FIGS/tile_plot.png",device = "png", plot = g2, width = 3, height = 2.5, units = "in", scale = 2.5)

g3 <-  all_cases4 %>%
  dplyr::filter(tFinal == tau, case.id > 417, R0==1.5, i0 == 1e-4) %>%
  mutate(EpiSize = 1-s) %>%
  ggplot(data = ., mapping = aes(x = tFinal, y = c, fill = EpiSize)) +
  geom_tile() + scale_fill_viridis_c()
ggsave(filename = "tile_plot_1pt5.png",device = "png", plot = g3, width = 3, height = 2.5, units = "in", scale = 2.5)


aux <- subset(all_cases4,R0 == 6 & c == 0.1 & i0 ==0.005)

g3 <- aux %>% ggplot(data = .,mapping = aes(x = tau, y = s,group = tFinal,color = tFinal)) + scale_colour_viridis_c()+
  geom_point()

ggsave(filename = "problems.png",device = "png", plot = g3, width = 3, height = 2.5, units = "in", scale = 2.5)


# Defining the Contour function

loess_countour_plot = function(results, x_variable = "tFinal",
                               y_variable = "c",
                               dependent_variable = "EpiSize",
                               facet_variable = "R0",
                               facet_vector = c(2, 4),
                               binwidth = 100,
                               nudge_y = 0,
                               nudge_x = 0,
                               skip = 0) {
  
  
  
  if(is.null(facet_variable)) {
    selected_results = results[,c(x_variable, y_variable, dependent_variable)]
    names(selected_results) = c("x","y", "z")
    loess_model = loess(z ~ x + y, data = selected_results)
  } else {
    selected_results = results[,c(x_variable, y_variable, facet_variable, dependent_variable)] 
    names(selected_results) = c("x","y","facet", "z")
    loess_model = loess(z ~ x + y + facet, data = selected_results)
  }
  
  x_vector= seq.default(from = min(selected_results$x), to = max(selected_results$x), length.out = 100)
  y_vector= seq.default(from = min(selected_results$y), to = max(selected_results$y), length.out = 100)
  
  
  if(is.null(facet_variable)) {
    meta_model_results <-  expand.grid(x_vector, y_vector)
    names(meta_model_results) = c("x","y")
    
    mtrx3d =  predict(loess_model, meta_model_results)
    
    # Transform data to long form
    mtrx.melt <- melt(mtrx3d, id.vars = c("x", "y"), measure.vars = "z")
    names(mtrx.melt) <- c("x", "y", "z")
    # Return data to numeric form
    mtrx.melt$x <- as.numeric(str_sub(mtrx.melt$x, str_locate(mtrx.melt$x, "=")[1,1] + 1))
    mtrx.melt$y <- as.numeric(str_sub(mtrx.melt$y, str_locate(mtrx.melt$y, "=")[1,1] + 1))
    
  } else {
    meta_model_results <-  expand.grid(x_vector, y_vector, facet_vector)
    names(meta_model_results) = c("x","y","facet")
    mtrx3d =  predict(loess_model, meta_model_results)
    
    # Transform data to long form
    mtrx.melt <- melt(mtrx3d, id.vars = c("x", "y", "facet"), measure.vars = "z")
    names(mtrx.melt) <- c("x", "y", "facet", "z")
    # Return data to numeric form
    mtrx.melt$x <- as.numeric(str_sub(mtrx.melt$x, str_locate(mtrx.melt$x, "=")[1,1] + 1))
    mtrx.melt$y <- as.numeric(str_sub(mtrx.melt$y, str_locate(mtrx.melt$y, "=")[1,1] + 1))
    mtrx.melt$facet <- as.numeric(str_sub(mtrx.melt$facet, str_locate(mtrx.melt$facet, "=")[1,1] + 1))
  }
  
  
  # We should not find NAs in the mtrx.melt dataframe:
  if(any(is.na(mtrx.melt))){
    # Found NAs in the mtrx.melt data.frame. Check inputs.
    browser()
  }
  
  
  ### Results from the metamodel:
  
  plot_function = function(data, facet_variable) {
    
    if(is.null(facet_variable)) {
      names(data) = c("x","y","z")
    } else {
      names(data) = c("x","y",facet_variable, "z")
    }
    
    plot <- ggplot(data, aes(x, y)) +
      geom_tile(aes(fill=z)) + 
      ggplot2::geom_contour(aes(z = z), color = "black",  binwidth = binwidth, size = 1)
    
    
    if(!is.null(facet_variable)) {
      plot = plot + facet_wrap(facets = facet_variable, labeller = label_both)
    } 
    
    # computing breaks:
    bins = (max(data$z) - min(data$z))/binwidth
    breaks = pretty(x = data$z,n = bins)
    
    plot = plot +
      scale_fill_gradientn(colours=tim.colors(128), breaks = breaks) + 
      # Adjusts Contour text positions
      metR::geom_text_contour(aes(z = z),show.legend = T,stroke.color = "white", binwidth = binwidth, skip = skip, nudge_y = nudge_y, nudge_x = nudge_x) + 
      #hrbrthemes::theme_ipsum_ps(axis_title_just = "c") + 
      ggpubr::theme_pubclean() + 
      theme(legend.position="bottom") + 
      theme(axis.title = element_text(face=1))
    
    plot
    
  }
  
  mtrx.melt %>%
    plot_function(data = ., facet_variable = facet_variable)
  
} 

loess_plot = loess_countour_plot(results = final_data,facet_vector = c(2,3), binwidth = 0.2, skip = 0, nudge_y = 0.02) + 
  ylab("Perceived social distancing cost") + 
  xlab("Epidemic time-frame") + 
  hrbrthemes::theme_ipsum_ps(axis_title_just = "c") + 
  theme(legend.position="bottom") + 
  labs(fill = "Epidemic size") + 
  theme(legend.key.width = unit(1, "in"))


loess_plot

ggsave(filename = "contour_plot.png",device = "png", plot = loess_plot, width = 3, height = 2.5, units = "in", scale = 2.5)


