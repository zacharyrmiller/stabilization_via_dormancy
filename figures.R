# Contact: zachary.miller@yale.edu
# Simulation and visualization code for
# "Stabilization of fluctuating population dynamics via the evolution of dormancy"
# by Z.R. Miller, D. Vasseur, and P.M. Hull

library(tidyverse)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(deSolve)

viridis_pal <- viridis(10) # define a palette to use throughout

##### Simulation functions #####

general_dormancy_dynamics <- function(r1, r2 = r1, # intrinsic growth rates 
                                      a1, a2 = a1, # active fractions (controlling fraction dormant)
                                      m1, m2 = m1, # mortality rates in dormancy
                                      x1_init, x2_init = 0, # initial population sizes
                                      t_max, # number of time steps to iterate dynamics
                                      dd_function = logistic_map # function defining density-dependence of net growth rates
                                      ){
  
  # General function for numerically iterating population dynamics of two ecotypes with dormancy
  # Dynamics correspond to Eq. S11
  # For dynamics of one ecotype, use x2_init = 0

  # Initialize populations
  x1 <- x2 <- vector(length = t_max)
  x1[1] <- x1_init
  x2[1] <- x2_init
  
  # Iterate dynamics
  for(i in 1:(t_max-1)){
    
    x1[i+1] <- r1 * a1 * x1[i] * dd_function(a1 * x1[i] + a2 * x2[i]) + (1 - a1) * (1 - m1) * x1[i]
    x2[i+1] <- r2 * a2 * x2[i] * dd_function(a1 * x1[i] + a2 * x2[i]) + (1 - a2) * (1 - m2) * x2[i]
    
  }
  
  # Organize output
  ts <- tibble(t = 1:t_max, x1 = x1, x2 = x2) # convert to tibble
  ts <- ts %>% pivot_longer(cols = -t, names_to = "type") # pivot to tidy form (for convenient plotting)
  if(x2_init == 0) ts <- ts %>% filter(type != "x2") # if second ecotype not present, remove from the output
  
  return(ts)
}


glv_dormancy <- function(t, x, parameters) {
  with(as.list(c(x, parameters)), {
    
    # general function for generalized Lotka-Volterra dynamics with dormancy
    # (for numerical integration with deSolve)
    
    A <- matrix(A, nrow = 3 * n) # interaction matrix
    B <- matrix(B, nrow = 3 * n) # dormancy transition matrix
    
    dx <- x * (r + A %*% x) + B %*% x # after Eq. 6
    
    return(list(dx))
  })
}

### Common forms of density-dependence (with arbitrary parameterizations)

logistic_map <- function(x){
  1 - x
}

ricker_map <- function(x){
  exp(-x)
}

hassell_map <- function(x){
  1 / (1 + x)^4
}

maynardsmith_map <- function(x){
  1 / (1 + x^4)
}


##### Manuscript figures #####

### Conceptual figure (Fig. 1)

## Panel A: Illustration of common (separable) density-dependence functions
r <- 3
df_models <- tibble(N = seq(0, 1.5, by = 0.01),
                    logistic = r * logistic_map(N),
                    ricker = r * ricker_map(N),
                    hassell = r * hassell_map(N),
                    maynardsmith = r * maynardsmith_map(N))

pa <- df_models %>% 
  pivot_longer(cols = -N, names_to = "model", values_to = "growth_rate") %>% # pivot to tidy form
  ggplot() + 
  aes(x = N, y = growth_rate, group = model, color = model) +
  geom_line(size = 1.2) + 
  annotate("text", x = c(0.3, 0.95, 0.95, 0.4), y = c(0.6, 0.7, 2.5, 2.4), 
           label = c("Hassell", "Logistic", "Maynard-Smith", "Ricker"),
           color = brewer.pal(4, "Dark2"),
           size = 4) + 
  xlab(expression(N[t])) + ylab(expression(f(N[t]))) + 
  scale_x_continuous(expand = c(0, 0.01)) + 
  scale_y_continuous(expand = c(0, 0.05), limits = c(0, 3)) + 
  scale_color_manual(values = brewer.pal(4, "Dark2")) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 14, face = "bold"))

# ggsave(filename = "./figures/models_conceptual.png", plot = pa, device = "png", 
#       width = 4, height = 3, units = "in", dpi = 400)

## Panel B: Illustration of population map with dormancy (using logistic density-dependence)
df_alphas <- tibble(N = seq(0, 1.7, by = 0.001), # generate maps for different choices of alpha
                    a_1 = N * r * logistic_map(N),
                    a_0.9 = 0.9 * N * r * logistic_map(0.9 * N) + (1 - 0.9) * (1 - 0.1) * N,
                    a_0.7 = 0.7 * N * r * logistic_map(0.7 * N) + (1 - 0.7) * (1 - 0.1) * N,
                    a_0.8 = 0.8 * N * r * logistic_map(0.8 * N) + (1 - 0.8) * (1 - 0.1) * N)
  
pb <- df_alphas %>% 
  pivot_longer(cols = -N, names_to = "alpha", values_to = "N_plus") %>% # pivot to tidy form
  ggplot() + 
  aes(x = N, y = N_plus, group = alpha, color = alpha) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + # add one-to-one line
  geom_line(size = 1.2) + 
  annotate("text", 
           x = c(1.45, 1.28, 1.15, 0.85), 
           y = c(0.5, 0.35, 0.2, 0.08), 
           label = c(0.7, 0.8, 0.9, "paste(alpha, \" = 1\")"),
           parse = TRUE,
           color = colorRampPalette(c("#F2C9AA", "#D95F02"))(4)) + # match scale_color palette (below)
  xlab(expression(N[t])) + ylab(expression(N[t+1])) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0.01), limits = c(0, 1)) + 
  scale_color_manual(values = colorRampPalette(c("#F2C9AA", "#D95F02"))(4)) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.title=element_text(size = 14, face = "bold"))

# ggsave(filename = "./figures/alphas_conceptual.png", plot = pb, device = "png", 
#       width = 4, height = 3, units = "in", dpi = 400)


### Logistic model bifurcation diagram (Fig. 2)

r_seq <- seq(0.01, 5, 0.001)
bif_df <- tibble(r = numeric(),
                 values = numeric())

# Iterate logistic model without dormancy to generate standard bifurcation diagram (for panel A)
for(r in r_seq){
  
  # Iterate dynamics for 1000 timesteps to determine attractor
  ts <- general_dormancy_dynamics(r1 = r, a1 = 1, m1 = 0, 
                            x1_init = max(0, (r - 1) / r + 0.01), # start near equilibrium for faster convergence
                            x2_init = 0, 
                            t_max = 1000, 
                            dd_function = logistic_map)
  
  ts <- ts %>% replace(is.na(.), 0)
  bif_df <- bif_df %>% add_row(r = r, values = unique(tail(ts$value, 100))) # save unique values from last 100 time steps
}

# Now use r_eff (Eq. 3) to classify dynamics for different combinations of r and alpha
df <- tibble(expand.grid(a = seq(0.01, 1, by = 0.001), 
                         r = r_seq))

m <- 0.01 # fix m at a small value

df <- df %>% 
  mutate(r_eff = a * r + (1 - a) * (1 - m), # calculate r_eff (Eq. 3)
         dyn = ifelse(r_eff < 1, 0, # classify dynamics using known bifurcation points for the logistic map (see e.g. Devenay 1986/2019)
                      ifelse(r_eff < 3, 1,
                             ifelse(r_eff < 3.4495, 2,
                                    ifelse(r_eff < 3.544, 4,
                                           ifelse(r_eff < 3.56995, 8,
                                                  ifelse(r_eff <= 4, Inf, 0)))))))

## Panel A: Plot standard bifurcation diagram

pa <- bif_df %>% 
  left_join(., df %>% filter(a == 1), by = "r") %>% # join data frames (to color by qualitative dynaimcs)
  filter(r >= 1) %>% # only plot r values where the population can grow
  ggplot() + 
  aes(x = r, y = values, color = as.factor(dyn)) + 
  geom_point(size = 0.1) + 
  xlab("Growth rate (r)") + 
  ylab("Population size") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c(viridis_pal[1],
                               viridis_pal[4],
                               viridis_pal[7:10])) +
  guides(color = FALSE) +
  theme_classic()

## Panel B: Qualitative dynamics as a function of alpha and r

pb <- df %>% 
  filter(a > 0, r >= 1) %>%
  ggplot() + 
  aes(x = r, y = a, fill = as.factor(dyn)) + 
  geom_raster() + 
  xlab("Growth rate (r)") + ylab(expression(Active~fraction~(alpha))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = c(viridis_pal[1],
                               viridis_pal[4],
                               viridis_pal[7:10]),
                    labels = c("Extinction",
                               "Stable",
                               "Period 2",
                               "Period 4",
                               "Period 8+",
                               "Chaos")) + 
  guides(fill = guide_legend(override.aes = list(size = 0.4))) + 
  theme_bw() + 
  theme(legend.position = c(0.65, 0.1), 
        legend.direction = "horizontal",
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.margin = margin(c(2,2,2,2)))

# combine panels
p <- ggarrange(plotlist = list(pa, pb), nrow = 2, heights = c(0.5, 1), 
          labels = "AUTO", font.label = list(face = "plain"))

# ggsave(filename = "./figures/bifurcation_diagram.png", plot = p, device = "png", 
#      width = 4.5, height = 6, units = "in", dpi = 400)


### Mutual invasion time series (Fig. 3)

## Main panel: invasion by high-dormancy ecotype

# set parameters (discussed in text)
r_res <- 3.9
r_inv <- 3.8
a_res <- 1
a_inv <- 0.7
m <- 0.05

t_est <- 50 # iterate resident dynamics for t_est steps so that invader enters with resident at steady state
t_max <- 175 # number of time steps for invasion time series

# iterate resident dynamics to reach steady state
establish_res <- general_dormancy_dynamics(r1 = r_res, r2 = r_inv, a1 = a_res, a2 = a_inv, m1 = m,
                                           x1_init = (a_res * r_res + (1 - a_res) * (1 - m) - 1) / (r_res * a_res^2) + 0.01, # start near equilibrium for faster convergence
                                           x2_init = 0, t_max = t_est)

# iterate dynamics with invader at low initial abundance
introduce_inv <- general_dormancy_dynamics(r1 = r_res, r2 = r_inv, a1 = a_res, a2 = a_inv, m1 = m,
                                           x1_init = max(tail(establish_res, 2)$value), 
                                           x2_init = 0.01, t_max = t_max)
ts <- introduce_inv

p_main <- ts %>% 
  ggplot() + 
  aes(x = t, y = value, group = type, color = type) + 
  geom_line() + 
  scale_y_log10() + 
  scale_color_manual(values = c(viridis_pal[10],
                     viridis_pal[4])) + 
  xlab("Time") + ylab("Population size") + 
  theme_classic() + 
  theme(legend.position = "none")

## Inset panel: Invasion by no-dormancy ecotype

# set parameters (reversed from main panel)
r_res <- 3.8
r_inv <- 3.9
a_res <- 0.7
a_inv <- 1
m <- 0.05

t_est <- 50 # iterate resident dynamics for t_est steps so that invader enters with resident at steady state
t_max <- 100# number of time steps for invasion time series

# iterate resident dynamics to reach steady state
establish_res <- general_dormancy_dynamics(r1 = r_res, r2 = r_inv, a1 = a_res, a2 = a_inv, m1 = m,
                                           x1_init = (a_res * r_res + (1 - a_res) * (1 - m) - 1) / (r_res * a_res^2) + 0.01, # start near equilibrium for faster convergence
                                           x2_init = 0, t_max = t_est)

# iterate dynamics with invader at low initial abundance
introduce_inv <- general_dormancy_dynamics(r1 = r_res, r2 = r_inv, a1 = a_res, a2 = a_inv, m1 = m,
                                           x1_init = max(tail(establish_res, 2)$value), 
                                           x2_init = 0.01, t_max = t_max)
ts <- introduce_inv

p_inset <- ts %>% 
  ggplot() + 
  aes(x = t, y = value, group = type, color = type) + 
  geom_line() + 
  scale_y_log10() + 
  scale_color_manual(values = c(viridis_pal[4],
                                viridis_pal[10])) + 
  xlab("Time") + ylab("Population size") + 
  theme_classic() + 
  theme(legend.position = "none")

# combine panels
p <- p_main + annotation_custom(ggplotGrob(p_inset), xmin = 75, xmax = 165, 
                       ymin = -2.4, ymax = -1.1)

# ggsave(filename = "./figures/invasion.png", plot = p, device = "png", 
#      width = 6, height = 4, units = "in", dpi = 400)


### Pairwise invasibility plots with logistic density-dependence (Fig. 4)                      

# Calculate IGRs for many combinations of resident and invader alpha values.
# To numerically estimate IGR, iterate resident dynamics for many time steps, then 
# identify a time point with population size as close as possible to the time point
# where the invasion begins -- the resident dynamics are approximately stationary
# over this sequence of time steps. Calculate IGR based on this sequence.

# produce PIPs for 4 combinations of r and m
r_vec <- c(3.3, 3.8)
m_vec <- c(0.01, 0.05)

# how many alpha values to test for resident and invader
res_n <- 500
inv_n <- 500

# generate sequence of alpha values, beginning with the smallest value for which the population can grow (r_eff > 1)
res_a_vec <- seq((1 - min(m_vec)) / (min(r_vec) - min(m_vec)) + 10/res_n, 1, length.out = res_n)
inv_a_vec <- seq((1 - min(m_vec)) / (min(r_vec) - min(m_vec)) + 10/inv_n, 1, length.out = inv_n)

t_max <- 1000 # number of timesteps for numerical iteration
t_inv <- t_max - 400 # timepoint at which to begin calculating IGR (allowing resident to reach steady state)

df <- tibble(res_a = numeric(), inv_a = numeric(), 
             r = numeric(), m = numeric(), 
             IGR = numeric())

for(r in r_vec){
  for(m in m_vec){
    
    for(i in 1:length(res_a_vec)){
      
      res_a <- res_a_vec[i]
      
      # iterate resident dynamics
      ts <- general_dormancy_dynamics(r1 = r, a1 = res_a, m1 = m, 
                                      x1_init = (res_a * r_res + (1 - res_a) * (1 - m) - 1) / (r * res_a^2) + 0.01, # start near equilibrium for faster convergence
                                      x2_init = 0, t_max = t_max)
      ts <- ts %>% pull(value)
      ts <- ts[t_inv:t_max] # pull late values (after transient dynamics have elapsed)
      
      # find a time point such that net population change is minimized
      end_pt <- which.min(abs(cumsum(log(res_a * r * (1 - res_a * ts) + (1 - res_a) * (1 - m)))))
      x <- ts[2:(end_pt + 1)] # extract stationary sequence of time points
      
      for(j in 1:length(inv_a_vec)){
        
        inv_a <- inv_a_vec[j]
        
        # calculate IGR over stationary resident dynamics
        IGR <- sum(log(r * inv_a * (1 - res_a * x) + (1 - inv_a) * (1 - m)))
        df <- df %>% add_row(res_a = res_a, inv_a = inv_a, r = r, m = m, IGR = IGR)
      }
    }
  }
}

p <- df %>% 
  filter(res_a != inv_a) %>%
  ggplot() + 
  aes(x = res_a, y = inv_a, fill = IGR > 0) + 
  geom_raster(alpha = 0.75) + 
  geom_abline(slope = 1, intercept = 0, size = 1, color = "gray") + # indicate neutral line
  geom_vline(data = . %>% group_by(res_a, r, m) %>% # add a vertical line for ESS
               summarize(max_IGR = max(IGR)) %>%
               group_by(r, m) %>% 
               filter(max_IGR == min(max_IGR)), 
             aes(xintercept = res_a),
              linetype = "dashed", size = 0.8) +
  geom_hline(data = . %>% group_by(res_a, r, m) %>% # add a horizontal line for ESS 
               summarize(max_IGR = max(IGR)) %>%
               group_by(r, m) %>% 
               filter(max_IGR == min(max_IGR)), 
             aes(yintercept = res_a),
             linetype = "dashed", size = 0.8) +
  geom_vline(aes(xintercept = (3 - (1 - m))/(r - (1 - m))), # add a vertical line for bifurcation point
             linetype = "11", color = "red", size = 1) + 
  annotate(geom = "text", x = 0.6, y = 0.45, label = "\u2013", size = 8) + 
  annotate(geom = "text", x = 0.45, y = 0.65, label = "+", size = 8) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = c(viridis_pal[10],
                                viridis_pal[4])) + 
  facet_grid(m ~ r, 
             labeller = labeller(m = function(x) paste0("m = ", x),
                                 r = function(x) paste0("r = ", x))) + 
  theme_classic() +
  xlab("Resident active fraction") + ylab("Invader active fraction") + 
  theme(legend.position = "none")

#ggsave(filename = "./figures/pip.png", plot = p, device = "png", 
#       width = 5, height = 5, units = "in", dpi = 400)


### Invasion of simple, chaotic food web by dormancy in continuous-time (Fig. 5, Box 1)

# Set up parameters for 3 spp food from Gilpin (1979) -- see also Robey et al. 2024

n <- 3 # number of species
A <- -1000 * matrix(c(0.001, 0.001, 0.01, 0.0015, 0.001, 0.001, -0.005, -0.0005, 0), 
                    byrow = TRUE, nrow = n) # 3 spp interaction matrix
r <- c(10, 10, -10) # growth rates

m <- 0.01 # morality rate in dormancy
p <- c(0, 0, 1.5) # transition rates (into dormancy)
q <- c(0, 0, 0.015) # transition rates (out of dormancy)

# build 3n x 3n interaction matrix (corresponding to ecotypes without (n) and with (n) dormancy + dormant states (n))
big_A <- rbind(cbind(A, A, matrix(0, n, n)),
               cbind(A, A, matrix(0, n, n)),
               cbind(matrix(0, n, 2 * n), -m * diag(n)))

# build 3n x 3n matrix of dormancy transition rates
big_B <- rbind(matrix(0, n, 3 * n),
               cbind(matrix(0, n, n), -diag(p), diag(q)),
               cbind(matrix(0, n, n), diag(p), -diag(q)))

big_r <- c(rep(r, times = 2), rep(0, n))


# initial population sizes for resident community (without dormancy)
x <- c(1, 10, 0.1,
       0, 0, 0,
       0, 0, 0)

parameters <- list(A = as.vector(big_A), B = as.vector(big_B), r = big_r, n = n)
times <- seq(0, 60, by = 0.1)

# integrate resident dynamics to reach steady state
out1 <- ode(y = x, time = times, func = glv_dormancy, parms = parameters, method = "ode45")

# re-initialize with current resident population sizes, but adding rare predator with dormancy
x <- c(out1[nrow(out1), c(2, 3, 4)],
       0, 0, 10^-20,
       0, 0, 0)

parameters <- list(A = as.vector(big_A), B = as.vector(big_B), r = big_r, n = n)
times <- seq(0, 100, by = 0.1)

# integrate dynamics with invader
out2 <- ode(y = x, time = times, func = glv_dormancy, parms = parameters, method = "ode45")

# combine time-series
out2[, 1] <- out2[, 1] + 60
out <- rbind(out1, out2)

df <- out[, c(1, 2, 3, 4, 7, 10)] # keep only non-zero populations for plotting

p <- df %>% 
  as_tibble() %>%
  pivot_longer(cols = -time) %>%
  filter(value > 0) %>%
  mutate(type = ifelse(name > 2 * n, "dormant", "active")) %>% # classify populations as dormant or active
  mutate(name = ifelse(name == 1, "V1", # assign species labels
                       ifelse(name == 2, "V2",
                              ifelse(name == 3, "P", 
                                     "P'")))) %>%
  ggplot() + 
  aes(x = time, y = log10(value), # plot log population dynamics
      group = interaction(name, type), 
      color = name, alpha = type) + 
  annotate(geom = "rect", xmin = 60, xmax = 95, ymin = -Inf, ymax = Inf, # annotate distinct intervals of the dynamics (see text)
           fill = "gray95") + 
  annotate(geom = "rect", xmin = 95, xmax = 120, ymin = -Inf, ymax = Inf,
           fill = "gray80") + 
  geom_vline(xintercept = 60, linetype = "dashed") + # indicate time point where invader is introduced
  geom_line(size = 1) + 
  xlab("Time") + ylab("Population size (log)") +
  facet_grid(name~., scales = "free_y") + 
  scale_y_continuous(breaks = scales::extended_breaks(n = 4)) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_alpha_discrete(range = c(1, 0.4)) + 
  theme_classic() + 
  theme(legend.position = "none")

# ggsave(filename = "./figures/box.png", plot = p, device = "png", 
#       width = 4, height = 3, units = "in", dpi = 400)


##### SI figures #####

### IGR as a function of alpha and m (for different choices of r) in the logistic model (Fig. S1)

# calculations as in PIPs (Fig. 4)

r_vec <- c(3.3, 3.5, 3.7, 3.9)

# how many alpha and m values to test
a_n <- 300
m_n <- 300

a_vec <- seq(0.001, 1, length.out = a_n)
m_vec <- seq(0.001, 0.5, length.out = m_n)

t_max <- 1000 # number of timesteps for numerical iteration
t_inv <- t_max - 400 # timepoint at which to begin calculating IGR (allowing resident to reach steady state)

df <- tibble(a = numeric(), 
             r = numeric(), 
             m = numeric(), 
             IGR = numeric())

for(r in r_vec){
  
  # iterate resident dynamics (always using resident with no dormancy)
  ts <- general_dormancy_dynamics(r1 = r, a1 = 1, m1 = 0, 
                                  x1_init = (r - 1)/r + 0.01, # start near equilibrium for fast convergence
                                  x2_init = 0, t_max = t_max)
  ts <- ts %>% pull(value)
  ts <- ts[t_inv:t_max] # pull late values (after transient dynamics have elapsed)
  
  # find a time point such that net population change is minimized
  end_pt <- which.min(abs(cumsum(log(r * (1 - ts)))))
  x <- ts[2:(end_pt + 1)] # extract stationary sequence of time points
  
  for(i in 1:length(a_vec)){
    
    inv_a <- a_vec[i]
    
    for(j in 1:length(m_vec)){
      
      m <- m_vec[j]
      
      # calculate IGR over stationary resident dynamics
      IGR <- (1 / length(x)) * sum(log(r * inv_a * (1 - x) + (1 - inv_a) * (1 - m))) 
      df <- df %>% add_row(a = inv_a, r = r, m = m, IGR = IGR)
    }
  }
}

p <- df %>%
  ggplot() + 
  aes(x = a, y = m, 
      fill = pmax(IGR, -0.1), color = pmax(IGR, -0.1)) + # color by log IGR, but truncate very small values for clearer visualization
  geom_tile(linewidth = 0.1) + 
  geom_line(data = df %>% filter(a < 1) %>% # indicate the threshold between positive and negative log IGR
              group_by(r, a) %>% 
              mutate(minIGR = min(abs(IGR))) %>% 
              filter(IGR == minIGR), 
            aes(x = a, y = m), color = "black") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = c(0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) + 
  facet_wrap(.~r,
             labeller = labeller(r = function(x) paste0("r = ", x))) + 
  theme_classic() + 
  xlab("Invader active fraction") + ylab("Mortality in dormancy") + 
  scale_fill_gradient2(high = viridis_pal[4], low = viridis_pal[10], 
                       midpoint = 0, name = "log IGR") + 
  scale_color_gradient2(high = viridis_pal[4], low = viridis_pal[10], 
                       midpoint = 0, name = "log IGR")

#ggsave(filename = "./figures/si_IGR_varying_m.png", plot = p, device = "png", 
#       width = 6, height = 5, units = "in", dpi = 400)


### Logistic model bifurcation diagram (Fig. S2)

# As Fig. 2, but with Ricker density-dependence

r_seq <- seq(0.1, 20, 0.01)
bif_df <- tibble(r = numeric(),
                 values = numeric())

for(r in r_seq){
  
  # Iterate Ricker model without dormancy to generate standard bifurcation diagram (for panel A)
  ts <- general_dormancy_dynamics(r1 = r, a1 = 1, m1 = 0, 
                                  x1_init = max(0, log(r) + 0.01), # start near equilibrium for faster convergence
                                  x2_init = 0, 
                                  t_max = 1000, 
                                  dd_function = ricker_map)
  
  ts <- ts %>% replace(is.na(.), 0)
  bif_df <- bif_df %>% add_row(r = r, values = unique(tail(ts$value, 100))) # save unique values from last 100 time steps
}

# Now, use numerical iteration to classify dynamics for different combinations of r and alpha
a_seq <- seq(0.1, 1, by = 0.005)
df <- tibble(a = numeric(), r = numeric(), dyn = numeric())

m <- 0.01 # fix m at a small value

for(r in r_seq){
  for(a in a_seq){
    
    # iterate resident dynamics
    ts <- general_dormancy_dynamics(r1 = r, a1 = a, m1 = 0, 
                                    x1_init = max(0, -log((1 - (1 - a)*(1-m)) / (a * r)) + 0.01), # start near equilibrium for faster convergence
                                    x2_init = 0, t_max = 2500, 
                                    dd_function = ricker_map)
    
    ts <- ts %>% replace(is.na(.), 0)
    x <- unique(round(tail(ts$value, 100), 2)) # save unique values from last 100 time steps (round to identify the limit set)
    
    # classify dynamics according to number of points in the limit set
    dyn <- ifelse(max(x) < 10^-6, 0,
                  ifelse(length(x) %in% c(1,2,4,8), length(x),
                          Inf))
    
    df <- df %>% add_row(a = a, r = r, dyn = dyn)
  }
}

## Panel A: Plot standard bifurcation diagram

pa <- bif_df %>% 
  left_join(., df %>% filter(a == 1), by = "r") %>% # join data frames (to color by qualitative dynaimcs)
  filter(r >= 1) %>% # only plot r values where the population can grow
  ggplot() + 
  aes(x = r, y = values, color = as.factor(dyn)) + 
  geom_point(size = 0.1) + 
  xlab("Growth rate (r)") + ylab("Population size") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) + 
  scale_color_manual(values = c(viridis_pal[1],
                                viridis_pal[4],
                                viridis_pal[7:10])) +
  guides(color = FALSE) +
  theme_classic()

## Panel B: Qualitative dynamics as a function of alpha and r

pb <- df %>% 
  filter(a > 0, r >= 1) %>%
  ggplot() + 
  aes(x = r, y = a, fill = as.factor(dyn)) + 
  geom_raster() + 
  xlab("Growth rate (r)") + 
  ylab(expression(Active~fraction~(alpha))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = c(viridis_pal[1],
                              viridis_pal[4],
                              viridis_pal[7:10]),
                   labels = c("Extinction",
                              "Stable",
                              "Period 2",
                              "Period 4",
                              "Period 8+",
                              "Chaos")) +
  guides(fill = guide_legend(override.aes = list(size = 0.4))) + 
  theme_bw() + 
  theme(legend.position = c(0.65, 0.1), 
        legend.direction = "horizontal",
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.margin = margin(c(2,2,2,2)))

# combine panels
p <- ggarrange(plotlist = list(pa, pb), nrow = 2, heights = c(0.5, 1), 
               labels = "AUTO", font.label = list(face = "plain"))

# ggsave(filename = "./figures/si_bifurcation_ricker.png", plot = p, device = "png", 
#       width = 4.5, height = 6, units = "in", dpi = 400)


### Pairwise invasibility plots with Ricker density-dependence (Fig. S3) 

# after Fig. 4

# produce PIPs for 4 combinations of r and m
r_vec <- c(10, 20)
m_vec <- c(0.01, 0.05)

# how many alpha values to test for resident and invader
res_n <- 500
inv_n <- 500

# generate sequence of alpha values
res_a_vec <- seq(10/res_n, 1, length.out = res_n)
inv_a_vec <- seq(10/inv_n, 1, length.out = inv_n)

t_max <- 1000 # number of timesteps for numerical iteration
t_inv <- t_max - 400 # timepoint at which to begin calculating IGR (allowing resident to reach steady state)

df <- tibble(res_a = numeric(), inv_a = numeric(), 
             r = numeric(), m = numeric(), 
             IGR = numeric())

# initialize a tibble to record bifurcation points as a function of alpha (for each value of r)
bif <- tibble(r = numeric(), a_c = numeric())

for(r in r_vec){
  
  bif_flag <- FALSE
  
  for(m in m_vec){
    
    for(i in 1:length(res_a_vec)){
      
      res_a <- res_a_vec[i]
      
      # iterate resident dynamics
      ts <- general_dormancy_dynamics(r1 = r, a1 = res_a, m1 = m, 
                                      x1_init = 1.01,
                                      x2_init = 0, t_max = t_max, 
                                      dd_function = function(x) ricker_map(x))
      ts <- ts %>% pull(value)
      ts <- ts[t_inv:t_max] # pull late values (after transient dynamics have elapsed)
      
      # find a time point such that net population change is minimized
      end_pt <- which.min(abs(cumsum(log(res_a * r * ricker_map(res_a * ts) + (1 - res_a) * (1 - m)))))
      x <- ts[2:(end_pt + 1)] # extract stationary sequence of time points
      
      # where the resident population first begins to oscillate, record the bifurcation point
      if((ifelse(length(x) > 1, var(x), 0) > 10^-6) & !bif_flag){ 
        bif_flag <- TRUE
        bif <- bif %>% add_row(r = r, a_c = res_a)
      }
      
      for(j in 1:length(inv_a_vec)){
        
        inv_a <- inv_a_vec[j]
        
        # calculate IGR over stationary resident dynamics
        IGR <- sum(log(r * inv_a * ricker_map(res_a * x ) + (1 - inv_a) * (1 - m)))
        df <- df %>% add_row(res_a = res_a, inv_a = inv_a, r = r, m = m, IGR = IGR)
      }
    }
  }
}

p <- df %>% 
  filter(res_a != inv_a) %>%
  ggplot() + 
  aes(x = res_a, y = inv_a, fill = IGR > 0) + 
  geom_raster(alpha = 0.75) + 
  geom_abline(slope = 1, intercept = 0, size = 1, color = "gray") + # indicate neutral line
  geom_vline(data = . %>% group_by(res_a, r, m) %>% # add a vertical line for ESS
               summarize(max_IGR = max(IGR)) %>%
               group_by(r, m) %>% 
               filter(max_IGR == min(max_IGR)), 
             aes(xintercept = res_a),
             linetype = "dashed", size = 0.8) +
  geom_hline(data = . %>% group_by(res_a, r, m) %>% # add a horizontal line for ESS
               summarize(max_IGR = max(IGR)) %>%
               group_by(r, m) %>% 
               filter(max_IGR == min(max_IGR)), 
             aes(yintercept = res_a),
             linetype = "dashed", size = 0.8) +
  geom_vline(data = bif, aes(xintercept = a_c), # add a vertical line for bifurcation point
             linetype = "11", color = "red", size = 1) + 
  annotate(geom = "text", x = 0.5, y = 0.3, label = "\u2013", size = 8) + 
  annotate(geom = "text", x = 0.3, y = 0.5, label = "+", size = 8) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values = c(viridis_pal[10],
                               viridis_pal[4])) + 
  facet_grid(m ~ r, 
             labeller = labeller(m = function(x) paste0("m = ", x),
                                 r = function(x) paste0("r = ", x))) + 
  theme_classic() +
  xlab("Resident active fraction") + ylab("Invader active fraction") + 
  theme(legend.position = "none")

# ggsave(filename = "./figures/si_pip_ricker.png", plot = p, device = "png", 
#       width = 5, height = 5, units = "in", dpi = 400)