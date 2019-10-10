################################################################################################
#
# Examen Estadística Comp
#
################################################################################################


# Pregunta 4 Resuelta por Prof. Teresa
################################################################################################
lambda <- 2.5
calcula_intervalos <- function(n = 60, B = 10000) {
  x <- rpois(n, lambda)
  theta <- exp(-2 * mean(x))
  theta_b <- rerun(B, sample(x, size = n, replace = TRUE)) %>% 
    map_dbl(~exp(-2 * mean(.)))
  bca <- bootstrap::bcanon(x, nboot = B, theta = function(y) exp(-2 * mean(y)), 
                           alpha = c(0.025, 0.975))$confpoints[, 2]
  intervalos <- data_frame(metodo = c("normal", "percent", "BC_a"), 
                           izq = c(theta - 1.96 * sd(theta_b), quantile(theta_b, probs = 0.025), 
                                   bca[1]),
                           der = c(theta + 1.96 * sd(theta_b), quantile(theta_b, probs = 0.975), 
                                   bca[2])
  )
  list(theta = theta, intervalos = intervalos)
}

set.seed(83789173)
n_sims <- 5000
sims_intervalos_60 <- rerun(n_sims, calcula_intervalos()) 
write_rds(sims_intervalos_60, path = "sims_intervalos_60.rds") 
sims_intervalos_60 <- read_rds("data/sims_intervalos_60.rds")
sims_intervalos_60 %>% 
  map_df(~.$intervalos) %>% 
  group_by(metodo) %>%
  summarise(
    falla_izq = 100 * sum(izq > exp(-2 * lambda)) / n_sims, 
    falla_der = 100 * sum(der < exp(-2 * lambda)) / n_sims, 
    cobertura = 100 - falla_izq - falla_der, 
    long_media = mean(der - izq),
    long_min = min(der - izq),
    long_max = max(der - izq)
  )

intervalos_muestra <- sims_intervalos_60 %>% 
  map_df(~.$intervalos) %>% 
  mutate(sim = rep(1:n_sims, each = 3)) %>% 
  filter(sim <= 500) %>% 
  mutate(
    sim_factor = reorder(sim, der - izq), 
    sim = as.numeric(sim_factor)
  )
thetas <- sims_intervalos_60 %>% 
  map_dbl(~.$theta) 

thetas_df <- data_frame(thetas = thetas, sim = 1:n_sims) %>% 
  mutate(
    sim_factor = factor(sim, 
                        levels = levels(intervalos_muestra$sim_factor)), 
    sim = as.numeric(sim_factor)
  ) %>% 
  dplyr::filter(sim <= 500) 

ggplot(intervalos_muestra) +
  geom_hline(yintercept = exp(-2 * 2.5), alpha = 0.6) +
  geom_line(aes(x = sim, y = izq), color = "red", alpha = 0.5) +
  geom_line(aes(x = sim, y = der), color = "red", alpha = 0.5) +
  geom_line(data = thetas_df, aes(x = sim, y = thetas), color = "blue", 
            alpha = 0.5) +
  facet_wrap(~ metodo, ncol = 1)


##############################################################
#     Borrador examen de Estadística Computacional           #
##############################################################

# Pregunta 4 (por mi)
# Los intervalos de ocnfianza tienen problemas

set.seed(8261285)

n <- 60 
lamda_1 <- 2.5
x <- rpois(n, lambda = 2.5)

x_hat <- mean(x)

thetaBoot <- function(){
  x_boot <- rpois(2000, lambda = x_hat)
  lamda_boot <- mean(x_boot)
  exp(-2 * lamda_boot)
}

sims_boot <- rerun(500, thetaBoot()) %>% flatten_dbl()
head(sims_boot, n = 30)  
hist(sims_boot)

# Intervalos de Confianza 

# Normal

li_normal <- round(x_hat - 1.96 * sd(sims_boot), 3)
ls_normal <- round(x_hat + 1.96 * sd(sims_boot), 3)


# Percentil

ls_percentil <- round(quantile(sims_boot, prob = 0.975), 3)
li_percentil <- round(quantile(sims_boot, prob = 0.025), 3)

stringr::str_c(li_normal, ls_normal, sep = ",")
stringr::str_c(li_percentil, ls_percentil, sep = ",")
# No dan datos consistentes entre ellos. Revisar

# Gráfica #

################################################################################################
# Pregunta 5

set.seed(6221285)

# con P_i+1 y i = 0,1,2,...
rBinNegI <- function(n = 10, pe = 0.3){
  U <- runif(1)
  i <- 0
  p <- pe^n
  P <- p
  while (U >= P) {
    p <- p * (1 - pe) * (i + n / i + 1)
    P <- P + p
    i <- i + 1
  }
  i
}
sims_binNeg <- rerun(10000, rBinNegI()) %>% flatten_dbl()
qplot(sims_binNeg, binwidth = 1)


# Prueba con P_i y i = 1,2,3,...
# Son los resultados más consistentes pero no sé si esté bien
rBinNegII <- function(n = 10, pe = 0.3){
  U <- runif(1)
  i <- 1
  p <- pe^n
  P <- p
  while (U >= P) {
    p <- p * (1 - pe) * ((i + n - 1) / (i))
    P <- P + p
    i <- i + 1
  }
  i
}
sims_binNeg_II <- rerun(10000, rBinNegII()) %>% flatten_dbl()
head(sims_binNeg_II, n = 5)

qplot(sims_binNeg_II, binwidth = 1)

################################################
# Prueba de gráfica

negbin <- rnbinom(10000, 10, .3)
qplot(negbin,binwidth = 1)

x_nbin <- seq(0,80)
y_nbin <- dnbinom(x_nbin, size = 10, prob = .3)
plot(x_nbin, y_nbin, type = "l")


################################################################################################
# Reproducción pregunta poisson de León
################################################################################################
rpoisI <- function(lambda = 1){
  U <- runif(1)
  i <- 0
  p <- exp(-lambda)
  P <- p
  while(U >= P){
    p <- lambda * p / (i + 1)
    P <- P + p
    i <- i + 1
  }
  i
}
sims_pois <- rerun(200, rpoisI()) %>% flatten_dbl()
qplot(sims_pois, binwidth = 1)


