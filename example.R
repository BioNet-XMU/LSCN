# Date: 20251102
# Author: Fanjing Guo
# Title: An example to apply the LSCN

library(here)
source(here("LSCN_Function", "Funcs_LSCN_Algorithm.R"))

# Read the input data
data <- read.csv("example_dataset1.csv", header = T)
data <- as.matrix(data)

# Normalize data using LSCN
Norm.Result <- LSCN(data, ita0 = 0.4, max_ite = 50, print.txt = T)
data.normalized <- Norm.Result$dat.norm
norm.factor <- Norm.Result$NORM.COEF
conv <- Norm.Result$conv.vals
runtime <- Norm.Result$RunTime

# View the normalized data matrix
View(data.normalized)

# View the normalization factors
print(norm.factor)

# Plot the convergence curve
library(ggplot2)
conv.df <- data.frame(iteration = c(1: length(conv)), 
                      value = conv)

ggplot(conv.df) + geom_line(aes(x = iteration, y = value)) + 
  geom_point(aes(x = iteration, y = value), shape = 21) + 
  theme_bw() + 
  xlab('Iteration') + 
  ylab('Value') + 
  theme(axis.title = element_text(size = 25, family = 'sans', face = 'bold'), 
        axis.text = element_text(size = 22, family = 'sans', face = "bold", color = 'black'), 
        legend.position = 'none', 
        panel.border = element_rect(linewidth = 2), 
        panel.grid = element_blank())

# Print the run time
print(runtime)
