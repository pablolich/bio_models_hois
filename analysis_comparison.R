library(tidyverse)

data = read.table("../compare_greedy_brute.csv")
colnames(data) = c("n", "sim", "same", "flips")

toplot = data %>% 
  group_by(n) %>% 
  summarise(avflips = mean(flips)) %>% 
  ggplot(aes(x=n, y=avflips))+
  geom_point()

toplot
