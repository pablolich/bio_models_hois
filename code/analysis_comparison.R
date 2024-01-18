library(tidyverse)

data = read.table("comparedists_greedy_brute.csv")
colnames(data) = c("n", "sim", "flips", "dist", "best")

toplot = data %>% 
  group_by(n) %>% 
  summarize(avflips = mean(flips)) %>% 
  mutate(total = n*(n+1)) %>% 
  pivot_longer(n)
  ggplot(aes(x=n, y=avflips))+
  geom_bar()

toplot
