library(tidyverse)
library(ggplot2)



dataref = read.table(file = "../data/feas_results_ref.csv")
colnames(dataref) = c("n", "nsols", "nsolsreal", "nsolsfeas")
dataref["model"] = "ref"
data = read.table(file = "../data/feas_results_null.csv")
colnames(data) = c("n", "nsols", "nsolsreal", "nsolsfeas")
data["model"] = "null"
dataopt = read.table(file = "../data/feas_results_opt.csv")
colnames(dataopt) = c("n", "nsols", "nsolsreal", "nsolsfeas")
dataopt["model"] = "opt"

datatot = rbind(dataref, data, dataopt)

datapf = datatot %>% 
  group_by(n, model) %>% 
  summarise(nsolav = mean(nsols),
            nrealav = mean(nsolsreal),
            nposav = mean(nsolsfeas),
            pf = mean(nsolsfeas>0)) 

ggplot(datapf, aes(x = n, y = pf, color = as.factor(model)))+
  geom_point()+
  geom_line()

