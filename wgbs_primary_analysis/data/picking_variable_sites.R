rm(list=ls())
setwd("Desktop/zaitlen_lab_desktop/")
top_var = read.table("training_with_var.txt")


all_greatest = subset(top_var, top_var$V18==0.2489)
random_10000_greatest = all_greatest[sample(nrow(all_greatest), 10000), ]

just_top= top_var[1:10000, ]
random_all = top_var[sample(nrow(top_var), 10000), ]
q
write.table(greatest, file="top_var_new_roadmap.txt")
