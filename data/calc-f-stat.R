
rm(list=ls())
setwd("~/Documents/UCSF_year1/Zaitlen-rotation1/Zaitlen_lab/data/")

m = read.csv("file3.csv", na.strings=c("","NA"))
# write.csv(m, "file3_with_NAs.csv")

key = read.csv("reference_key.csv")
m_relevant = m[, which(names(m) %in% key$Name)]
m_relevant[1:75] <- lapply(m_relevant[1:75], as.character)
m_relevant[1:75] <- lapply(m_relevant[1:75], as.numeric)

m_relevant_NA <- replace(m_relevant, m_relevant == -1, NA)

m_cov = cov(m_relevant_NA, use = "complete.obs")

muscle = c("GSM1179509", "GSM1179510", "GSM1179511", "GSM1179512", "GSM1179513")
brain = c("GSM992014",	"GSM992015",	"GSM992016",	"GSM992017",	"GSM992018", "GSM992019",	"GSM992020",	"GSM992021","GSM992022",	"GSM992023",	"GSM992024",	"GSM992025",	
          "GSM992026",		"GSM992028",	"GSM992030",	"GSM992031",	"GSM992032",	"methylation_level_E050_RO_01549",	"methylation_level_E053",	"methylation_level_E054",
          "methylation_level_E054")

muscle_subset <- m_relevant_NA[muscle]
brain_subset <- m_relevant_NA[brain]
var.test(muscle_subset, brain_subset)

for (row in 1:nrow(m_relevant_NA)) { 
  var.test(as.vector(muscle_subset[row,]), as.vector(brain_subset[row,]))
  
}

typeof(as.vector(brain_subset[row,]))
       