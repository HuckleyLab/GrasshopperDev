#Functions/code that I commonly use, in nearly every script
#Author: Rory S. Telemeco
#Last updated: 06/20/2016


#ggplot2 theme:
rory_theme <- theme_classic(base_size = 28) + 
  theme(axis.title.y=element_text(vjust=1.5), axis.title.x=element_text(vjust=0.2)) + #adjust axis title position
  theme(plot.margin = unit(c(.3,.3,.6,.6), "cm"), line = element_line(size = 1.25)) + #adjust plot margins and line element size
  theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black")) + #draw x and y axes
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) + #put margins around axis labels so that nothing overlaps
  theme(axis.ticks.length =unit(-0.3, "cm")) + # move tickmarks inside the axes
  theme(panel.margin = unit(2, units = "lines")) + #spread out facets
  theme(strip.background = element_blank()) #remove border from facet labels

#rory_theme_OLD <- theme_classic(base_size = 28) + theme(axis.title.y=element_text(vjust=1.5), axis.title.x=element_text(vjust=0.2)) + #adjust axis and text size and position
#  theme(axis.line=element_line(size = 1.25), plot.margin = unit(c(.3,.3,.6,.6), "cm"), axis.ticks=element_line(size = 1.25)) + #adjust plot margins and line element size
#  theme(axis.ticks.length =unit(-0.3, "cm"), axis.ticks.margin = unit(0.7, "cm")) + #move tick marks to the inside
#  theme(panel.margin = unit(2, units = "lines")) #+

#function for graphical model validation
Validate <- function(Mod, main) {
  if(missing(main)) main <- deparse(substitute(Mod))
  E<-resid(Mod)
  Fit <- fitted(Mod)
  op <- par(mfrow =c(1,2))
  plot(x=Fit, y=E,
       xlab = "Fitted values", ylab = "Residuals",
       main = paste("Residuals vs Fitted Values from", main))
  hist(E, nclass = 15, 
       xlab = "Residuals",
       main = paste("Histrogram of residuals from", main))
  par(op) 
}


#Function for creating a csv file with pairwise test results from lsmeans
#Inputs are a standard pairwise contrast object with tukey corrected pvalues (PT), and same object but original code had adjust = "none" (PN), and a filename (Fname)
#Creates a csv file in the working directory and returns data fame with the values
LSwrite <- function (PT, PN, Fname = "PairwiseResults.csv") {
  PT_db <- data.frame(summary(PT[[2]]))
  PT_db$Pcorrection <- "Tukey"
  PN_db <- data.frame(summary(PN[[2]]))
  PN_db$Pcorrection <- "None"
  Pall <- rbind(PT_db, PN_db)
  write.csv(Pall, file = Fname)
  return(Pall)
}

std.err <- function(x) sd(x)/sqrt(length(x))
