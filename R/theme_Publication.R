
#' Plot variable dependence using an oblique random survival forest
#' @param base_size how big to make the text


theme_Publication <- function(base_size=16) {
  
  (theme_foundation(base_size=base_size)+ 
     theme(plot.title = element_text(face = "bold",
                                     size = rel(1), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.border = element_rect(colour = 'black'),
           legend.key.size = unit(3,"line"),
           legend.key.height = unit(3,"line"),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size=rel(1)),
           axis.line = element_blank(), 
           axis.ticks = element_line(),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.spacing = unit(2, "lines"),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.title = element_text(face="italic"),
           legend.text = element_text(size=rel(1)),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="black",fill="#f0f0f0"),
           strip.text = element_text(size=rel(1), face="bold")
     ))
  
}