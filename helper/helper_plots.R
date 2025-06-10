axis_title <- 16
title_size <- 18
axis_text_size <- 14
legend_text <- 14
line_size <- 1



theme_single <- function(){ 
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.spacing.x = unit(1, "lines"),
      
      #text elements
      plot.title = element_text(size = title_size),
      axis.text = element_text(size = axis_text_size - 2),
      axis.title = element_text(size = axis_title),
      legend.text = element_text(size = legend_text),
      strip.text = element_text(size = axis_title - 5, 
                                margin = margin(1.5, 1.5, 2, 1.5)),
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      legend.position = "none"
    )
}

theme_single_legend <- function(){ 
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.spacing.x = unit(1, "lines"),
      
      #text elements
      plot.title = element_text(size = title_size),
      axis.text = element_text(size = axis_text_size - 2),
      axis.title = element_text(size = axis_title),
      legend.text = element_text(size = legend_text),
      strip.text = element_text(size = axis_title - 5, 
                                margin = margin(1.5, 1.5, 2, 1.5)),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom"
    )
}


theme_single_y <- function(){ 
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.spacing.x = unit(1, "lines"),
      
      #text elements
      plot.title = element_text(size = title_size),
      axis.text = element_text(size = axis_text_size - 2),
      axis.title = element_text(size = axis_title),
      legend.text = element_text(size = legend_text),
      strip.text = element_text(size = axis_title - 5, 
                                margin = margin(1.5, 1.5, 2, 1.5)),
      legend.position = "none"
    )
}


theme_single_grid <- function(){ 
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.spacing.x = unit(1, "lines"),
      
      #text elements
      plot.title = element_text(size = title_size),
      axis.text = element_text(size = axis_text_size - 4),
      axis.title = element_blank(),
      legend.text = element_text(size = legend_text),
      strip.text = element_text(size = axis_title - 5, 
                                margin = margin(1.5, 1.5, 2, 1.5)),
      legend.position = "none"
    )
}
