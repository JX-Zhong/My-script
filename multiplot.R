multiplot <- function(...,plotlist = NULL,cols = 1,layout = NULL) {
  # layout : a vector of the layout. 
  #          Example :
  #          c(111211131114,3) ->            [,1] [,2] [,3] [,4]
  #                                    [1,]    1    1    1    2
  #                                    [2,]    1    1    1    3
  #                                    [3,]    1    1    1    4
  library(grid)
  library(dplyr)
  
  plots <- c(list(...),plotlist)
  
  numPlots <- length(plots)
  
  if(is.null(layout)) {
    layout_mtx <- matrix(seq(1,cols * ceiling(numPlots/cols)),
                     ncol = cols, byrow = T)
  } else {
    layout_mtx <- substring(layout[1],1:nchar(layout[1]),1:nchar(layout[1])) %>%
                  as.numeric() %>% 
                  matrix(nrow = layout[2],byrow = T)
    
  }
    
  if(numPlots == 1 ) {
    print(plots[[1]])
  } else {

    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout_mtx), 
                                               ncol(layout_mtx))))
    
    for(i in 1:numPlots) {
       matchidx <- as.data.frame(which(layout_mtx == i, arr.ind = TRUE))
       
       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                       layout.pos.col = matchidx$col))
    }
  }
}
