
#' ggPlotEmpty
#' @description Version 1.0 20190110
#' @family plot functions
#' @param main plot main
#' @return empty ggplot object
#' @export
ggPlotEmpty <- function(main=NULL){
   ggplot()+
         geom_point(aes(1,1), colour="white") +
         ggtitle(label = main) +
         #geom_text(x=0, y=0, label=label, colour=color) +
         theme(axis.ticks=element_blank(),
         panel.background=element_blank(),
         axis.text.x=element_blank(), axis.text.y=element_blank(),
         axis.title.x=element_blank(), axis.title.y=element_blank())
}



#' Match up a color key to given a series of given unique annotation
#' Add colors if annots are missing
#' @family plot wraps
#' @family color
#' @param color.key
#' @param annots
#' @return color_key vector where all annotatoions have a color
#' @export
colorKeyFix <- function(color.key = NULL, annots =NULL, na.color = "light gray"){
   if(is.null(annots)) stop("character vector of all annotation values must be provided")
   if(is.null(color.key)){
      color_key <- rainbow(length(annots))
   }else{
      u <- match(annots, names(color.key))
      u.missing <- annots[is.na(u)]

      if(all(!is.na(u))){
           color_key <- color.key[u]
      }else{
         temp_key <- rainbow(length(u.missing))
         names(temp_key) <- u.missing
         color_key <- c(color.key[u[!is.na(u)]], temp_key)
         }
   }
   return(color_key)

} # end colorKeyFix


#' Match up a shape key to given a series of given unique annotation
#' Add shapes if annots are missing
#' @family plot wraps
#' @family color
#' @param shape.key
#' @param annots
#' @return shape_key vector where all annotatoions have a color
#' @export
shapeKeyFix <- function(shape.key = NULL, annots = NULL, add.shape = "21"){
   annots = sort(unique(as.character(annots)))
   if(is.null(annots)) stop("character vector of all annotation values must be provided")
   if(is.null(shape.key)){
      shape_key <- rep(add.shape, length(annots))
   }else{
      u <- match(annots, names(shape.key))
      u.missing <- annots[is.na(u)]

      if(all(!is.na(u))){
           shape_key <- shape.key[u]
      }else{
         temp_key <- rep(add.shape,length(u.missing))
         names(temp_key) <- u.missing
         shape_key <- c(shape.key[u[!is.na(u)]], temp_key)
         }
   }
   # if(re.sort) shape_key <- shape_key[order(names(shape_key))]

   return(shape_key)

} # end shapeKeyFix









#' @param palette.sat what saturation levels to use for palette. defaults to c(-2,2)
#' @param pal.no select from list of pre-defoined palettes. 1= 'c(rev(brewer.pal(9,"Blues")[1:7])', "white", brewer.pal(9,"YlOrRd")[1:7]), 2= c(grDevices::rainbow(n=21)), and 3 = c(rev(brewer.pal(11,"Spectral"))
#' @param pal.vec if supplying a palette and not use pre-defined
#' @param na.col, color for NA-values. defaults to 'gray'
#' @param n.col, how menay steps in palette
#' @return a list object with 3 slots:  legend, col and val
#' @export
palette_dl_heatmap2 <- function(x,
                           palette.sat = c(-2,2),
                           pal.no = 1,
                           pal.vec = NULL,
                           na.col = 'gray',
                           n.col = 50
                           ){ # Start function
   # Nov 2016
   requireNamespace("grDevices")
   requireNamespace("RColorBrewer")
   out.list <- list()
   # palettes
   my_pals <- list(
      c(rev(RColorBrewer::brewer.pal(9,"Blues")[1:7]), "white", RColorBrewer::brewer.pal(9,"YlOrRd")[1:7]),
      c(grDevices::rainbow(n=21)),
      c(rev(RColorBrewer::brewer.pal(11,"Spectral")))
   )
   if(is.null(pal.vec)) my.pal <- my_pals[[pal.no]]
   if(!is.null(pal.vec)) my.pal<-pal.vec

            # if x is a matrix - flag and generate a matrix for output
            if(is.matrix(x)){
               matrix.details <- list( ncol=ncol(x),
                                       nrow=nrow(x),
                                       colnames=colnames(x),
                                       rownames=rownames(x)
                                       )
                  x <- as.numeric(x)
                  }else{matrix.details=NULL}
            if(is.integer(x)) x <- as.numeric(x)
            stopifnot(is.numeric(x))
            # b = seq(from=palette.sat[1], to=palette.sat[2], length.out=101)

            my.ramp <- colorRampPalette(my.pal)(n.col)
            b <- seq(from=palette.sat[1], to=palette.sat[2], length.out = n.col+1)
            y <- x
            y[y < palette.sat[1]] <- palette.sat[1]
            y[y > palette.sat[2]] <- palette.sat[2]
            xx <- .bincode(y, breaks=b, right=T, include.lowest = T)
            my.colors <- my.ramp[xx]
            my.colors[is.na(my.colors)] <- na.col
         #create color ramp starting from blue to red
            out.list$col <- my.colors
            out.list$val <- y

            if(!is.null(matrix.details)){
               out.list$col <- matrix(ncol=matrix.details$ncol, data=out.list$col, byrow=F, dimnames = list(matrix.details$rownames, matrix.details$colnames))
               out.list$val <- matrix(ncol=matrix.details$ncol, data=out.list$val, byrow=F, dimnames = list(matrix.details$rownames, matrix.details$colnames))
               }
            out.list$legend <- list(legend=b, col=my.ramp)
            return(out.list)
}


#' Generate ramped color palettes - values and ramps described by vectors with llength n
#' This version (3) is different from palette_dl_heatmap2 since it uses palette numbers defined outside of function
#' @family plot functions
#' @family misc
#' @param x a matrix containg expression values
#' @param palette.sat what saturation levels to use for palette. defaults to c(-2,2)
#' @param my.pal character vector of colors
#' @param pal.vec if supplying a palette and not use pre-defined
#' @param na.col, color for NA-values. defaults to 'gray'
#' @param n.col, how menay steps in palette
#' @return a list object with 3 slots:  legend, col and val
#' @export
palette_gradientRamp <- function(x,
                           palette.sat = c(-2,2),
                           my.pal = dlfoo2::palette_gradients[["blue2white2red"]],
                           na.col = "gray",
                           n.col = 50
                           ){ # Start function
            # Nov 2016
            requireNamespace("grDevices")
            requireNamespace("RColorBrewer")
            out.list <- list()
            # palettes

            # if x is a matrix - flag and generate a matrix for output
            if(is.matrix(x)){
               matrix.details <- list( ncol=ncol(x),
                                       nrow=nrow(x),
                                       colnames=colnames(x),
                                       rownames=rownames(x)
                                       )
                  x <- as.numeric(x)
                  }else{matrix.details=NULL}
            if(is.integer(x)) x <- as.numeric(x)
            stopifnot(is.numeric(x))

            ## If na-color is na & if x is not a matrix, then remove all NA vlaues
            if(is.na(na.col) & !is.matrix(x)){
               x <- x[!is.na(x)]
            }

            # b = seq(from=palette.sat[1], to=palette.sat[2], length.out=101)

            my.ramp <- colorRampPalette(my.pal)(n.col)
            b <- seq(from=palette.sat[1], to=palette.sat[2], length.out = n.col+1)
            y <- x
            y[y < palette.sat[1]] <- palette.sat[1]
            y[y > palette.sat[2]] <- palette.sat[2]
            xx <- .bincode(y, breaks=b, right=T, include.lowest = T)
            my.colors <- my.ramp[xx]
            my.colors[is.na(my.colors)] <- na.col
         #create color ramp starting from blue to red
            out.list$col <- my.colors
            out.list$val <- y

            if(!is.null(matrix.details)){
               out.list$col <- matrix(ncol=matrix.details$ncol, data=out.list$col, byrow=F, dimnames = list(matrix.details$rownames, matrix.details$colnames))
               out.list$val <- matrix(ncol=matrix.details$ncol, data=out.list$val, byrow=F, dimnames = list(matrix.details$rownames, matrix.details$colnames))
               }
            out.list$legend <- list(legend=b, col=my.ramp)

            return(out.list)
            }


#' Small function to extract x-left and x-right values for when plotting rectangles, primarily using circlize package (segments within one sector specified by x-coordinates)
#'
#' @family plots
#' @family rect
#' @family circlize
#' @family circos.rect
#' @param x A vector vith integer values. Each value repesent the size of each reactangle.
#' @param x.offset if to adjust the x-value. Defaults to 0, which means rectangles start at x=0 and end at x=1 for input value 1
#' @return Matrix with The cumulative sum, i.e. x-axis positions  (x-axis positions for x right and x left). A third column with mid-values is provided
#' @export
#'
rectMatrix <- function(x, x.offset= 0 ){
   if(!is.numeric(x)) stop("x must be a numeric vector")
   # if(length(x)<2) stop("no point if only one x-value")
   xx <- cumsum(x)

   y <- matrix(
      data = c(0, xx[-length(xx)], xx),
         #c(1, xr[-length(x)]) + x.offset,
         #xl + x.offset),
      nrow=length(x), byrow=F)

   y <- cbind(y, (y[,2]-y[,1])/2+y[,1])

   return(y)
}


#' Small function to extract x (left and right) values between two groups for plotting in circlize (cicos.link) functions.
#' if multiple links from (and to) segments - then create x values as cumsum so links do not start/end at same positions
#' Used as auxillary funtion to be able to 'zoom' which is not possible in chordDiagram()
#' @family plots
#' @family circlize
#' @family circos.link
#' @param links.df A data frame with connections between two groups.
#' @param scale.x1 if x-axis (no connections should be scaled)
#' @param scale.x2 if y-axis (no connections should be scaled)
#' @param x.offset if to offset x-values
#' @return Data frame with connected x-values to use in links function (x1 and x2 slots)
#' @export
#'
links2sectors <- function(links.df, scale.x1=1, scale.x2=1,  x.offset = 0){
   links.df <- data.frame(links.df)
   links.df$order <- 1:nrow(links.df)
   # fix if na's
   if(any(is.na(links.df$s2))) links.df[is.na(links.df$s2),]$s2 <- "NA"

   x <- split(links.df[,],links.df[,"s1"])
   x1 <- lapply(x, function(y){
         ym <- rectMatrix(y$n*scale.x1, x.offset=x.offset)[,1:2]
         if(is.null(nrow(ym))) ym=matrix(nrow=1, data=ym)
         y$s1x1 <- ym[,1]
         y$s1x2 <- ym[,2]
         return(y)
      })
   xx <- split(links.df[,],links.df[,"s2"])
   x2 <- lapply(xx, function(y){
         ym <- rectMatrix(y$n*scale.x2, x.offset = x.offset)[,1:2]
         if(is.null(nrow(ym))) ym=matrix(nrow=1, data=ym)
         y$s2x1 <- ym[,1]
         y$s2x2 <- ym[,2]
         return(y)
      })
   x1 <- do.call("rbind",x1)
   x2 <- do.call("rbind",x2)
   x1 <- x1[order(x1$order),]
   x2 <- x2[order(x2$order),]
   stopifnot(identical(links.df$s1, x1$s1))
   stopifnot(identical(x2$s1, x1$s1))
   stopifnot(identical(x2$s2, x1$s2))

   x3 <- cbind(x1, x2[,c("s2x1","s2x2")])
   return(x3)
}

#' list of different colors used for heatmaps
#' Object of type List. Vectors with different colors used for plotting
#' @family plot
#' @family misc
#' @export
palette_gradients <- list(
   blue2white2red = c(rev(RColorBrewer::brewer.pal(9,"Blues")[1:7]), "white", RColorBrewer::brewer.pal(9,"YlOrRd")[1:7]),
   blue2white2red_2 = rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4")),
   blue2gray2red = c("#D73027", "#FC8D59", "#FEE090", "#CACACA", "#E0F3F8", "#91BFDB", "#4575B4"),
   blue2white2yellow = c("#FFD800","#FFEDA0","#FFFFFF","#91BFDB","#08519C"),
   blue2black2yellow = c("#FFD800","black","#08519C"),

   green2beige2red = c("#3A5F0B","#FEE090", "#D73027"),
   green2black2red = c("#3A5F0B", "black", "#D73027"),

   red2white2blue = rev(c( rev(RColorBrewer::brewer.pal(9,"Blues")[1:7]), "white", RColorBrewer::brewer.pal(9,"YlOrRd")[1:7])),
   red2white2blue_2 = c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4"),
   rainbow = grDevices::rainbow(n=21),
   spectral = c(rev(RColorBrewer::brewer.pal(11,"Spectral"))),
   spectral_2 = c(colorRamps::matlab.like(9)),
   spectral_dl =  c(rev(RColorBrewer::brewer.pal(9,"GnBu"))[-1], RColorBrewer::brewer.pal(9,"YlOrRd"), "black")

)



#' list of different colors used for heatmaps
#' Object of type List. Vectors with different colors used for plotting
#' @family plot
#' @family misc
#' @export
palette.dl.heatcol <- list(
   blue.white.red = list(c.low=c("#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4"),
                         c.high=c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF")),
   blue.gray.red = list(c.low=c("#CACACA", "#E0F3F8", "#91BFDB", "#4575B4"),
                        c.high=c("#D73027", "#FC8D59", "#FEE090", "#CACACA")),
   blue.white.yellow = list(c.low=c("#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4"),
                            c.high=c("#FFD800", "#FFFFFF")),
   blue.black.red = list(c.low=c("#000000", "#4575B4"),
                         c.high=c("#D73027", "#000000")),
   blue.black.yellow = list(c.low=c("#000000", "#0033FF"),
                            c.high=c("#FFD800", "#000000")),
   green.yellow.red = list(c.low=c("#FEE090", "#3A5F0B"),
                           c.high=c("#FEE090", "#D73027")),
   green.white.red = list(c.low=c("#FFFFFF","#CDDB9D","#F4F7EC","#3A5F0B"),
                          c.high= c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF")),

   white.red = list(c.one=rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF"))),
   white.red.2 = list(c.one=rev(c("#990033","#FFFFFF"))),
   white.blue = list(c.one=c("#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4")),
   white.brown = list(c.one=c("#FFFFFF","#EFE0B9","#E4B04A","#B67721","#643B0F")),
   white.redbrown = list(c.one=c("#FFFFFF","#FDE8D7","#DEADA1","#B85750","#7F171F")),
   white.green = list(c.one=c("#FFFFFF","#CDDB9D","#3A5F0B")),
   white.orange = list(c.one=c("#FFFFFF","#FFDE00","#CC3300")),
   white.black = list(c.one=c("#FFFFFF","#000000"))
)


#' list of different colors used for plotting
#'
#' \code{palette.dl.annot}
#' Object of type List. Vectors with different colors used for plotting
#' @export
palette.dl.annot <- list(
   c("#C0C0C0","#000000"),
   c("#BFCFCC","#3A5F0B"),
   c("#FDE8D7","#7F171F"),
   c("#E0F3F8","#006699"),
   c("#CDDB9D","#3A5F0B"),
   c("#CCCC99","#006699"),
   c("#C3C3E5","#8C489F"),
   c("#EFE0B9","#643B0F"),
   c("#CCCC99","#CC6666","#663333"),
   c("#F4F0CB","#B7C68B","#685642"),
   c("#F4F0CB","#6E7587","#806641"),
   c("#983265","#979735","#006699"),
   c("#154890","#6699FF","#FF6600"),
   c("#EFE0B9","#E4B04A","#643B0F"),
   c("#666633","#CCCC99","#CC6666","#663333"), #15
   c("#B7C68B","#F4F0CB","#B3A580","#685642"),
   c("#FDE8D7","#DEADA1","#B85750","#7F171F"),
   c("#020731","#3862C6","#6E7587","#806641","#F4F0CB"),
   #c("#983265","#986598","#979735","#339798","#006699","#CCCC99"),

   c("#006699","#666633","#CC6666","#663333","#B7C68B","#F4F0CB", "#B3A580","#020731","#3862C6","#6E7587","#806641","#983265","#986598","#339798",
     "#6A3D9A", "#FF7F00", "#E31A1C", "#FB9A99", "#CAB2D6", "#A6CEE3", "#33A02C", "#FDBF6F", "#B2DF8A", "#1F78B4","#7F171F"),
   c("#CA278C","#006699","#FB9A99","#CC6666","#663333","#FFDB00FF","#CC99FF","#CDDB9D","#86942A","#A6CB45","#FFC56C","#B3A580","#986598","#6A3D9A","#CC9900","#0092FFFF","#6E7587","#BB772E","#983265","#3862C6","#339798","#806641","#FF7F00","#FF0000","#AA0114"),

   c(
      as.character(ghibli::ghibli_palette("MarnieMedium2")[c(4,2,5,3,6,7)]),
      as.character(ghibli::ghibli_palette("PonyoMedium")[2:7]),
      as.character(ghibli::ghibli_palette("MarnieMedium1")[3:7]),
      as.character(ghibli::ghibli_palette("LaputaMedium")[3:7]),
      as.character(ghibli::ghibli_palette("MononokeMedium")[4:7]),
      wesanderson::wes_palette("Zissou1")[c(2,3,5,4,1)],
      wesanderson::wes_palette("Darjeeling1")[c(2)],
      wesanderson::wes_palette("Darjeeling2")[c(1,2,4)],
      wesanderson::wes_palette("Chevalier1")[c(1,3,4)],
      wesanderson::wes_palette("Moonrise2")[c(1,3)],
      wesanderson::wes_palette("Moonrise3")[c(1,2,3,5)],
      wesanderson::wes_palette("GrandBudapest1")[c(1,2,3)],
      wesanderson::wes_palette("GrandBudapest2")[c(1,2,3,4)]
   )
)




#' Plot colors
#'
#' \code{plot_colors}
#' Function to view different (hexadecimal) colors given as a character vector
#' @param col.vec, vector with colors
#' @param text.vec, optional. define text to be plotted for each color
#' @return a plot
#' @export

plot_colors <- function(col.vec, text.vec=NULL){
   n<-length(col.vec)
   my.mat<-matrix(ncol=n, data=1)
   colnames(my.mat)<-paste(c(1:n))

   if(!is.null(text.vec)){
      stopifnot(length(col.vec)==length(text.vec))
      #require(colortools)
      #col.vec.comp = sapply(col.vec, complementary, plot=F)
      #text(x=0.25, y=seq(1*2-0.5,n*2-0.5,2), labels=text.vec, col=col.vec.comp[2,])
      colnames(my.mat) = text.vec
   }
   # par(mar=c(12,8,2,8))
   # quartz(width=8.2, height=11.7, title="color test")
   barplot(my.mat, col=col.vec, horiz=T, beside=T, axes=F, las=2)

   if(is.null(text.vec)){
      text(x=0.75, y=seq(1*2-0.5,n*2-0.5,2), labels=col.vec, col="black")
      text(x=0.25, y=seq(1*2-0.5,n*2-0.5,2), labels=col.vec, col="white")
   }

} # end function plot_colors



#' list of different colors used for heatmaps
#'
#' \code{palette.dl.heatcol}
#' Object of type List. Vectors with different colors used for plotting
#' @export
palette.dl.heatcol <- list(
   blue.white.red = list(c.low=c("#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4"),
                         c.high=c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF")),
   blue.gray.red = list(c.low=c("#CACACA", "#E0F3F8", "#91BFDB", "#4575B4"),
                        c.high=c("#D73027", "#FC8D59", "#FEE090", "#CACACA")),
   blue.white.yellow = list(c.low=c("#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4"),
                            c.high=c("#FFD800", "#FFFFFF")),
   blue.black.red = list(c.low=c("#000000", "#4575B4"),
                         c.high=c("#D73027", "#000000")),
   blue.black.yellow = list(c.low=c("#000000", "#0033FF"),
                            c.high=c("#FFD800", "#000000")),
   green.yellow.red = list(c.low=c("#FEE090", "#3A5F0B"),
                           c.high=c("#FEE090", "#D73027")),
   green.white.red = list(c.low=c("#FFFFFF","#CDDB9D","#F4F7EC","#3A5F0B"),
                          c.high= c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF")),

   white.red = list(c.one=rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF"))),
   white.red.2 = list(c.one=rev(c("#990033","#FFFFFF"))),
   white.blue = list(c.one=c("#FFFFFF", "#E0F3F8", "#91BFDB", "#4575B4")),
   white.brown = list(c.one=c("#FFFFFF","#EFE0B9","#E4B04A","#B67721","#643B0F")),
   white.redbrown = list(c.one=c("#FFFFFF","#FDE8D7","#DEADA1","#B85750","#7F171F")),
   white.green = list(c.one=c("#FFFFFF","#CDDB9D","#3A5F0B")),
   white.orange = list(c.one=c("#FFFFFF","#FFDE00","#CC3300")),
   white.black = list(c.one=c("#FFFFFF","#000000"))
)


#' list of different colors used for plotting
#'
#' \code{palette.dl.annot}
#' Object of type List. Vectors with different colors used for plotting
#' @export
palette.dl.annot <- list(
   c("#C0C0C0","#000000"),
   c("#BFCFCC","#3A5F0B"),
   c("#FDE8D7","#7F171F"),
   c("#E0F3F8","#006699"),
   c("#CDDB9D","#3A5F0B"),
   c("#CCCC99","#006699"),
   c("#C3C3E5","#8C489F"),
   c("#EFE0B9","#643B0F"),
   c("#CCCC99","#CC6666","#663333"),
   c("#F4F0CB","#B7C68B","#685642"),
   c("#F4F0CB","#6E7587","#806641"),
   c("#983265","#979735","#006699"),
   c("#154890","#6699FF","#FF6600"),
   c("#EFE0B9","#E4B04A","#643B0F"),
   c("#666633","#CCCC99","#CC6666","#663333"), #15
   c("#B7C68B","#F4F0CB","#B3A580","#685642"),
   c("#FDE8D7","#DEADA1","#B85750","#7F171F"),
   c("#020731","#3862C6","#6E7587","#806641","#F4F0CB"),
   #c("#983265","#986598","#979735","#339798","#006699","#CCCC99"),

   c("#006699","#666633","#CC6666","#663333","#B7C68B","#F4F0CB", "#B3A580","#020731","#3862C6","#6E7587","#806641","#983265","#986598","#339798",
     "#6A3D9A", "#FF7F00", "#E31A1C", "#FB9A99", "#CAB2D6", "#A6CEE3", "#33A02C", "#FDBF6F", "#B2DF8A", "#1F78B4","#7F171F"),
   c("#CA278C","#006699","#FB9A99","#CC6666","#663333","#FFDB00FF","#CC99FF","#CDDB9D","#86942A","#A6CB45","#FFC56C","#B3A580","#986598","#6A3D9A","#CC9900","#0092FFFF","#6E7587","#BB772E","#983265","#3862C6","#339798","#806641","#FF7F00","#FF0000","#AA0114")
)



#' wrapper for easy A4/A5 png plots
#'
#' \code{png_dl}
#' wrapper for easy A4/A5 png plots
#' @param ax,what type of paper size (A4 or A5)
#' @param orientation, prortrait (P) or Landscape (L)
#' @export

png_dl <- function(file.name=NULL, Ax = "A4", orientation = "P"){
   requireNamespace("png")
   if(Ax=="A4" & orientation=="P"){
      width = 596
      height = 842
   }
   if(Ax=="A4" & orientation=="L"){
      width = 842
      height =  596
   }
   if(Ax=="A5" & orientation=="P"){
      width = 421
      height =  596
   }
   if(Ax=="A5" & orientation=="L"){
      width = 596
      height = 421
   }
   png(filename = file.name, width=width, height = height)
} # end png_dl






#' create heat-color matrix from eset or matrix
#'
#' \code{heatcol_creator}
#' Function to create a matrix with hexadecimal colors used for heatmaps
#' @param x, an expressionset or a matrix
#' @param palette.sat, vector with two values that define saturation levels c(min, max)
#' @param pal, define what palette to be used. Defined from object "palette.dl.heatcol"
#' @param na.col, color for na-values
#' @param low.mid.high, if to skew midpoint of color scale. vector with three values (percent), c(1,50,100).
#' @return a list with two slots; col, a matrix with colors; legend, vector with colors mapped to their respecive values
#' @export
heatcol_creator <- function(
   x,
   palette.sat = c(-1,1),
   pal = c('blue.white.red','blue.white.yellow','blue.black.red','blue.black.yellow',
           "white.red","white.blue","white.brown","white.redbrown","white.green","white.black"),
   na.col = 'gray',
   low.mid.high=c(1,50,100)
   #my.pal.list=pal.list
){ # Start function
   # Nov 2015
   requireNamespace("grDevices")

   # if x is a matrix - flag and generate a matrix for output
   #x.in <- x
   if(is.matrix(x)){
      matrix.details <- list( ncol=ncol(x),
                              nrow=nrow(x),
                              colnames=colnames(x),
                              rownames=rownames(x)
      )
   }else{matrix.details=NULL}
   if(is.integer(x)) x <- as.numeric(x)
   stopifnot(is.numeric(x))

   pal <- match.arg(arg = pal, choices = names(palette.dl.heatcol))
   my.pal <- palette.dl.heatcol[pal][[1]]

   #             if(is.null(palette.max) & is.null(palette.min)){
   #                   palette.min <- quantile(x, probs = c(0.1,0.9), na.rm = T)[1]
   #                   palette.max <- quantile(x, probs = c(0.1,0.9), na.rm= T )[2]
   #                   }
   b = seq(from=palette.sat[1], to=palette.sat[2], length.out=101)

   # TWO COLOR PALETTE
   if(length(my.pal)==2){

      c.low <- colorRampPalette(rev(my.pal$c.low))(low.mid.high[2]-low.mid.high[1]+1)
      c.high = colorRampPalette(rev(my.pal$c.high))(low.mid.high[3]-low.mid.high[2]+1)

      c = rep(NA, 100)
      c[low.mid.high[1]:low.mid.high[2]] = c.low
      c[1:low.mid.high[1]] = c.low[1]
      c[low.mid.high[2]:low.mid.high[3]] = c.high
      c[low.mid.high[3]:100] = c.high[length(c.high)]

      foo <- function(xx, b,c){
         u <- which(b >= xx)[1]
         if(is.na(u)) u <- length(b)-1
         if(u==length(b)) u <- length(b)-1
         return(c[u])
      }

      c2 <- sapply(x, foo, b=b, c=c)
      c2[is.na(x)] <- na.col

   } # end if TWO COLORS

   # ONE COLOR PALETTE
   if(length(my.pal)==1){
      #c = colorRampPalette(rev(my.pal$c.one))(low.mid.high[3]-low.mid.high[1]+1)
      c = colorRampPalette(my.pal$c.one)(low.mid.high[3]-low.mid.high[1]+1)
      foo <- function(xx, b,c){
         u <- which(b >= xx)[1]
         if(is.na(u)) u <- length(b)-1
         if(u==length(b)) u <- length(b)-1
         return(c[u])
      }
      c2 <- sapply(x, foo, b=b, c=c)
      c2[is.na(x)] <- na.col

   } # END ONE COLOR

   my.legend <- list(legend=b, col=c)

   # If matrix mode - output is a matrix
   if(!is.null(matrix.details)){
      c2 <- matrix(ncol=matrix.details$ncol, nrow=matrix.details$nrow, data=c2, byrow = F, dimnames = list(matrix.details$rownames, matrix.details$colnames))
   }

   return(list(col=c2,legend=my.legend))
} # END heatcol.pal.creator.dl



#' create heat-color matrix from eset or matrix - alternate version
#'
#' \code{heatcol_creator2}
#' Function to create a matrix with hexadecimal colors used for heatmaps
#' @param x, an expressionset or a matrix
#' @param palette.sat, vector with two values that define saturation levels c(min, max)
#' @param my.pal, a vector withcolors. The palette will be split on these in "equal" parts
#' @param na.col, color for na-values
#' @param n.col, the number of color steps
#' @return a list with two slots; col, a vector/matrix with colors; legend, vector with colors mapped to their respecive values
#' @export
heatcol_creator2 <- function(
   x,
   palette.sat = c(-1,1),
   my.pal = RColorBrewer::brewer.pal(9, "YlOrRd"),
   na.col = 'gray',
   n.col = 50
){ # Start function
   # Nov 2016
   requireNamespace("grDevices")
   requireNamespace("RColorBrewer")
   out.list <- list()

   # if x is a matrix - flag and generate a matrix for output
   if(is.matrix(x)){
      matrix.details <- list( ncol=ncol(x),
                              nrow=nrow(x),
                              colnames=colnames(x),
                              rownames=rownames(x)
      )
      x <- as.numeric(x)
   }else{matrix.details=NULL}
   if(is.integer(x)) x <- as.numeric(x)
   stopifnot(is.numeric(x))
   # b = seq(from=palette.sat[1], to=palette.sat[2], length.out=101)

   my.ramp <- colorRampPalette(my.pal)(n.col)
   b <- seq(from=palette.sat[1], to=palette.sat[2], length.out = n.col+1)
   y <- x
   y[y < palette.sat[1]] <- palette.sat[1]
   y[y > palette.sat[2]] <- palette.sat[2]
   xx <- .bincode(y, breaks=b, right=T, include.lowest = T)
   my.colors <- my.ramp[xx]
   my.colors[is.na(my.colors)] <- na.col
   #create color ramp starting from blue to red
   out.list$col <- my.colors

   if(!is.null(matrix.details)){
      out.list$col <- matrix(ncol=matrix.details$ncol, data=out.list$col, byrow=F, dimnames = list(matrix.details$rownames, matrix.details$colnames))
   }
   out.list$legend <- list(legend=b, col=my.ramp)
   return(out.list)
} # END heatcol_creator2



#' Rbind or cbind color matrixes
#'  Rbind or cbind color matrixes and optinally add a neutral color in between
#'
#' @param my.mat1, a color matrix used for heatmaps
#' @param my.mat2, a color matrix used for heatmaps
#' @param c_or_r, shoud be "c" (cbind) or "r" (rbinds)
#' @param add.color.space, define how many rows/columns if to add a "neutral" color between the matrix
#' @param add.color, what color to add
#' @return a matrix
#' @export
bind.color.matrixes <- function(my.mat1, my.mat2, c_or_r="c", add.color.space=0, add.color="white"){
   if(c_or_r == "c"){
      cat(" \n ... cbinding matrixes")
      stopifnot(nrow(my.mat1)==nrow(my.mat2))
      if(add.color.space>0){
         cat(" \n ... adding add.color.space: ",add.color)
         white.mat <- matrix(nrow=nrow(my.mat1), ncol=add.color.space, data=add.color)
         my.mat1<-cbind(my.mat1, white.mat)
      }
      return(cbind(my.mat1, my.mat2))
   }
   if(c_or_r == "r"){
      cat(" \n ... rbinding matrixes")
      stopifnot(ncol(my.mat1)==ncol(my.mat2))
      if(add.color.space>0){
         cat(" \n ... adding add.color.space: ",add.color)
         white.mat <- matrix(ncol=ncol(my.mat1), nrow=add.color.space, data=add.color)
         my.mat1<-rbind(my.mat1, white.mat)
      }
      return(rbind(my.mat1, my.mat2))
   }
} #end bind color matrixes



#' Widen a color matrix
#' @param my.mat, a matrix
#' @param x, how many times each column shoud be widen
#' @return a matrix
#' @export
widen.color.matrix <- function(my.mat, x=1){
   wide.mat <- matrix(nrow=nrow(my.mat), ncol=x*ncol(my.mat), data=rep(as.vector(t(my.mat)), each=x),byrow=T)
   return(wide.mat)
} # end widen.color.matrix



#' Wrap Text in ggplots.
#' Function to wrap text with newline. Used in ggplots (e.g. enrichmentPlot_bubble)
#' @family plotmisc
#' @family ggplot
#' @param value, a matrix
#' @param width, how many characters are accepted per line
#' @return, text wrapped with newline
#' @export
label_wrap_mod <- function(value, width = 15) {
         sapply(strwrap(as.character(value), width=width, simplify=FALSE),
         paste, collapse="\n")
         }


