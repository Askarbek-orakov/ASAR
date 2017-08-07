mg.key <- function(fname,
                   names,
                   low = "green",
                   mid = "gray",
                   high = "red",
                   node.size,
                   bins = 10,
                   res = 300,
                   key.pos = "topright",
                   off.sets = c(x = 0, y = 0),
                   align = "n",
                   cex = 1,
                   lwd = 1,
                   tstep=0.1,
                   rstep=0.1) {
  col.pos <-
    function(bins = 10,
             graph.size,
             node.size,
             size.by.graph = TRUE,
             key.pos = "topright",
             off.sets = c(x = 0, y = 0),
             align = "n") {
      width = graph.size[1]
      height = graph.size[2]
      if (size.by.graph == T) {
        xs = width / 80
        ys = height / 40
      } else if (!missing(node.sizes)) {
        xs = node.size[1] / bins
        ys = node.size[2]
      } else{
        message("Note: ",
                "color key not plotted, node.size is needed\n when size.by.graph=FALSE!")
        return(off.sets)
      }

      if (align == "x") {
        off.sets['x'] = 4 * xs
        off.sets['y'] = off.sets['y'] + 3 * ys
      }
      if (align == "y")
        off.sets = off.sets + c(x = 3 * xs, y = 0)
      if (align == "n")
        off.sets = off.sets + c(x = 2 * xs, y = 2 * ys)
      if (length(grep('right', key.pos)) == 1) {
        off.sets['x'] = off.sets['x'] + bins * xs
        x = width - off.sets['x']
      } else {
        x = off.sets['x']
        off.sets['x'] = off.sets['x'] + bins * xs
      }
      if (length(grep('top', key.pos)) == 1)
        y = height - off.sets['y']
      else
        y = off.sets['y']

      return(list(
        off = off.sets,
        x = x,
        y = y,
        xs = xs,
        ys = ys
      ))
    }
# the code was taken from wordcloud package https://github.com/ifellows/wordcloud
# it was slightly modified for the purposes of KEGG legend, so it won't be used from the 
# package above
  overlap <- function(x1, y1, sw1, sh1,boxes) {
    s <- 0
    if (length(boxes) == 0)
      return(FALSE)
    for (i in 1:length(boxes)) {
      bnds <- boxes[[i]]
      x2 <- bnds[1]
      y2 <- bnds[2]
      sw2 <- bnds[3]
      sh2 <- bnds[4]
      if (x1 < x2)
        overlap <- x1 + sw1 > x2-s
      else
        overlap <- x2 + sw2 > x1-s

      if (y1 < y2)
        overlap <- overlap && (y1 + sh1 > y2-s)
      else
        overlap <- overlap && (y2 + sh2 > y1-s)
      if(overlap){
        return(TRUE)
      }
    }
    FALSE
  }

  # the code was taken from wordcloud package https://github.com/ifellows/wordcloud
  # it was slightly modified for the purposes of KEGG legend, so it won't be used from the 
  # package above
  wordlayout <- function(x, y, words, cex=1, rotate90 = FALSE,
                         xlim=c(-Inf,Inf), ylim=c(-Inf,Inf), tstep=.1, rstep=.1, ...){
    tails <- "g|j|p|q|y"
    n <- length(words)
    sdx <- sd(x,na.rm=TRUE)
    sdy <- sd(y,na.rm=TRUE)
    if(sdx==0)
      sdx <- 1
    if(sdy==0)
      sdy <- 1
    if(length(cex)==1)
      cex <- rep(cex,n)
    if(length(rotate90)==1)
      rotate90 <- rep(rotate90,n)

    dx<-min(x)-min(xlim)
    dy<-max(y)-max(ylim)
    boxes <- list()
    for(i in 1:length(words)){
      rotWord <- rotate90[i]
      r <-0
      theta <- ((i-1)*2*pi/length(words)-pi)#runif(1,0,2*pi)
      x1 <- xo <- x[i]-dx
      y1 <- yo <- y[i]-dy
      wid <- strwidth(words[i],cex=cex[i],...)
      ht <- strheight(words[i],cex=cex[i],...)
      #mind your ps and qs
      if(grepl(tails,words[i]))
        ht <- ht + ht*.2
      if(rotWord){
        tmp <- ht
        ht <- wid
        wid <- tmp
      }
      isOverlaped <- TRUE
      while(isOverlaped){
        if(!overlap(x1-.5*wid,y1-.5*ht,wid,ht,boxes) &&
           x1-.5*wid>xlim[1] && y1-.5*ht>ylim[1] &&
           x1+.5*wid<xlim[2] && y1+.5*ht<ylim[2]){
          boxes[[length(boxes)+1]] <- c(x1-.5*wid,y1-.5*ht,wid,ht)
          isOverlaped <- FALSE
        }else{
          theta <- theta+tstep
          r <- r + rstep*tstep/(2*pi)
          x1 <- xo+sdx*r*cos(theta)
          y1 <- yo+sdy*r*sin(theta)
        }
      }
    }
    result <- do.call(rbind,boxes)
    colnames(result) <- c("x","y","width","ht")
    rownames(result) <- words
    result
  }


    wp <- function (x,
                  y,
                  words,
                  cex = 0.5,
                  lwd = 0.5,tstep=tstep,rstep=rstep,
                  ...)
  {
    lay <- wordlayout(x, y, words,tstep=tstep,rstep=rstep, cex, ...)
      for (i in 1:length(x)) {
        xl <- lay[i, 1]
        yl <- lay[i, 2]
        w <- lay[i, 3]
        h <- lay[i, 4]
        if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] >
            yl + h) {
          nx <- xl + 0.5 * w
          ny <- yl + h
          lines(c(x[i], nx), c(y[i], ny), col = "grey",lwd=lwd)
        }
    }
    text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4],
         words, cex = cex, ...)
  }


  img <- png::readPNG(fname)
  width <- ncol(img)
  height <- nrow(img)
  of1 <- col.pos(
    bins = bins,
    graph.size = c(width, height),
    node.size = c(node.size[1] * 3, node.size[2]),
    align = 'n'
  )
  lbins <- length(names)
  of2 <- col.pos(
    bins = lbins,
    graph.size = c(width, height),
    node.size = node.size,
    align = 'y',
    off.sets = of1$off
  )
  #names<-gsub(' +','\n',names)
  img.file <- paste0(gsub('.png$', '', fname), '.leg.png')
  png(img.file,
      width = width,
      height = height,
      res = res)
  op = par(mar = c(0, 0, 0, 0))
  plot(
    c(0, width),
    c(0, height),
    type = "n",
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  rasterImage(img, 0, 0, width, height, interpolate = F)
  x <- of2$x
  y <- of2$y
  xs <- of2$xs
  ys <- of2$ys
  ckx = seq(x, x + lbins * xs, length = lbins + 1)
  cky = c(y, y + ys)
  cols = gplots::colorpanel(length(names),
                            low = low,
                            mid = mid,
                            high = high)
  data.cuts = seq_along(names)
  image(
    x = ckx,
    y = cky,
    z = cbind(data.cuts),
    col = cols,
    axes = FALSE,
    add = T
  )
  names<-gsub(' +','\n',names)
  totstr<-paste(names,collapse = ' ')
  wp(
    x = ckx[-1] - xs / 2,
    y = rep(y+ys/2, length(names)),
    words = names,
    new = FALSE,
    cex = cex,
    ylim = c(-4*max(strheight(names),cex = cex), -ys/2)+y,
    xlim = c(-1.2*sum(strwidth(names,cex = cex)),xs)+width -  of1$off['x'],
    lwd=lwd,tstep=tstep,rstep=rstep
  )

  # xy<-maptools::pointLabel(x=ckx[-1]-xs/2,
  #                          y=rep(y-ys, length(names)),
  #                          labels =names, cex=cex,doPlot = T,method = 'GA',allowSmallOverlap = FALSE)
  # for(i in 1:lbins){
  #   lines(c(xy$x[i],ckx[i]+xs/2),c(xy$y[i],y),lwd=lwd)
  # }
  # text(x=seq(x,x+lbins*xs,
  #            length=length(names)),
  #      y=rep(y-ys,
  #            length(names)),
  #      label=names, cex=cex)
  par(op)
  dev.off()

}
