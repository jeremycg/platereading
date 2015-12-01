#a hideous hack to allow me to use ggplot2  and dplyr in this package, and not destroy R CMD check
#see http://stackoverflow.com/a/12429344/3760920 for details
utils::globalVariables(c("well","plate","strain","temperature"))

#' The Buchanan model for three phase growth
#'
#' The model as described by  Buchanan, R. L., Whiting, R. C. & Damert,
#' W. C. 1997 When is simple good enough: A comparison of the Gompertz, Baranyi,
#' and three-phase linear models for fitting bacterial growth curves.
#' Food Microbiol. 14, 313-326.
#' @usage buchworm
#' @format an equation, with od modelled as a response of od0, lag, mumax and odmax
buchworm <- as.formula(od ~ od0 + (time >= lag) * (time <= (lag + (odmax - od0) *
                                                              log(10)/mumax)) * mumax * (time - lag)/log(10) + (time >= lag) *
                         (time > (lag + (odmax - od0) * log(10)/mumax)) * (odmax - od0))


#' Read in a single plate directory of data
#'
#' @param x A Directory containing one plate run
#' @return A dataframe containing the data sorted by well and time
#' @importFrom reshape2 melt
#' @export
readonedir <- function(x){
  x3 <- lapply(paste0(x, "/", list.files(x)),
    function(file){
       x1 <- melt(read.csv(file, check.names = FALSE)[1:13], id.vars = "")
       names(x1) <- c("column", "row", "od")
       x1$well <- paste0(x1$column, x1$row)
       x1$time <- gsub("^.*/|.csv$", "", file)
       x1
       }
  )
  x3 <- do.call(rbind, x3)
  x3$time <- as.numeric(x3$time)
  x3$plate <- x
  return(x3)
}


#' Fit a single directory of data (output by readonedir) with the buch function
#'
#' @param x A Dataframe produced by readonedir
#' @param y an indication of whether to omit t=0 or not, defaults to 1, remove
#' @return A dataframe containing the fitted data

fitbuch<-function(x,y=1, lag1 = 35, mumax1 = 0.025, od01 = 0.25, odmax1 = 0.95){
  x$od<-1-x$od
  if(y==1){x=x[x$time!=0,]}
  z<-NULL
  try(z<-nls(buchworm,x,list(lag = lag1, mumax = mumax1, od0 = od01, odmax = odmax1),control=nls.control(minFactor=1/4096,warnOnly=T)))
  if(is.null(z)){
    z3<-data.frame(t(rep(NA,8)))
    names(z3)=c("lag","mumax","od0","odmax","residual","well","plate","run")
    return(z3)
  }
  z3=data.frame(t(matrix(coef(z))),sum(resid(z)^2)/nrow(x),x$well[1],
                x$plate[1],as.numeric(tail(strsplit(as.character(x$plate[1]), split=' ', fixed=TRUE)[[1]],1)))
  names(z3)=c("lag","mumax","od0","odmax","residual","well","plate","run")
  return(z3)
}

#' Fit multiple directories of data (output by readonedir) with the buch function
#'
#' @param x A Dataframe produced by outputs of readonedir
#' @return A dataframe containing the fitted data
#' @importFrom plyr ddply
#' @importFrom plyr "."
plyrfit<-function(x, lag1 = 35, mumax1 = 0.025, od01 = 0.25, odmax1 = 0.95){
  ddply(x,.(plate,well),function(df){fitbuch(df, lag = lag1, mumax = mumax1, od0 = od01, odmax = odmax1)})
}

#' Fit multiple directories of data (output by readonedir) with the buch function
#'
#' @param x A Dataframe produced by outputs of readonedir
#' @return A dataframe containing the fitted data
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by
#' @importFrom dplyr do
#' @importFrom plyr "."
dplyrfit<-function(x, lag1 = 35, mumax1 = 0.025, od01 = 0.25, odmax1 = 0.95){x %>%
                        group_by(plate, well) %>%
                        do(fitbuch(., lag = lag1, mumax = mumax1, od0 = od01, odmax = odmax1))
}



#' Read in a directory containing multiple plates of data
#'
#' @export
#' @param x A Directory containing all data and required outputfiles
#' @return A dataframe containing the data fit and sorted by well
#' @importFrom plyr rbind.fill
looper <- function(x, lag1 = 35, mumax1 = 0.025, od01 = 0.25, odmax1 = 0.95){
  startingdir = getwd()
  z3 = data.frame()
  setwd(x)
  dirlist = list.files()
  for(zz in 1:(length(dirlist)-1)){
    z3 = rbind.fill(z3,dplyrfit(readonedir(dirlist[zz]), lag = lag1, mumax = mumax1, od0 = od01, odmax = odmax1))
  }
  setwd(startingdir)
  return(z3)
}


#' Name output from a directory containing multiple plates of data
#'
#' @param x A dataframe containing the output of the looper function
#' @export
#' @param y A file containing all names of wells in plates
#' @return A dataframe containing the data fit and sorted by well

namer<-function(x,y){
  strainlist=read.csv(y)
  x$column=as.numeric(substring(as.character(x$well), 2, nchar(as.character(x$well))))
  x$row=substring(as.character(x$well), 1,1)
  z2=merge(x, strainlist, by=c("run","column"))
  return(z2)
}

#put it all together!
#setwd(choose.dir())
#if(file.exists("grouped.csv")){file.remove("grouped.csv")}
#if(file.exists("outputfits.csv")){file.remove("outputfits.csv")}
#x=namer(looper(getwd()),"strainlist.csv")


#' Filter a data frame by residuals, and combine passing values into means and sds
#'
#' @param x A dataframe containing the named output of the looper function
#' @param y A residual cutoff, defaults to 0.01
#' @return A dataframe containing the data filtered and summarised
#' @importFrom plyr ddply
#' @importFrom plyr "."
#' @export
compileall<-function(x,y=0.01){
  x=x[which(x$residual<y),]
  z=ddply(x,.(strain,temperature),function(df){c(mean(df$lag),mean(df$mumax),mean(df$od0),mean(df$odmax),
                                                 sd(df$lag),sd(df$mumax),sd(df$od0),sd(df$odmax),length(df$lag))})
  names(z)=c("strain","temperature","lag","mumax","od0","odmax","sdlag","sdmumax","sdod0","sdodmax","N")
  return(z)
}




#' Plot all representatives of a strain with their fits and residuals
#'
#' @param x A directory containing multiple plate runs
#' @param y A strain to plot out
#' @param z The file containing strain names
#' @return plots of each well of the required strain
#' @importFrom plyr d_ply
#' @importFrom plyr "."
#' @export
plotlots<-function(x, y, z, lag1 = 35, mumax1 = 0.025, od01 = 0.25, odmax1 = 0.95){ #function of directory, strain, strainlist
  startdir = getwd()
  setwd(x)
  strainlist = read.csv(z)
  strainlist = strainlist[which(strainlist$strain == y), ]
  if(length(strainlist[ ,1]) == 0){return("strain not found")}
  for(i in 1:length(strainlist[ ,1])){
    singlehold = readonedir(paste("Plate", sprintf("%02d", strainlist[i, ]$run)))
    singlehold = singlehold[which(substring(singlehold$well, 2, nchar(x)) == strainlist[i, ]$column), ]
    singlehold$temperature = strainlist[i, ]$temperature
    d_ply(singlehold, .(well), plotter, lag = lag1, mumax = mumax1, od0 = od01, odmax = odmax1)
  }
  setwd(startdir)
}

#' Plot a single well with its fit
#'
#' @param x dataframe containing runs
#' @param well the well to choose
#' @export
#' @return plots of the well
plotter<-function(x,well, lag1 = 35, mumax1 = 0.025, od01 = 0.25, odmax1 = 0.95){
  x=x[x$time!=0,]
  x$od=1-x$od
  workingfit<-nls(buchworm, x,list(lag = lag1, mumax = mumax1, od0 = od01, odmax = odmax1),control=nls.control(minFactor=1/4096,warnOnly=T))
  plot(x$od~x$time,ylab="od",xlab="hours",main=c(paste("Plate:", x$plate[1],
                                                       "Temp:", x$temperature[1],"Well",x$well[1],"Residual:",round(sum(resid(workingfit)^2)/nrow(x),4))))
  points(x$time,predict(workingfit,list(time=x$time)),pch=3,col=2)
}

#' A shiny app for interactive production of plots and outputs
#'
#' @param directory the directory containing data and strain names
#' @return a shiny app
#' @export
#' @importFrom gplots plotmeans
#' @importFrom shiny shinyApp pageWithSidebar headerPanel sidebarPanel sliderInput
#' @importFrom shiny selectInput mainPanel tabsetPanel tabPanel renderPlot reactive
#' @importFrom shiny checkboxInput plotOutput tableOutput renderTable actionButton
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot geom_line geom_point facet_grid
plateshiny <- function(directory) {
  startingdir<-getwd()
  setwd(directory)
  strainlist<-read.csv("strainlist.csv")
  platelist <- list.files()
  platelist <- platelist[grepl("^[Pp]late", platelist)]
  stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
  if(!file.exists("outputfits.csv")){
    write.csv(namer(looper(getwd(), lag1 = 35, mumax1 = 0.025, od01 = 0.25, odmax1 = 0.95),"strainlist.csv"),file="outputfits.csv",row.names=F)}
  data<-read.csv("outputfits.csv")
  data$factortemp=as.factor(data$temperature)
  shinyApp(
    ui = pageWithSidebar(
      headerPanel("Plate Analysis"),
      sidebarPanel(
        sliderInput("cutoff","cutoff of residuals:",min = 0.0,
                    max = 1.0, value = 0.5, step = 0.001, ticks = T),
        selectInput("variable", "Choose a variable:",
                    choices = c("mumax", "lag", "od0", "odmax")),
        selectInput("query", "Strain:", levels(strainlist$strain)),
        sliderInput("initialmumax","initialmumax",min = -1.0,
                    max = 1.0,value = 0.025,step=0.001,ticks=T),
        sliderInput("initiallag","initial lag:",min = -20.0,
                    max = 150,value = 35,step=0.5,ticks=T),
        sliderInput("initialod0","initialod0:",min = -1.0,
                                max = 1.0,value = 0,step=0.001,ticks=T),
        sliderInput("initialodmax","initialodmax:",min = 0.0,
                                max = 2.0,value = 0.95,step=0.001,ticks=T),
        actionButton("do", "Click to remake output.csv")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("BoxPlots", checkboxInput("ordered1",
                                            "Ordered?", value = FALSE), plotOutput("Plot1", height = "100%")),
          tabPanel("Table", tableOutput("table1"), tableOutput("table2")),
          tabPanel("Mean and sd", checkboxInput("ordered2",
                                               "Ordered?", value = FALSE), plotOutput("Plot3", height= "100%")),
          tabPanel("Fitted Plots", plotOutput("Plotfitted", height = "100%")),
          tabPanel("Plate Plot", plotOutput("Plotstrain", height= "100%"))
        )
      )
    ),
    server = function(input, output) {
      on.exit(setwd(startingdir))
      output$Plot1 <- renderPlot({
        par(mfrow=c(length(levels(data$factortemp)),1))
        for(i in 1:length(levels(data$factortemp))){
          data2<-data[which(data$factortemp==levels(data$factortemp)[i]),]
          dataclean <- data2[data2$residual<input$cutoff,]
          if(input$ordered1==T){
            dataclean$strain<-factor(dataclean$strain,
                                     names(sort(rank(tapply(subset(dataclean,select=input$variable)[,1],
                                                            dataclean$strain,median),ties.method="first")))
            )
          }
          plot(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,
               xlab="",ylab=as.character(input$variable),
               main=c(levels(data$factortemp)[i]," degrees"),las=2)
        }
      }, height = 1000, width = 1500)
      output$table1 <- renderTable({
        compileall(data,input$cutoff)
      },digits=5)
      output$Plot3 <- renderPlot({
        par(mfrow=c(length(levels(data$factortemp)),1))
        for(i in 1:length(levels(data$factortemp))){
          data2<-data[which(data$factortemp==levels(data$factortemp)[i]),]
          dataclean <- data2[data2$residual<input$cutoff,]
          if(input$ordered2==T){
            dataclean$strain<-factor(dataclean$strain,
                                     names(sort(rank(tapply(subset(dataclean,select=input$variable)[,1],
                                                            dataclean$strain,mean),ties.method="first"))))
          }
          plotmeans(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,
                    n.label=F,barcol="black",connect=F,
                    main=c(levels(data$factortemp)[i]," degrees"),xlab="",
                    ylab=as.character(input$variable),p=0.68,par(las =2)
          )
        }
      }, height = 1000, width = 1500)
      output$Plotfitted<-renderPlot(
        height=reactive({length(strainlist[strainlist$strain==input$query,]$run)*3000}),
        units="px",{
          par(mfrow=c(length(strainlist[strainlist$strain==input$query,]$run)*8,1))
          plotlots(directory, input$query, "strainlist.csv", lag1 = input$initiallag, mumax1 = input$initialmumax, od01 = input$initialod0, odmax1 = input$initialodmax)
        }
      )
      output$Plotstrain <- renderPlot(
        {z <- do.call(rbind, lapply(list.files()[grepl("[pP]late", list.files())], readonedir))
        z$plate <- as.numeric(sub("[pP]late ", "", z$plate))
        z$row <- as.numeric(as.character(z$row))
        z <- left_join(z, strainlist, by = c(plate = "run", row = "column"))
        toplot <- z[z$strain == input$straintoplot, ]
        print(plot(mtcars))
        s <- ggplot(toplot, aes(x = time, y = 1-od, col = factor(plate))) + geom_point() + facet_grid(~temperature) + geom_line(aes(group = well))
        print(s)
        }, height = 1000, width = 1500
        )
      observeEvent(input$do, {
        write.csv(namer(looper(getwd(), lag1 = input$initiallag, mumax1 = input$initialmumax, od01 = input$initialod0, odmax1 = input$initialodmax),"strainlist.csv"),file="outputfits.csv",row.names=F)
        data<-read.csv("outputfits.csv")
      })

        }
  )
}
