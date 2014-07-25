#checking and installing libraries
libraries=c("nlstools","plyr","dplyr","shiny")
for(i in libraries){if(i %in% rownames(installed.packages()) == FALSE) {install.packages(i)}}
lapply(libraries, require, character.only=T)
#setting working directory
startingdir=choose.dir()
setwd(startingdir)
strainlist=read.csv("strainlist.csv")
#defining some functions!
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

buchworm= od ~ od0 + (time >= lag) * (time <= (lag + (odmax - od0) * 	
    log(10)/mumax)) * mumax * (time - lag)/log(10) + (time >= lag) * 	
    (time > (lag + (odmax - od0) * log(10)/mumax)) * (odmax - od0)
	
compileall<-function(x,y=0.1){
  x=x[which(x$residual<y),]
  z=ddply(x,.(strain,temperature),function(df){c(mean(df$lag),mean(df$mumax),mean(df$od0),mean(df$odmax),
                                             sd(df$lag),sd(df$mumax),sd(df$od0),sd(df$odmax),length(df$lag))})
  names(z)=c("strain","temperature","lag","mumax","od0","odmax","sdlag","sdmumax","sdod0","sdodmax","N")
  return(z)
}

readonedir<-function(x){
  startingdir=getwd()
  setwd(x)
  filelist<-list.files()  											
	x3<-data.frame()															
	for(y in 1:length(filelist)){		
    time=strsplit(filelist[y],".c",fixed=T)[[1]][1]
		x1=data.frame(stringsAsFactors = FALSE)
    t=read.csv(filelist[y],header=TRUE,check.names=FALSE) 											
		t=t[1:13]
    for(zz in 2:length(t[1,])){
      for(xx in 1:length(t[,1])){
        x1=rbind(x1,data.frame(paste(t[xx,1],names(t)[zz],sep=""),t[xx,zz],time,x,stringsAsFactors = FALSE))
        }
      }
    x3=rbind(x3,as.data.frame(x1))
    }
  setwd(startingdir)
  names(x3)=c("well","od","time","plate")
  x3$time=as.numeric(x3$time)
  return(x3)
}

fitbuch<-function(x,y=1){
  x$od<-1-x$od
  if(y==1){x=x[x$time!=0,]}
  z=nls(buchworm,x,list(lag=35, mumax=0.025, od0 = 0.25, odmax = 0.95),control=nls.control(minFactor=1/4096,warnOnly=T))
  z3=data.frame()
  z3=rbind(z3,c(as.data.frame(t(coef(z))),sum(resid(z)^2),x$well[1],x$plate[1],as.numeric(strsplit(as.character(x$plate), split=' ', fixed=TRUE)[[1]][2])))
  names(z3)=c("lag","mumax","od0","odmax","residual","well","plate","run")
  return(z3)
}

plyrfit<-function(x){
  ddply(x,.(plate,well),function(df){fitbuch(df)})
}

looper<-function(x){
  startingdir=getwd()
  z3=data.frame()
  setwd(x)  			
  dirlist=list.files()
  for(zz in 1:(length(dirlist)-1)){
    z3=rbind(z3,plyrfit(readonedir(dirlist[zz])))
    }
  setwd(startingdir)
  return(z3)
}

namer<-function(x,y){
  strainlist=read.csv(y)
  x$column=as.numeric(substring(as.character(x$well), 2, nchar(as.character(x$well))))
  x$row=substring(as.character(x$well), 1,1)
  z2=merge(x, strainlist, by=c("run","column"))
  asNumeric <- function(x){as.numeric(as.character(x))}
  return(z2)
}

plotlots<-function(x,y,z){ #function of directory, strain, strainlist
  startdir=getwd()
  setwd(x)
	strainlist=read.csv(z)
	strainlist=strainlist[which(strainlist$strain==y),]
  if(length(strainlist[,1])==0){return("strain not found")}
	for(i in 1:length(strainlist[,1])){
		singlehold=readonedir(paste("plate",sprintf("%02d", strainlist[i,]$run)))
		singlehold=singlehold[which(substring(singlehold$well, 2, nchar(x))==strainlist[i,]$column),]
		singlehold$temperature=strainlist[i,]$temperature
		d_ply(singlehold,.(well),plotter)
	}
	setwd(startdir)
}

plotter<-function(x,well,strain){
	x=x[x$time!=0,]
	x$od=1-x$od
	workingfit<-nls(buchworm, x,list(lag=35, mumax=0.025, od0 = 0.25, odmax = 0.95),control=nls.control(minFactor=1/4096,warnOnly=T))
	plot(x$od~x$time,ylab="od",xlab="hours",main=c(paste("Plate:", x$plate[1],"Temp:", x$temperature[1],"Well",x$well[1],"Residual:",round(sum(resid(workingfit)^2),4))))
	points(x$time,predict(workingfit,list(time=x$time)),pch=3,col=2)
}

#delete files from previous runs - deletes, otherwise makes it
if(file.exists("grouped.csv")){file.remove("grouped.csv")}
if(!file.exists("outputfits.csv")){write.csv(namer(looper(getwd()),"strainlist.csv"),file="outputfits.csv",row.names=F)}

#read the csv in
data<-read.csv("outputfits.csv")
write.csv(compileall(data,0.1),file="grouped.csv",row.names=F)
data$factortemp=as.factor(data$temperature)
#shiny!
shinyServer(function(input, output) {
	output$Plot1 <- renderPlot({
		par(mfrow=c(length(levels(data$factortemp)),1))
		for(i in 1:length(levels(data$factortemp))){
			data2<-data[which(data$factortemp==levels(data$factortemp)[i]),]
			dataclean <- data2[data2$residual<input$cutoff,]
			plot(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,xlab="strain",ylab=as.character(input$variable),main=c(levels(data$factortemp)[i]," degrees"))}
	}, height = 1000, width = 1500)

	output$table1 <- renderTable({
		compileall(data,input$cutoff)
	},digits=5)
	
	output$Plot3 <- renderPlot({	
		par(mfrow=c(length(levels(data$factortemp)),1))
		for(i in 1:length(levels(data$factortemp))){
			data2<-data[which(data$factortemp==levels(data$factortemp)[i]),]
			dataclean <- data2[data2$residual<input$cutoff,]
			plotmeans(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,n.label=F,barcol="black",connect=F,main=c(levels(data$factortemp)[i]," degrees"),xlab="Strain",ylab=as.character(input$variable),p=0.68)
			}
	}, height = 1000, width = 1500)
	
	output$Plotfitted <- renderPlot(height = reactive({length(strainlist[strainlist$strain==input$query,]$run)*3000}), units="px",{
		par(mfrow=c(length(strainlist[strainlist$strain==input$query,]$run)*8,1))
		plotlots(startingdir,input$query,"strainlist.csv")
	})	
})
