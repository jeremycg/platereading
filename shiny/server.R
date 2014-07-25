library(shiny)
library(gplots)
library(nlstools)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
startingdir=choose.dir() #change for your directory
setwd(startingdir)
data<-read.csv("outputfits.csv")
strainlist=read.csv("strainlist.csv")
x3=data.frame()
buchworm= od ~ od0 + (time >= lag) * (time <= (lag + (odmax - od0) * 	#this is my modification for my factors
    log(10)/mumax)) * mumax * (time - lag)/log(10) + (time >= lag) * 	#not sure about the log(10)s
    (time > (lag + (odmax - od0) * log(10)/mumax)) * (odmax - od0)



shinyServer(function(input, output) {
  output$Plot1 <- renderPlot({
	data2<-data[data$temperature==20,]
    dataclean <- data2[data2$residual<input$cutoff,]
    plot(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,xlab="strain",ylab=as.character(input$variable),main="20 degrees")
  })
  output$Plot2 <- renderPlot({
	data2<-data[data$temperature==29,]
    dataclean <- data2[data2$residual<input$cutoff,]
    plot(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,xlab="strain",ylab=as.character(input$variable),main="29 degrees")
  })
  output$table1 <- renderTable({
		data2<-data[data$temperature==20,]
		dataclean <- data2[data2$residual<input$cutoff,]
		table11=c()
		table11$strain=names(as.data.frame(t(tapply(dataclean[,which(names(dataclean)==input$variable)],dataclean$strain,mean))))
		table11$mean20=unname(tapply(dataclean[,which(names(dataclean)==input$variable)],dataclean$strain,mean))
		table12=c()
		table12$strain=names(as.data.frame(t(tapply(dataclean[,which(names(dataclean)==input$variable)],dataclean$strain,stderr))))
		table12$stderr20=unname(tapply(dataclean[,which(names(dataclean)==input$variable)],dataclean$strain,stderr))
		table3=merge(table11,table12,by="strain")
		data2<-data[data$temperature==29,]
		dataclean2 <- data2[data2$residual<input$cutoff,]
		table13=c()
		table13$strain=names(as.data.frame(t(tapply(dataclean2[,which(names(dataclean2)==input$variable)],dataclean2$strain,mean))))
		table13$mean29=unname(tapply(dataclean2[,which(names(dataclean2)==input$variable)],dataclean2$strain,mean))
		table14=c()
		table14$strain=names(as.data.frame(t(tapply(dataclean2[,which(names(dataclean2)==input$variable)],dataclean2$strain,stderr))))
		table14$stderr29=unname(tapply(dataclean2[,which(names(dataclean2)==input$variable)],dataclean2$strain,stderr))
		table4=merge(table13,table14,by="strain")
		merge(table3,table4,by="strain")
		
	},digits=5)
	output$Plot3 <- renderPlot({	
	data2<-data[data$temperature==20,]
    dataclean <- data2[data2$residual<input$cutoff,]
	plotmeans(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,n.label=F,barcol="black",connect=F,main="20 degrees",xlab="Strain",ylab=as.character(input$variable),p=0.68)
	})
	output$Plot4 <- renderPlot({	
	data2<-data[data$temperature==29,]
    dataclean <- data2[data2$residual<input$cutoff,]
	plotmeans(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,n.label=F,barcol="black",connect=F,main="29 degrees",xlab="Strain",ylab=as.character(input$variable),p=0.68)
	})
	output$Plotfitted <- renderPlot(height = reactive({length(strainlist[strainlist$strain==input$query,]$run)*3000}), units="px",{

xx=strainlist[strainlist$strain==input$query,]
xx$run=sprintf("%02d", xx$run)
x3=data.frame()
for(i in 1:length(xx$run)){
	setwd(startingdir)
	setwd(paste("plate", xx$run[i]))
	#do the stuff here!!!
	filelist <- list.files()												#read file list
				
	for (y in 1:length(filelist)){											#now the loop. it loops over all the files in the dir
		t=read.csv(filelist[y]) 											#reads the csvs
		t=t[,2:13] 													#removes the labels and od600
		t=t[,xx$column[i]]											#takes the column we want
		
		x2=c()																#new vector to collate with well names and time
		x2$letter=c("a","b","c","d","e","f","g","h")	
		x2$od=t															#moves in the od variables
		x2$time=rep(strsplit(filelist[y],".c",fixed=T)[[1]][1],8)		#this takes the filename and makes it the time column
		x2$plate=rep(xx$run[i],8)
		x2$column=rep(xx$column[i],8)
		x2$temperature=rep(xx$temperature[i],8)
		x2=as.data.frame(x2)												#puts into a dataframe to use next
		x3=rbind(x3,x2)														#collates in global list
	}
	setwd(startingdir)
}
for(i in 1:length(x3$plate)){
	x3$group[i]=paste(x3$letter[i],x3$plate[i],x3$column[i],x3$temperature[i],sep=" ")
	}
x3$group=as.factor(x3$group)
par(mfrow=c(length(levels(x3$group)),1))
for(i in 1:length(levels(x3$group))){
	workings<-x3[x3$group==levels(x3$group)[i],]
	workings$od<-1-workings$od
	workings$time<-as.numeric(as.character(workings$time))
	workings<-workings[workings$time!=0,]
	workingfit=nls(buchworm, workings,list(lag=40, mumax=0.02, od0 = 0.2, odmax = 0.9),control=nls.control(minFactor=1/4096,warnOnly=T))
	plot(workings$od~workings$time,ylab="od",xlab="hours",main=c(paste("Temp:", workings$temperature[1],"Plate:", workings$plate[1],"Row:",workings$letter[1],"Residual:",round(sum(resid(workingfit)^2),4))))
	points(workings$time,predict(workingfit,list(time=workings$time)),pch=3,col=2)
	}
	})	
})
