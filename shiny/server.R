#checking and installing libraries
library(platereading)

startingdir=getwd()
strainlist=read.csv("strainlist.csv")
#defining some functions!
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

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
			if(input$ordered1==T){dataclean$strain<-factor(dataclean$strain,names(sort(rank(tapply(subset(dataclean,select=input$variable)[,1],dataclean$strain,median),ties.method="first"))))}
			plot(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,xlab="",ylab=as.character(input$variable),main=c(levels(data$factortemp)[i]," degrees"),las=2)}
	}, height = 1000, width = 1500)

	output$table1 <- renderTable({
		compileall(data,input$cutoff)
	},digits=5)

	output$Plot3 <- renderPlot({
		par(mfrow=c(length(levels(data$factortemp)),1))
		for(i in 1:length(levels(data$factortemp))){
			data2<-data[which(data$factortemp==levels(data$factortemp)[i]),]
			dataclean <- data2[data2$residual<input$cutoff,]
			if(input$ordered2==T){dataclean$strain<-factor(dataclean$strain,names(sort(rank(tapply(subset(dataclean,select=input$variable)[,1],dataclean$strain,mean),ties.method="first"))))}
			plotmeans(dataclean[,which(names(dataclean)==input$variable)]~dataclean$strain,n.label=F,barcol="black",connect=F,main=c(levels(data$factortemp)[i]," degrees"),xlab="",ylab=as.character(input$variable),p=0.68,par(las =2))
			}
	}, height = 1000, width = 1500)

	output$Plotfitted <- renderPlot(height = reactive({length(strainlist[strainlist$strain==input$query,]$run)*3000}), units="px",{
		par(mfrow=c(length(strainlist[strainlist$strain==input$query,]$run)*8,1))
		plotlots(startingdir,input$query,"strainlist.csv")
	})
})
