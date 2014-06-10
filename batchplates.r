z=c()
startingdir=choose.dir() #change for your directory
setwd(startingdir)				
dirlist=list.files()
for(zz in 1:(length(dirlist)-1)){
	setwd(dirlist[zz])
	filelist <- list.files()												#read file list
	x3=data.frame()															#set an empty dataframe for results
	for (y in 1:length(filelist)){											#now the loop. it loops over all the files in the dir
		t=read.csv(filelist[y]) 											#reads the csvs
		t=t[,2:13] 															#removes the labels and od600
		x1=c() 																#clears a vector for temp results
		for(ii in 1:12){														#this loop takes it from the plate layout to flat
			x1=c(x1,t[,ii])													#adds each column seperately to make a long vector
		}
		x2=c()																#new vector to collate with well names and time
		x2$well=paste(rep (letters[1:8], 12),rep (1:12, each = 8),sep="")	#makes a well column with 96 well names
		x2$od=x1															#moves in the od variables
		x2$time=rep(strsplit(filelist[y],".c",fixed=T)[[1]][1],96)			#this takes the filename and makes it the time column
		x2=as.data.frame(x2)												#puts into a dataframe to use next
		x3=rbind(x3,x2)														#collates in global list
	}
	x3$time=as.numeric(as.character(x3$time))

	library(nlstools)														#loads the need library
	buchworm= od ~ od0 + (time >= lag) * (time <= (lag + (odmax - od0) * 	#this is my modifiction for my factors
		log(10)/mumax)) * mumax * (time - lag)/log(10) + (time >= lag) * 	#not sure about the log(10)s
		(time > (lag + (odmax - od0) * log(10)/mumax)) * (odmax - od0)		#look into this before publishing
	for(ij in 1:length(levels(x3$well))){									#loops over all levels of wells
		workings=x3[x3$well==levels(x3$well)[ij],]							#takes out each well
		workings$od<-1-workings$od											#reverses the OD curve
		workings=workings[workings$time!=0,]
		z=rbind(z,c(coef(nls(buchworm, workings,list(lag=40, mumax=0.02, od0 = 0.2, odmax = 0.9),control=nls.control(minFactor=1/4096,warnOnly=T))),sum(resid(nls(buchworm, workings,list(lag=40, mumax=0.02, od0 = 0.2, odmax = 0.9),control=nls.control(minFactor=1/4096,warnOnly=T)))^2),levels(x3$well)[ij],zz))
	}
	setwd(startingdir)
}
z2=as.data.frame(z)
names(z2)=c("lag","mumax","od0","odmax","residual","well","run")
z2$column=rep(c(1,10,11,12,2,3,4,5,6,7,8,9),times=length(z2$well)/12)
strainlist=read.csv("strainlist.csv")
z2=merge(z2, strainlist, by=c("run","column"))
setwd(startingdir)
write.table(z2,file="grouped.csv",sep=",",row.names=F)