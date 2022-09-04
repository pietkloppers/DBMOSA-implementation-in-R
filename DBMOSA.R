
#Histogram diversity implementation
hypergrid<-function(dataframe){
  #Assigning each entry of the archive to a neighborhood
  for(i in 1:nrow(dataframe)){
    if(dataframe[i,1] < 2 & dataframe[i,2] < 2){
      dataframe[i,3] <- 1
    }
    else if(dataframe[i,1] < 2 & dataframe[i,2] >= 2){
      dataframe[i,3] <- 3
    }
    else if(dataframe[i,1] >= 2 & dataframe[i,2] < 2){
      dataframe[i,3] <- 2
    }
    else if(dataframe[i,1] >= 2 & dataframe[i,2] >= 2){
      dataframe[i,3] <- 4
    } 
  }
  #Just ensuring the subset function uses the correct column
  colnames(dataframe)<-c("f1","f2","Grid_number")
  
  #Only start calculating the size of the neighborhood as a percentage of the 11th entry of the archive as the 
  #first few entries will have very high percentages
  if(nrow(dataframe) > 10 ){
    #calculating the percentage
    percentage<-nrow(subset(dataframe, Grid_number == dataframe[nrow(dataframe),3]))/nrow(dataframe)*100 
    #rejecting if the percentage is too high
    if(percentage > 60){
      returnv = 0
      #accepting if the percentage is lower than the stipulated threshold  
    }else{
      returnv = 1
    }
    #accepting the first 10 entries  
  }else{
    returnv = 1
  }
  
  return(returnv)
}

#crowding diversity measure

crowding<-function(tada){
  
  
  #Inserting condition to start after 3 entries in the Archive to ensure that cubicles can be formed
  if(nrow(tada) > 3){
    #Saving the latest entry of the archive to a different variable to extract information later
    finalrow<-tada[nrow(tada),]
    
    #ordering the dataframe
    tada<-tada[order(tada[,1]),]
    
    #assigning "infinite" value to the first circumference
    tada$circumference[1]<- 1000000
    
    #Calculating the circumference of the cubicles
    for(i in 1:(nrow(tada)-2)){
      tada[i+1,3] <- (abs(tada[i,1]-tada[i+2,1])*2)+(abs(tada[i,2]-tada[i+1,2])*2)
    }
    
    #assigning "infinite" value to the last circumference
    tada$circumference[nrow(tada)]<-100000
    
    #Setting default return value
    returnv = 1
    
    #Stipulating that if the circumference is smaller than 0.15 the entry will be rejected
    for(i in 1:nrow(tada)){
      if(finalrow[1,1] == tada[i,1]  &  finalrow[1,2] == tada[i,2] & tada[i,3] < 0.25){
        returnv = 0
        break
      }
    } 
  }else{
    returnv = 1
  }
  
  return(returnv)
}

#Dominance function
dominates<-function(df, Temperature,t,d,c){
  countR<-data.frame(x=0)
  countA<-data.frame(x=0)
  counter = 1
  counterA = 1
  #getting the count of times the new solution dominates or is dominated
  for(i in 1:(nrow(df)-1)){
    if(df[nrow(df),1] >= df[i,1]  & df[nrow(df),2] >= df[i,2]){
      countR[counter,1]<-counter
      counter= counter+1
    }
    if(df[nrow(df),1] <= df[i,1]  & df[nrow(df),2] <= df[i,2]){
      countA[counter,1]<-counter
      counterA= counterA+1
    }
  }
  #Calculating the Energy 
  Energy <-(max(countR)- max(countA))/nrow(unique(df))
  #Calculating the probability
  prob<-c(1,0)
  prob[2]<- exp(-Energy/Temperature)
  
  probability<-min(prob)
  
  #Rejecting the solution if the random number is bigger
  if(runif(1) > probability){
    df <- df[-nrow(df),]
    t<-t+1
    d<-d+1
  }
  #Diversity-based criterion (histogram)
  #else if(hypergrid(df) == 0){
  #  df <- df[-nrow(df),]
  #  t<-t+1
  #  d<-d+1
  #}
  #Diversity-based criterion (Crowding)
  else if(crowding(df) == 0){
    df <- df[-nrow(df),]
    t<-t+1
    d<-d+1
  }
  else{
    for(i in 1:(nrow(df)-1)){
      if(df[nrow(df),1] <= df[i,1]  &  df[nrow(df),2] <= df[i,2]){
        df[i,1] <- -1
        df[i,2] <- -1
      }
    }
    #Adding the solution
    c=c+1
    df <- subset(df, f1 > -0.95 & f2 > -0.95 )
    t=t+1
  }
  
  return(list(df,t,d,c))
}


#generating new x neighbours
generate_neigh<-function(my_ex){
  #min value of either 0.4 , 10^5 minus x or x is accepted- both the objective function values 
  #are quadratic and thus negative values are not required
  x1<-my_ex[nrow(my_ex),1] + (0.5-runif(1))*min(0.4,10^5-my_ex[nrow(my_ex),1],my_ex[nrow(my_ex),1])
  
  dief<-rbind(my_ex,x1)
  
  return(dief)
}

#obtaining the objective space values
gen_f<-function(xe){
  output<-data.frame(f1=0,f2=0)
  #f1
  output[1,1]<-(xe[nrow(xe),1])^2
  #f2
  output[1,2]<-((xe[nrow(xe),1])-2)^2
  
  return(output)
}


#Initializing parameters
initial=data.frame(x1=1)
Archive = data.frame(f1=1,f2=1)

#from preliminary experiments the standard deviation of the objective function is around 0.25. 
#For reproducibility is was used 

#Acceptance deviation for initial Temp
Temperature <- (-3*0.25)/log(0.8, base = exp(1))

#Accept all
#Temperature <- 5

epochs <- 1
accept<- 0
reject<- 0
iterations <- 1

epochs_max <- 10
rejects_max <- 20
accept_max <- 20

#Starting DBMOSA

#Reaching final temperature
#while(Temperature < 10 & Temperature > 0.1){

#Predefined number of epochs
while(epochs < epochs_max){
  #Adaptive paradigm
  if(reject == rejects_max){
    #Reheating schedule variation
    Temperature = Temperature * 1.2
    #Temperature = Temperature + 0.2
    epochs = epochs + 1
    accept = 0
    reject = 0
  }
  if(accept == accept_max){
    #Cooling schedule variation
    Temperature = Temperature * 0.9
    #Temperature = Temperature - 0.1
    epochs = epochs + 1
    accept = 0
    reject = 0
  }
  
  #Static paradigm
  #if(iterations%%5 == 0){
  #  epochs = epochs + 1
  #  accept = 0
  #  reject = 0
  #  if(iterations %% 2 ==0){
  #    Temperature = Temperature*1.1
  #  }
  #  if(iterations %% 2 ==1){
  #    Temperature = Temperature*0.9
  #  }
  #}
  
  #generating more x values
  initial<-generate_neigh(initial)
  #generating the neighboring f values from the updated x values
  fprime<-gen_f(initial)
  #updating the archive before the dominance check
  Archive<-rbind(Archive,fprime)
  #Dominance check
  placeholder<-dominates(Archive, Temperature , iterations, reject ,accept)
  
  #Updating archive after the dominance check
  Archive <-placeholder[[1]]
  #Updating the iterations after dominance check
  iterations <-placeholder[[2]]
  #Updating the rejections after dominance check
  reject <- placeholder[[3]]
  #Updating the acceptances after dominance check
  accept <- placeholder[[4]]

}

#F1 vs F2
plot(Archive,xlim=c(0,4), ylim=c(0,4), main = "Crowding Pareto Front")

#True pareto front
x <- seq(0,2,0.05)
y1<- x^2
y2 <- (x-2)^2

a<-initial[1:172,]

plot(y1,y2, xlim=c(0,4), ylim = c(0,4) , main = "True Pareto Front")

plot(x,y1, xlim = c(0,2), ylim = c(0,4),main = "Pareto Optimal Set" , xlab = "Optimal set values", ylab = "Objective space values" )
points(x, y2 , xlim =c(0,2), col = 2)
legend("topright",c("f1", "f2"), pch =15,col=c(1,2), title = "Objective space",)

#write for results
#write.csv(Archive,"Crowding.csv",row.names=FALSE)

#F1 vs F2
plot(Archive, xlim=c(0,4), ylim=c(0,4), main = "Approximate Pareto Front")




