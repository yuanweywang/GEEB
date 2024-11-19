#Generalized Estimating Equations Boosting (GEEB) Machine for Correlated Data

#Authors: Yuan-Wey Wang, Hsin-Chou Yang, Yi-Hau Chen, Chao-Yu Guo (Corresponding author)

#Corresponding author:Chao-Yu Guo, Ph. D., Professor, Division of Biostatistics and Data Science
#Institute of Public Health, College of Medicine, National Yang Ming Chiao Tung University, Taipei, Taiwan, ROC

#Email address:
#Yuan-Wey Wang: vivian4012@gmail.com
#Hsin-Chou Yang: hsinchou@stat.sinica.edu.tw
#Yi-Hau Chen: yhchen@stat.sinica.edu.tw
#Chao-Yu Guo: cyguo@nycu.edu.tw


##source of cround: http://andrewlandgraf.com/2012/06/15/rounding-in-r/
cround <- function(x, digits=0) {
  posneg = sign(x)
  z = abs(x)*10^digits
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^digits
  z*posneg
}


geebm <- function(formula, id, iteration=100, feature_rate=0.5, lrate=0.1, standardize=F, data){
  library(geepack)
  Features <- as.vector(all.vars(formula))[-1]  
  outcome  <- as.character(formula[[2]])        
  data <- subset(data,select=(c(Features,id,outcome)))  
  data_d <- as.data.frame(model.matrix(~., data=data[Features]))  
  Features_d <- names(data_d)                                     
  if(standardize==T){
    data_d[, !(names(data_d) %in% c(id,"(Intercept)"))]<-scale(data_d[, !(names(data_d) %in% c(id,"(Intercept)"))])
  }
  data_d <- cbind(data_d,data[id],data[outcome])               
  Coefficient <- setNames(data.frame(matrix(0,ncol = length(Features_d), nrow = iteration)),Features_d)  
  PredValue <- mean(data_d[[outcome]]) 
  data_d <- cbind(data_d, PredValue)
  for (m in (1:iteration)){
    R <- as.numeric(sample(1:length(Features),size=(cround(feature_rate*length(Features))),replace=F))
    R_name_d <- grep(paste(Features[R], collapse = "|"), colnames(data_d), value = TRUE)
    data_rf_d <- data_d[,R_name_d]
    data_rf_d["Residual"] <- data_d[[outcome]]-data_d[["PredValue"]]
    data_rf_d <- cbind(data_rf_d,data_d[id[1]])
    modelGEEBM <- eval(parse(text = paste0("geeglm(Residual ~ . - ", id[1], ", id =", id[1],", data = data_rf_d)")))
    for (i in seq_along(Coefficient)) {
      if (colnames(Coefficient)[i] %in% names(coef(modelGEEBM))) {
        Coefficient[m, i] <- coef(modelGEEBM)[colnames(Coefficient)[i]] 
      }
    }
    PredResidual <- predict(modelGEEBM, newdata = data_d)
    data_d["PredResidual"] <- PredResidual
    data_d['PredValue'] <- data_d$PredValue+lrate*data_d$PredResidual
  }
  model <- structure(list(lrate=lrate, PredValue=PredValue, Features=Features, Features_d=Features_d, Coefficient=Coefficient, standardize = standardize), class = "geebm" )
  return(model)
}


predict.geebm <- function(model, data){
  data_d <- model.matrix(~., data = data[model$Features])
  if(model$standardize==T){
    data_d[,!(apply(data_d, 2, function(col) all(col == col[1])))]<-scale(data_d[,!(apply(data_d, 2, function(col) all(col == col[1])))])
  }
  prediction  <- model$PredValue+model$lrate*rowSums(data_d %*% t(model$Coefficient))
  return(as.vector(prediction))
}


geebm.importance <- function(model){
  importance <- setNames(data.frame(matrix(0,ncol = 2, nrow = length(model$Features_d))),c("Feature","Importance"))
  importance[,2] <- apply(model$Coefficient, 2, function(col) mean(col[col != 0]))
  importance[,1] <- model$Features_d
  return(importance)
}


geebm.importance.plot <- function(importance,filename="FeatureImportance.jpg"){
  library(ggplot2)
  plot <- ggplot(importance, aes(x = Importance, y = Feature)) +
    geom_bar(stat = "identity", fill = "#003060") +
    labs(x = NULL, y = NULL) +
    coord_flip()+
    theme(axis.line = element_line(color = "black"), 
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(), 
          axis.text.y = element_text(color = "black", size = 13),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 14),
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.y = element_line(color = "black", size = 0.5))+ 
    scale_x_continuous(breaks = seq(-0.03, max(importance[["Importance"]]), by = 0.005))
  print(plot)
  ggsave(filename,height=3,width=5,dpi=300,units="in",scale=2, plot=plot)
}









