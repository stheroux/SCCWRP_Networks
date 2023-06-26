setwd("/Users/Mel/Desktop/SCCWRP")
library(leaps)
library(tidyverse)
library(caret)
library(MASS)

## This script takes the output of SpiecEasi (network topology table, produced from provided SCCWRP data) and determines which network topology variables are the best predictors of reference status.
## This analysis employs ordered logistic regression models with the 'polr' command in the 'MASS' package.
## Below, we use two methods (stepwise (AIC) and leaps) to determine the best model for each taxonomic level and primer. 

###################
###
### Load and process data
###
###################

#Load data
topology = read.csv("network_topology.csv", row.names = 1)
topology <- subset(topology, Num.Subsampled != "All Samples")
# remove 'inf' values
topology = topology[!grepl("Inf", topology$Pos.Neg),]
# convert reference status to ordered factor
topology$Reference.Status <- ordered(topology$Reference.Status, 
                                     levels = c("physeq_Reference", "physeq_Intermediate", "physeq_Stressed"))

###################
###
###  Record topology summary and prepare data for model selection
###
###################

primers = c("asv16S","asv18S","asvrbcL")
levels = c("Family","Genus","Species")
topsum = data.frame()
report = data.frame()

for(primer in primers){
  
  topology_primer <- subset(topology, Primer == primer)
  
  for(level in levels){
    
    topology_sub <- subset(topology_primer, Taxonomic.Level == level)
    
    # Record summary of topology variables
    vars = as.vector(colnames(topology)[-c(1:4)])
    for(var in vars){
      min = summary(topology[[var]])[1]
      max = summary(topology[[var]])[6]
      mean = summary(topology[[var]])[3]
      sd = sd(topology[[var]])
      
      temp <- c(primer, level, var, min, max, mean, sd)
      topsum = rbind(topsum, temp)
    }
    colnames(topsum) <- c("Primer","Taxonomic Level", "Network Topology Variable", "Minimum","Maximum","Mean","SD")
    
    # Prepare data for model selection
    # Remove multicollinear variables, size-dependent variables, constants, primer, taxonomic level
    subtop = subset(topology_sub, select=-c(Primer, Taxonomic.Level, Num.Subsampled, Positive.Edges,Negative.Edges,Pos.Total,Neg.Total, Avg.Degree, Total.Edges, Total.Nodes))
    
    
    ###################
    ###
    ### AIC/steps method of determining best fit
    ### http://www.sthda.com/english/articles/36-classification-methods-essentials/150-stepwise-logistic-regression-essentials-in-r/
    ###
    ###################
    
    # Split the data into training and test set
    training.samples <- subtop$Reference.Status %>% 
      caret::createDataPartition(p = 0.8, list = FALSE)
    train.data  <- subtop[training.samples, ]
    test.data <- subtop[-training.samples, ]
    
    ## Remove two outliers from Family asv18S data
    if(level == "Family" & primer == "asv18S"){
      train.data <- subset(train.data, Pos.Neg < 100)
    }
    
    ## Full logistic regression model. (Incorporating all predictor variables)
    full.model <- MASS::polr(Reference.Status ~., data = train.data)
    #coef(full.model)
    predictors_full = as.character(names(coef(full.model)))
    predictors_full = c(predictors_full, "Reference.Status")
    predictors_full
    subtop_full = subtop[predictors_full]
    
    ## Perform stepwise variable selection. (Select the most contributive variables)
    step.model <- full.model %>% MASS::stepAIC(trace = FALSE)
    #coef(step.model)
    predictors_step = as.character(names(coef(step.model)))
    predictors_step = c(predictors_step, "Reference.Status")
    predictors_step
    subtop_step = subtop[predictors_step]
    
    ###################
    ###
    ### leaps method of determining best fit
    ### https://towardsdatascience.com/selecting-the-best-predictors-for-linear-regression-in-r-f385bf3d93e9
    ###
    ###################
    
    # Run the regsubsets() function on all variables.
    Best_Subset <-
      regsubsets(Reference.Status~.,
                 data =subtop,
                 nbest = 1,      # 1 best model for each number of predictors
                 nvmax = NULL,    # NULL for no limit on number of variables
                 force.in = NULL, force.out = NULL,
                 method = "exhaustive")
    summary_best_subset <- summary(Best_Subset)
    as.data.frame(summary_best_subset$outmat)
    which.max(summary_best_subset$adjr2) #see what the package recommends in terms of the number of predictors to use for our dataset
    predictors = summary_best_subset$which[which.max(summary_best_subset$adjr2),] #What are the best predictors? best predictors are indicated by ‘TRUE’.
    predictors = names(predictors)[predictors]
    predictors = predictors[-1]
    predictors_leaps = c(predictors, "Reference.Status")
    rm(predictors)
    predictors_leaps #print
    subtop_leaps = subtop[predictors_leaps]
    
    ## leaps model
    leaps.model <- MASS::polr(Reference.Status ~., data = subtop_leaps)

    ###################
    ###
    ### Compare full, stepwise, and leaps models
    ###
    ###################
    
    # The best model is defined as the model uses the fewest predictor variables without compromising accuracy. 
    
    # Prediction accuracy of the full logistic regression model:
    # Make predictions
    probabilities <- full.model %>% predict(test.data, type = "probs")
    predicted.classes <- colnames(probabilities)[apply(probabilities,1,which.max)]
    # Model accuracy (Full)
    observed.classes <- test.data$Reference.Status
    full.accuracy <- mean(predicted.classes == observed.classes)
    
    # Prediction accuracy of the stepwise logistic regression model:
    # Make predictions
    probabilities <- step.model %>% predict(test.data, type = "probs")
    predicted.classes <- colnames(probabilities)[apply(probabilities,1,which.max)]
    # Model accuracy (Stepwise)
    observed.classes <- test.data$Reference.Status
    step.accuracy <- mean(predicted.classes == observed.classes)
    
    # Prediction accuracy of the leaps logistic regression model:
    # Make predictions
    probabilities <- leaps.model %>% predict(test.data, type = "probs")
    predicted.classes <- colnames(probabilities)[apply(probabilities,1,which.max)]
    # Model accuracy (Leaps)
    observed.classes <- test.data$Reference.Status
    leaps.accuracy <- mean(predicted.classes == observed.classes)
    
    ###################
    ###
    ### Select which model to use.
    ###
    ###################   
    
    # Which model has the best accuracy?
    accuracy_list <- c("full.model"=full.accuracy, "step.model"=step.accuracy, "leaps.model"=leaps.accuracy)
    max_models <- names(which(accuracy_list==max(accuracy_list))) # returns indices with max value (including ties)
    
    # If there's a tie, which model uses fewest predictors?
    if(length(max_models) > 1){
      
      mod1 <- get(paste("predictors_",sub(".model","",max_models[1]), sep = ""))
      mod2 <- get(paste("predictors_",sub(".model","",max_models[2]), sep = ""))
      
      if(length(mod1)<length(mod2)){
        best_model <- get(max_models[1])
        print(paste(max_models[1], "has fewer variables than",max_models[2]))
        print(paste("best model is",max_models[1]))
      } else if(length(mod2)<length(mod1)){
        best_model <- get(max_models[2])
        print(paste(max_models[2], "has fewer variables than",max_models[1]))
        print(paste("best model is",max_models[2]))
      } else {
        best_model <- get(max_models[1])
      }
      
    } else if (length(max_models) == 1){
      best_model <- get(max_models[1])
      print(paste("best model is",max_models[1]))
    } 

    
    ###################
    ###
    ### Summarize and record odds ratio for best model. 
    ###
    ###################    
    
    # get summary
    summary(best_model)
    # get coefficients (it's in matrix form)
    coefficients <- summary(best_model)$coefficients
    # calculate p-values
    p_value <- (1 - pnorm(abs(coefficients[ ,"t value"]), 0, 1))*2
    # bind back to coefficients
    coefficients <- cbind(coefficients, p_value)
    
    # calculate odds ratios
    odds_ratio <- exp(coefficients[ ,"Value"])
    # combine with coefficient and p_value
    # (coefficients <- cbind(
    #   coefficients[ ,c("Value", "p_value")],
    #   odds_ratio
    # ))
    
    table <- (coefficients <- cbind(
      coefficients[ ,c("Value", "p_value")],
      odds_ratio
    ))
    
    table <- table[1:length(best_model$coefficients),]
    table <- cbind(topology = rownames(table), table)
    table <- cbind(TaxLevel = level, table)
    table <- cbind(ASVprimer = primer, table)
    row.names(table) <- NULL
    
    report <- rbind(report, table)
    
    # #Calculate percent lower/higher odds 
    # (percent_odds <- ifelse(p_value < 0.05, paste0(as.character(round((coefficients[,"odds_ratio"]-1)*100, 2)), "%"), "N/A (p > 0.05)"))
    
  } #for level in levels
} #for primer in primers



#Remove rows with p value > 0.05
report$p_value <- as.numeric(report$p_value)
report <- subset(report, p_value < 0.05)

#Write statements that interpret odds ratio
#NOTE: Odds ratio values can be challenging to interpret. 
#For ease of interpretation, the reciprocal for values with negative exponents is calculated 
#to represent how “less likely” the odds of increased site stress is with each one unit increase 
#in the corresponding network topology factor and is recorded in the "Statement" column. 
#Values with positive exponents are interpreted as that much “more likely” to have increased 
#stress with each one unit increase in the corresponding network topology factor.
report$odds_ratio <- as.numeric(report$odds_ratio)
report <- report %>% 
  mutate(Interpretation = if_else(report$odds_ratio>=1, "leave", "reciprocal"))
report <- report %>% 
  mutate(Statement = if_else(report$Interpretation=="leave", 
                             paste("For each one unit increase in",report$topology,"increased stress is",report$odds_ratio,"times more likely"), 
                             paste("For each one unit increase in",report$topology,"increased stress is",1/report$odds_ratio,"times less likely")))
report <- subset(report, select = -Interpretation)

write.csv(topsum, "/Users/Mel/Desktop/topology_summary.csv")
write.csv(report, "/Users/Mel/Desktop/prediction_summary.csv")



