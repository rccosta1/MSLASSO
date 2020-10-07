
########################
# The code below accompanies Costa, R. and Baker, J. (2020) Factors Influencing Business Recovery After the 2011 Tohoku Earthquake
# Synthethic data is provided to test the code. This is NOT the data used for the paper.
#
# Directions:
# 1. Download the synthetic data from: http://web.stanford.edu/~bakerjw/gm_selection.html
# 2. Install all packages required by the MSLASSOAlgorithm function
# 3. Load the MSLASSOAlgorithm function
# 4. Run the analysis for B random splits, this will automatically generate three plots
        B = 10
        results <- MSLASSOAlgorithm(B) 
        
# 5. Plot the model summary
        summary(final.fit)
        
# Note: the synthetic data generates          
#######################

# Main function
MSLASSOAlgorithm <- function(n){
  library(visdat)
  library(tidyverse)
  library(tidyselect)
  library(DMwR)
  library(glmnet)
  library(MALDIquant)
  library(pROC)
  library(ggplot2)
  
  ##############################################
  # Get and clean the data
  ##############################################
  # Read the main file with the data (should be .csv)
  data = read.csv(file = "SyntheticData.csv",header = TRUE, na.strings = "NA")
  names(data)[ncol(data)] = "y"
  data$y = as.factor(data$y)
  
  # Reset the random number generator seed to have the same results consistently
  set.seed(13)
  
  # Separate the data into train and test sets
  train.rows <- sample(1:nrow(data),0.75*nrow(data))
  train <- data[train.rows,]
  test <- data[-train.rows,]
  
  # Get the odds ratio of the data to use in the SMOTE function
  rate_ones <- sum(data$y==1)
  rate_zeroes <- sum(data$y==0)
  oddsratio <- 0
  
  if(rate_zeroes > rate_ones){
    oddsratio <- rate_ones/rate_zeroes
  } else {
    oddsratio <- rate_zeroes/rate_ones
  }
  
  # Get the SMOTEd data
  smoted <- SMOTE(y ~., train, perc.over = 100/oddsratio, proc.under = 100)
  
  # Extract the X and Y vectors as dummy variables from the data
  train.factors <- as.matrix(model.matrix(~.,smoted)[,-1])
  x.train <- train.factors[,-ncol(train.factors)]
  y.train <- train.factors[,ncol(train.factors)]
  
  # Get variables names
  var.names <- colnames(x.train)
  
  # Variable importance measure data frame
  VIM <- data.frame("Variables" = var.names, "VIM" = 0)
  
  # Temporary matrices
  beta.lasso <- data.frame("Variables" = var.names, "Multiplier" = 0)
  beta.all <- data.frame("Variables" = c("(Intercept)",var.names))
  betas.union <- data.frame("Variables" = c("(Intercept)",var.names), "Mean" = 0, "Stdev" = 0, "Test" = 0)
  betas.convo <- data.frame("Variables" = c("(Intercept)",var.names), "Mean" = 0, "Stdev" = 0, "Test" = 0)
  p.final <- data.frame("Variables" = var.names)
  e.final <- data.frame("Variables" = var.names)
  s.final <- data.frame("Variables" = var.names)
  
  
  ##############################################
  # Multi-split Algorithm
  ##############################################
  for (i in 1:n){
    
    # Dataframe to hold the raw values of p-values, estimates and standard errors
    raw <- data.frame("Variables" = var.names, "p.raw" = 1, "estimate.raw" = 0, "error.raw" = 0)
    
    # 1. randomly split the data into screening and cleaning sets
    screen.rows <- sample(1:nrow(train.factors),0.5*nrow(train.factors))
    screen <- train.factors[screen.rows,]
    x.screen <- screen[,-ncol(screen)]
    y.screen <- screen[,ncol(screen)]
    
    # 2. Use LASSO on the screening set to find the meaning predictors, i.e., Beta LASSO <> 0
    fit <- cv.glmnet(x.screen, y.screen, alpha = 1, type.measure = "deviance", nfold=10, family = "binomial")
    
    # Get the best lambda
    min.lambda = fit$lambda.min
    best.lambda = fit$lambda.1se
    
    # Get the model coefficients
    tmp_coeffs <- coef(fit, s = best.lambda)
    coefs <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    
    # 3. Obtain the least square estimate of Beta for the predictors selected with LASSO for the cleaning set
    for (k in 1:nrow(beta.lasso)){
      for (j in 1:nrow(coefs)){
        if (coefs[j,1] == beta.lasso[k,1]){
          beta.lasso[k,2] <- 1
        }
      }
    }
    
    
    # Cleaning data set
    clean.all <- train.factors[-screen.rows,]
    clean <- clean.all
    
    # Select only the predictors with Beta.lasso > 0 for the next step
    for (j in 1:nrow(beta.lasso)){
      for (k in 1:ncol(clean.all)){
        if (beta.lasso[j,1] == colnames(clean.all)[[k]] & beta.lasso[j,2] == 0.0){
          index <- (grep(paste0("^",colnames(clean.all)[[k]],"$"),colnames(clean)))
          clean <- clean[,-index]
        }
      }
    }
    
    df.clean <- as.data.frame(clean)
    
    # Fit a logistic regression model using Maximum Likelihood Estimator
    mle.fit <- glm(y1 ~., data=df.clean, family = binomial)
    
    # Get the summary statistics for the model
    summary <- summary(mle.fit)$coefficients
    
    # 4. Obtain the raw p-values for the predictors selected with the MLE
    for (j in 1:nrow(raw)){
      for (k in 1:ncol(df.clean)){
        if (raw[[j,1]] == colnames(df.clean)[[k]]){
          
          for (l in 1:nrow(summary)){
            if (rownames(summary)[l] == colnames(df.clean)[[k]]){
              raw[[j,2]] = summary[l,4] # raw p-values
              raw[[j,3]] = summary[l,1] # raw estimate
              raw[[j,4]] = summary[l,2] # raw standard error
            }
          }
        }
      }
    }
    
    
    # 5. Done already
    # 6. Calculate the final p-value
    p.corr = c()
    for (j in 1:nrow(raw)){
      p.corr = c(p.corr,min(raw[j,2] * (ncol(clean)-1),1))
    }
    
    # Vectors containing all p-values, estimates, and standard errors
    p.final <- cbind(p.final,p.corr)
    e.final <- cbind(e.final, raw[,3])
    s.final <- cbind(s.final, raw[,4])
    
    
    # calculate the VIM  
    for (j in 1:nrow(coefs)){
      for (k in 1:nrow(VIM)){
        
        # Calculate the variable importance metric
        if (coefs[j,1] == VIM[k,1]){
          VIM[k,2] <- VIM[k,2] + 1/n
        }
      }
    }
  }
  
  # Calculate the p-values, estimates, and standard errors
  df.summary <- data.frame("Variables" = var.names, "p" = 1)
  
  for (i in 1:nrow(p.final)){
    # 1. obtain the empirical quantile function
    p.sorted <- sort(t(p.final[i,-1]), decreasing = FALSE) 
    deltas <- seq(0.05,1,0.01)
    quantiles <- quantile(p.sorted,probs = deltas)
    quantiles.df <- data.frame("d" = as.numeric(deltas), "q" = as.numeric(as.matrix(quantiles)), "q/d" = as.numeric(as.matrix(quantiles))/as.numeric(deltas))
    
    # 2. find the delta that minimizes q/d
    min.delta <- quantiles.df[match(min(quantiles.df$q.d),quantiles.df$q.d),1]
    q.min.delta <- quantiles.df[match(min(quantiles.df$q.d),quantiles.df$q.d),2]
    
    # 3. get the quantile that yelds the summary p-value
    summ.delta <- min(4 * q.min.delta/min.delta,1)
    summ.index <- match.closest(summ.delta,p.sorted,tolerance = Inf,nomatch = 0) + 1 # index of the summary p-value in the unsorted p-value list 
    p.summ <- p.final[[i,summ.index]] # need +1 here cause I removed one column earlier to sort
    e.summ <- e.final[[i,summ.index]]
    s.summ <- s.final[[i,summ.index]]
    df.summary$p[i] <- p.summ # summary p-value
  }
  
  # Results from the multi-split algorithm
  results <- cbind(VIM,df.summary[,-1])
  names(results)[2] <- "VIM"
  names(results)[3] <- "pValue"
  results$Code <- 0
  results$Check <- 0
  
  for (i in 1:nrow(results)){
    
    if(results[[i,1]]=="Dmg_HQYes"){
      results[[i,4]] <- "D1"
    }
    else if(results[[i,1]]=="Dmg_ProdPlantsYes"){
      results[[i,4]] <- "D2"
    }
    else if(results[[i,1]]=="Dmg_SalesOfficeYes"){
      results[[i,4]] <- "D3"
    }
    else if(results[[i,1]]=="Dmg_WarehouseYes"){
      results[[i,4]] <- "D4"
    }
    else if(results[[i,1]]=="Dmg_NoneYes"){
      results[[i,4]] <- "D5"
    }
    
    
    else if(results[[i,1]]=="SizeLarge"){
      results[[i,4]] <- "S1"
    }
    else if(results[[i,1]]=="SizeVeryLarge"){
      results[[i,4]] <- "S2"
    }
    
    else if(results[[i,1]]=="IndustryFinance_Insurance_Real_state"){
      results[[i,4]] <- "I1"
    }
    else if(results[[i,1]]=="IndustryWholesale_Retail"){
      results[[i,4]] <- "I2"
    }
    else if(results[[i,1]]=="IndustryServices"){
      results[[i,4]] <- "I3"
    }
    else if(results[[i,1]]=="IndustryManufacturing_Construction_Contracting"){
      results[[i,4]] <- "I4"
    }
    
    
    else if(results[[i,1]]=="BCP_CreationYes"){
      results[[i,4]] <- "P1"
    }
    else if(results[[i,1]]=="Action_PlanYes"){
      results[[i,4]] <- "P2"
    }
    else if(results[[i,1]]=="BCP_TrainingYes"){
      results[[i,4]] <- "P3"
    }
    else if(results[[i,1]]=="Alternative_CostumersYes"){
      results[[i,4]] <- "P4"
    }
    else if(results[[i,1]]=="CooperationYes"){
      results[[i,4]] <- "P5"
    }
    else if(results[[i,1]]=="Data_BackupYes"){
      results[[i,4]] <- "P6"
    }
    else if(results[[i,1]]=="OfficeYes"){
      results[[i,4]] <- "P7"
    }
    else if(results[[i,1]]=="Diversification_SuppliersYes"){
      results[[i,4]] <- "P8"
    }
    else if(results[[i,1]]=="NonSeismicYes"){
      results[[i,4]] <- "P9"
    }
    else if(results[[i,1]]=="Fixing_EquipmentYes"){
      results[[i,4]] <- "P10"
    }
    else if(results[[i,1]]=="SeismicYes"){
      results[[i,4]] <- "P11"
    }
    
    
    else if(results[[i,1]]=="Shortage_ElectricityYes"){
      results[[i,4]] <- "U1"
    }
    else if(results[[i,1]]=="Shortage_IndustrialWaterYes"){
      results[[i,4]] <- "U2"
    }
    else if(results[[i,1]]=="Shortage_WaterSupplyYes"){
      results[[i,4]] <- "U3"
    }
    else if(results[[i,1]]=="Shortage_SewageYes"){
      results[[i,4]] <- "U4"
    }
    else if(results[[i,1]]=="Shortage_GasYes"){
      results[[i,4]] <- "U5"
    }
    else if(results[[i,1]]=="Shortage_InformationYes"){
      results[[i,4]] <- "U6"
    }
    
    else if(results[[i,1]]=="Financing_LoanYes"){
      results[[i,4]] <- "F1"
    }
    else if(results[[i,1]]=="Financing_Self_fundYes"){
      results[[i,4]] <- "F2"
    }
    else if(results[[i,1]]=="Financing_InsuranceYes"){
      results[[i,4]] <- "F3"
    }
    else if(results[[i,1]]=="Financing_None"){
      results[[i,4]] <- "F4"
    }
  }
  
  
  for (i in 1:nrow(results)){
    results$Check[i] <- 0
    if (results[[i,2]] > 0.75 & results[[i,3]] < 0.05){
      results$Check[i] <- 1
    }
  }
  
  
  
  ##############################################
  # Prediction analysis
  ##############################################
  # Use the predictors found as important to predict the testing data
  # Get hold of the test data
  test.factors <- as.matrix(model.matrix(~.,test)[,-1])
  train.factors <- as.matrix(model.matrix(~.,smoted)[,-1])
  n <- ncol(train.factors)
  
  for (i in 1:(n-1)){
    ii = n-i 
    if ((results[[ii,2]] < 0.75 | results[[ii,3]] > 0.05)){
      #cat(ii,colnames(test.factors)[ii],"\n")
      train.factors <- train.factors[,-ii]
      test.factors <- test.factors[,-ii]
    }
  }

  
  # Fit the final MLE model
  cat("Final Model", "\n")
  final.fit <<- glm(y1 ~., data=as.data.frame(train.factors), family = binomial)
  final.coefficients <- coef(final.fit)
  
  df.test.factors <- as.data.frame(test.factors)
  x.test <- select(df.test.factors,-ncol(df.test.factors))
  predicted <- predict(final.fit, newdata = x.test, type = "response")
  

  
  ##############################################
  # Get the area under the curve
  ##############################################
  y.test <- test.factors[,ncol(test.factors)]
  plot.roc(y.test,predicted)
  theROC <<- roc(y.test,predicted)
  theAUC <<- auc(y.test,predicted)
  cat("The AUC is:",theAUC, "\n")
  

  
  ##############################################
  # Plot the predicted vs the observed values
  ##############################################
  df.S <- data.frame("XB" = 0, "Predicted" = predicted, "True" = y.test)
  
  for (i in 1:nrow(x.test)){
    sum = 0
    for (j in 1:ncol(x.test)){
      sum = sum + final.coefficients[j+1] * x.test[[i,j]]
    }
    df.S[[i,1]] <- sum
  }
  
  unique.XB <- unique(df.S$XB)
  
  # Plot the observed values that are 1 
  df.plot1 <- data.frame("XB" = unique.XB, "Value" = 1, "size" = 0)
  for (i in 1:nrow(df.S)){
    for (j in 1:nrow(df.plot1)){
      df.plot1[[j,3]] <- sum(df.S$True == 1 & df.S$XB == df.plot1[[j,1]])
    }
  }
  
  # Plot the observed values that are 0 
  df.plot0 <- data.frame("XB" = unique.XB, "Value" = 0, "size" = 0)
  for (i in 1:nrow(df.S)){
    for (j in 1:nrow(df.plot0)){
      df.plot0[[j,3]] <- sum(df.S$True == 0 & df.S$XB == df.plot0[[j,1]])
    }
  }
  
  # Plot the predicted values 
  df.plot.pred <- data.frame("XB" = unique.XB, "Value" = 0, "size" = 5)
  for (i in 1:nrow(df.S)){
    for (j in 1:nrow(df.plot0)){
      df.plot0[[j,3]] <- sum(df.S$True == 0 & df.S$XB == df.plot0[[j,1]])
    }
  }
  
  themax <- ceiling(max(df.S$XB))
  themin <- floor(min(df.S$XB))
  thedelta <- (themax - themin)/4
  breaks <- c(themin, themin+thedelta, themin+2*thedelta, themin+3*thedelta, themax)
  labels <- sprintf("%.2f",breaks)
  
  
  print( 
    ggplot() +
      geom_point(data = df.plot0, aes(x=XB, y=Value, size = size)) +
      geom_point(data = df.plot1, aes(x=XB, y=Value, size = size)) +
      geom_point(data = df.S, aes(x=XB, y=Predicted, size = 20, col = "red")) +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
      scale_x_continuous(breaks=breaks, labels=labels, name="X\u03b2", limits=c(themin, themax)) +
      scale_y_continuous(name = "P(R<t),observed values") + 
      #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      theme_bw() 
  )
  
  print(
  ggplot(results, aes(x=VIM, y=pValue)) +
    geom_point(col = ifelse(results$Check > 0,"red","black"),size=5) + 
    xlab('VIM') +
    ylab('p-values') + 
    geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
    geom_vline(xintercept=0.75, linetype="dashed", color = "red") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    theme_bw()
  )
  
  cat ("Analysis completed!", "\n")
  return(results)
}





