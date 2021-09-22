
# load libraries ----------------------------------------------------------

library(tidyverse)

# load data ---------------------------------------------------------------

com_lrc <- read.csv("Data/organized_ciras_lrc.csv")

PARlrc <- c(1500, 1000, 500, 200, 100, 50, 0)


# run LRC through colums --------------------------------------------------

results <-
  do.call(rbind,
          map(.x = colnames(com_lrc), ~{
            
            photo <- com_lrc%>% select(!!sym(.x))%>% purrr::reduce(c)
            
            model <-
              tryCatch(nls(photo ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photo)-min(photo)),
                                                                                                                             AQY=0.05,Rd=-min(photo),theta=1)), error = function(e){
                                                                                                                               NULL
                                                                                                                             })
            
            LCPT <-function(x)
            {
              (1/(2*summary(model)$coef[4,1]))*(summary(model)$coef[2,1]*x+summary(model)$coef[1,1]-sqrt((summary(model)$coef[2,1]*x+summary(model)$coef[1,1])^2-4*summary(model)$coef[2,1]*summary(model)$coef[4,1]*summary(model)$coef[1,1]*x))-summary(model)$coef[3,1]}
            
            
            LSP <-function(y) {(1/(2*summary(model)$coef[4,1]))*(summary(model)$coef[2,1]*y+summary(model)$coef[1,1]-sqrt((summary(model)$coef[2,1]*y+summary(model)$coef[1,1])^2-4*summary(model)$coef[2,1]*summary(model)$coef[4,1]*summary(model)$coef[1,1]*y))-summary(model)$coef[3,1]-(0.75*summary(model)$coef[1,1])+0.75*(summary(model)$coef[3,1])}
            
            
            if(is.null(model)){
              
              results <- data.frame(Variable = .x,
                                    LCPT = "NA",
                                    LSP = "NA")
              
            }
            
            
            else{
              
              results <- data.frame(Variable = .x,
                                    LCPT = uniroot(LCPT,c(0,50), extendInt = "yes")$root,
                                    LSP = uniroot(LSP,c(0,1000), extendInt = "yes")$root)
              
              
            }
            
            
          }
          )
  )


results



write.csv(results, "Data/lrc_results.csv") # save results
