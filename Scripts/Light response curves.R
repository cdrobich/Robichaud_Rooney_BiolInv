library(tidyverse)

ciras <- read.csv("data/CIRAS_sum_bothyears.csv", header = TRUE)

#### Light Response Curves from Heberling #####
# Heberling JM, Brouwer NL, Kalisz S. 2017. Effects of deer on the photosynthetic performance of invasive and native forest
# herbs. AoB PLANTS 9: plx011; doi:10.1093/aobpla/plx011

# 'leaf level photosynthesis' https://sites.google.com/site/fridleylab/home/protocols for code"

lrc <- read.csv("data/sample_lrc.txt",sep="",skip=16)

PARlrc <- lrc$PARi #PAR (aka PPFD or Q)
photolrc <- lrc$Photo #net photosynthetic rate (Anet)

curvelrc <- data.frame(PARlrc,photolrc)
curvelrc # *inspect raw data and check notebook (data reasonable or need edited/discarded?)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PARlrc,photolrc,xlab="", ylab="", ylim=c(-2,max(photolrc)+2),cex.lab=1.2,cex.axis=1.5,cex=2)
mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)

curve.nlslrc = nls(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=1)) 


par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PARlrc,photolrc,xlab="", ylab="", ylim=c(-2,max(photolrc)+2),cex.lab=1.2,cex.axis=1.5,cex=2)
mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=2)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2,cex=2)
curve((1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)

# ---Solve for light compensation point (LCPT), PPFD where Anet=0 ---
x <-function(x) {(1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1]}

uniroot(x,c(0,50),extendInt = "yes")$root #Light compensation point

# ---Solve for light saturation point (LSP), PPFD where 75% of Amax is achieved (75% is arbitrary - cutoff could be changed)
y <-function(y) {(1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*y+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*y+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*y))-summary(curve.nlslrc)$coef[3,1]-(0.75*summary(curve.nlslrc)$coef[1,1])+0.75*(summary(curve.nlslrc)$coef[3,1])}

uniroot(y,c(0,1000),extendInt = "yes")$root #Light saturation point




# My data -----------------------------------------------------------------

comp_lrc <- ciras %>% 
  select(ID_yr, light, carbon) %>% 
  pivot_wider(names_from = ID_yr, values_from = carbon)

PARlrc <- comp_lrc$light

com_lrc <- comp_lrc[,-1]


glimpse(com_lrc)


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


write.csv(results, "Data/lrc_results.csv")



# Import LRC values -------------------------------------------------------

lrc_data <- read.csv("Data/lrc_results_env.csv")

sum.LCPT <- lrc_data %>% 
  group_by(Species, Treatment) %>% 
  summarise(median = median(LCPT, na.rm = TRUE),
            mean = mean(LCPT, na.rm = TRUE),
            LCPT.sd = sd(LCPT, na.rm = TRUE),
            N = length(LCPT),
            sterr = (LCPT.sd/sqrt(N)))


phrag %>% 
  group_by(type, Treatment) %>% 
  summarise(median = median(LCPT, na.rm = TRUE),
            mean = mean(LCPT, na.rm = TRUE),
            LCPT.sd = sd(LCPT, na.rm = TRUE),
            N = length(LCPT),
            sterr = (LCPT.sd/sqrt(N)))

##  type                     Treatment      median  mean LCPT.sd     N sterr
##1 Phragmites_Calamagrostis Competition      30.6 17.6     46.9    12 13.5 
##2 Phragmites_Calamagrostis No competition   35.4 35.9     25.8    12  7.44
##3 Phragmites_Carex         Competition      38.4 31.5     23.8    12  6.87
##4 Phragmites_Carex         No competition   33.5 35.8     17.2    11  5.18
##5 Phragmites_Typha         Competition      11.1 -9.83    63.9    12 18.5 
##6 Phragmites_Typha         No competition   29.5 18.2     23.8    10  7.52


sum.LSP <- lrc_data %>% 
  group_by(Species, Treatment) %>% 
  summarise(median = median(LSP, na.rm = TRUE),
            mean = mean(LSP, na.rm = TRUE),
            LSP.sd = sd(LSP, na.rm = TRUE),
            N = length(LSP),
            sterr = (LSP.sd/sqrt(N)))


phrag %>% 
  group_by(type, Treatment) %>% 
  summarise(median = median(LSP, na.rm = TRUE),
            mean = mean(LSP, na.rm = TRUE),
            LSP.sd = sd(LSP, na.rm = TRUE),
            N = length(LSP),
            sterr = (LSP.sd/sqrt(N)))

#type                     Treatment      median  mean LSP.sd     N sterr
#1 Phragmites_Calamagrostis Competition      870.  861.   300.    12  86.5
#2 Phragmites_Calamagrostis No competition   827.  888.   277.    12  80.1
#3 Phragmites_Carex         Competition      766.  748.   139.    12  40.2
#4 Phragmites_Carex         No competition   886.  855.   172.    11  51.9
#5 Phragmites_Typha         Competition      870.  835.   332.    12  95.8
#6 Phragmites_Typha         No competition  1068.  980.   308.    10  97.3


LCPT <- ggplot(lrc_data, aes(x = Treatment,
                             y = LCPT)) +
  geom_boxplot(alpha = 0.7, lwd = 1, width = 0.5) +
  geom_jitter(aes(fill = Treatment,
                  shape = Treatment),
              size = 5,
              width = 0.08,
              alpha = 0.7) +
  facet_wrap("Species", scales = "free") + 
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=9), 
        legend.text=element_text(size=9),
        strip.text = element_text(size=12),
        axis.text.x = element_text(size = 12)) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(y = expression(paste("Light Compensation Point"," ", " (", "\u00B5mol photons ", " ", m^-2," ", s^-1, sep=")")),
       x = " ") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#24908C", "#3A518B")) +
  scale_shape_manual(values = c(21,24))




LSP <- ggplot(lrc_data, aes(x = Treatment,
                            y = LSP)) +
  geom_boxplot(alpha = 0.7, lwd = 1, width = 0.5) +
  geom_jitter(aes(fill = Treatment,
                  shape = Treatment),
              size = 5,
              width = 0.08,
              alpha = 0.7) +
  facet_wrap("Species", scales = "free") + 
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=9), 
        legend.text=element_text(size=9),
        strip.text = element_text(size=12),
        axis.text.x = element_text(size = 12)) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(y = expression(paste("Light Saturation Point"," ", " (", "\u00B5mol photons ", " ", m^-2," ", s^-1, sep=")")),
       x = " ") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#24908C", "#3A518B")) +
  scale_shape_manual(values = c(21,24))


review.plot <- LCPT + LSP + plot_annotation(tag_levels = "A")

ggsave("Figures/LSP_LCPT.jpeg", dpi = 300)


# anova -------------------------------------------------------------------

library(car)
library(agricolae)
library(performance)

lrc_data <- lrc_data %>% 
  unite("type", Species,Neighbour, remove = FALSE)


LSP_aov <- lm(LSP ~ type * Treatment, data = lrc_data, na.rm = TRUE)

Anova(LSP_aov, type = 3)

#Response: LSP
#                   Sum Sq  Df F value Pr(>F)
#(Intercept)    1.2181e+07   1  0.0015 0.9690
#type           1.3365e+06   5  0.0000 1.0000
#Treatment      9.3600e+03   1  0.0000 0.9991
#type:Treatment 4.0479e+10   5  1.0095 0.4143
#Residuals      1.1628e+12 145 

check_model(LSP_aov)


lrc_data <- lrc_data %>% mutate(logLCPT = log10(LCPT + 1))

LCPT_aov <- lm(logLCPT ~ type * Treatment, data = lrc_data, na.rm = TRUE)

Anova(LCPT_aov , type = 3)

#Response: LCPT
#               Sum Sq  Df F value  Pr(>F)   
#(Intercept)     13334   1  9.9062 0.00200 **
#type            14154   5  2.1030 0.06832 . 
#Treatment           0   1  0.0002 0.98840   
#type:Treatment  11954   5  1.7762 0.12130   
#Residuals      195178 145 


check_model(LCPT_aov)

res <- c("Carex", "Calamagrostis", "Typha")
resident <- lrc_data %>% filter(Species %in% res)
phrag <- lrc_data %>% filter(Species == "Phragmites")


LCPT_aovres <- lm(LCPT ~ type * Treatment, data = resident, na.rm = TRUE)

Anova(LCPT_aovres , type = 3)

#Response: LCPT
#                Sum Sq Df F value   Pr(>F)   
#(Intercept)     13334  1 10.4681 0.001752 **
#type              323  2  0.1268 0.881080   
#Treatment           0  1  0.0002 0.988089   
#type:Treatment   5781  2  2.2694 0.109832   
#Residuals      104451 82


LSP_aovres <- lm(LSP ~ type * Treatment, data = resident, na.rm = TRUE)

Anova(LSP_aovres , type = 3)

#Response: LSP
#Sum Sq Df F value Pr(>F)
#(Intercept)    1.2181e+07  1  0.0009 0.9767
#type           4.5375e+05  2  0.0000 1.0000
#Treatment      9.3600e+03  1  0.0000 0.9994
#type:Treatment 3.4073e+10  2  1.2014 0.3060
#Residuals      1.1628e+12 82 


# Phragmites
LCPT_aovph <- lm(LCPT ~ type * Treatment, data = phrag, na.rm = TRUE)

Anova(LCPT_aovph , type = 3)

#Response: LCPT
#Sum Sq Df F value  Pr(>F)  
#(Intercept)      3719  1  2.5822 0.11307  
#type            10629  2  3.6904 0.03051 *
#Treatment        2002  1  1.3904 0.24277  
#type:Treatment   1607  2  0.5578 0.57528  
#Residuals       90727 63  

type <- HSD.test(LCPT_aovph, "type")

#LCPT groups
#Phragmites_Carex         33.552141      a
#Phragmites_Calamagrostis 26.737704     ab
#Phragmites_Typha          2.907755      b


LSP_aovph <- lm(LSP ~ type * Treatment, data = phrag, na.rm = TRUE)

Anova(LSP_aovph , type = 3)

#response: LSP
#Sum Sq Df  F value Pr(>F)    
#(Intercept)    8901834  1 127.2261 <2e-16 ***
#  type             84561  2   0.6043 0.5496    
#Treatment         4295  1   0.0614 0.8051    
#type:Treatment   42054  2   0.3005 0.7415    
#Residuals      4408024 63