#Silver shiner bioenergetics analysis
#Jake Brownscombe


#source packages
setwd("~/github/Silver_Shiner_Bioenergetics")
source("worksheets/packages.R")




#metabolism data ----
#Field MMR 
MMR.field <- readRDS("data/MMR.field.rds") %>% as.data.frame()
MMR.field$V.ratio <- MMR.field$Volume/(MMR.field$weight.g/1000)
head(MMR.field)


#Chase time data 
ggplot(MMR.field, aes(TempMean, chase.sec/60, col=weight.g))+geom_point()+
  theme_bw()+scale_x_continuous(limits=c(0,22))+scale_y_continuous(limits=c(0,20))+
  scale_color_viridis_c(name="Weight (g)")+xlab("Temperature (ºC)")+ylab("Chase time (min)")+
  labs(title="Chase Time", subtitle="A)Temperature - Chase time")+
  
  ggplot(MMR.field, aes(weight.g, chase.sec/60,col=TempMean))+geom_point()+
  theme_bw()+scale_x_continuous(limits=c(0,15))+scale_y_continuous(limits=c(0,20))+
  scale_color_viridis_c(option="plasma", name="ºC")+xlab("Weight (g)")+ylab("Chase time (min)")+
  labs(subtitle="B) Weight - Chase time")+
  
  ggplot(MMR.field, aes(chase.sec/60, O2.Kg.L.h, col=weight.g))+geom_point()+
  theme_bw()+scale_x_continuous(limits=c(0,22))+scale_y_continuous(limits=c(0,500))+
  scale_color_viridis_c(name="Weight (g)")+ylab(bquote(~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))+xlab("Chase time (min)")+
  labs(subtitle="C) Weight - Metabolic Rate")+
  
  ggplot(MMR.field, aes(chase.sec/60, O2.Kg.L.h, col=TempMean))+geom_point()+
  theme_bw()+scale_x_continuous(limits=c(0,22))+scale_y_continuous(limits=c(0,500))+
  scale_color_viridis_c(option="plasma",name="ºC")+ylab(bquote(~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))+xlab("Chase time (min)")+
  labs(subtitle="D) Temperature - Metabolic Rate")
#appendix Fig


#model
#chase response
chase.lm <- lm(chase.sec ~ poly(TempMean,3) + weight.g, data=MMR.field)
summary(chase.lm)
drop1(chase.lm, test="Chisq")
#complex third order effect 



#mo2 response 
chase.mo2.lm <- lm(log(O2.Kg.L.h) ~ chase.sec + poly(TempMean, 2) + log(weight.g) + 
                     chase.sec*poly(TempMean, 2) + chase.sec*log(weight.g), data=MMR.field)
summary(chase.mo2.lm)
drop1(chase.mo2.lm, test="Chisq")

chase.mo2.lm <- lm(log(O2.Kg.L.h) ~ chase.sec + poly(TempMean, 2) + log(weight.g) + 
                     chase.sec*poly(TempMean, 2), data=MMR.field)
summary(chase.mo2.lm)
drop1(chase.mo2.lm, test="Chisq")
#interaction between chase time and poly temp 


#plot field MMR
ggplot(MMR.field, aes(TempMean, O2.Kg.L.h, col=weight.g))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2))+
  theme_bw()+scale_x_continuous(limits=c(0,22))+scale_y_continuous(limits=c(0,700))+
  scale_color_viridis_c(name="Weight (g)")+xlab("Temperature (ºC)")+
  ylab(bquote(~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))
#appendix Fig









#Field SMR 
SMR.field <- readRDS("data/smr.field.rds") %>% as.data.frame()
head(SMR.field)

#plot
ggplot(SMR.field, aes(Rest.Duration, O2.Kg.L.h, col=TempMean))+geom_point()+facet_wrap(~SS.ID)+
  labs(x="Time in chamber (hours)", y=bquote(~MO[2]* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))+
  scale_colour_viridis_c(option="plasma", name="ºC")+theme_bw()
#appendix Fig

#filter to trials with at least 3 data points
SMR.field.f <- SMR.field %>% filter(sample.n>=3)
ggplot(SMR.field.f, aes(Rest.Duration, O2.Kg.L.h, col=TempMean))+geom_point()+facet_wrap(~SS.ID)+
  scale_colour_viridis_c(option="plasma", name="ºC")+theme_bw()

MR.field.sum <- SMR.field.f %>% group_by(SS.ID) %>%  slice_min(order_by =O2.Kg.L.h, n = 3) %>%
  summarize(temp=mean(TempMean), weight.g=mean(weight.g), SMR = mean(O2.Kg.L.h)) %>% filter(!is.na(SMR))

MR.field.sum2 <- SMR.field.f  %>% group_by(SS.ID) %>% 
  summarize(RMR = mean(O2.Kg.L.h), MMR=max(O2.Kg.L.h)) 
head(MR.field.sum2)

MR.field.sum <- cbind(MR.field.sum, MR.field.sum2 %>% select(-SS.ID))
head(MR.field.sum)


ggplot(MR.field.sum, aes(temp, SMR))+geom_point()+geom_smooth(method="lm", col="black")+
  geom_point(aes(temp, RMR), col="blue")+
  geom_smooth(aes(temp, RMR), method="lm", col="blue")+
  geom_point(aes(temp, MMR), col="red")+
  geom_smooth(aes(temp, MMR), method="lm", col="red")
#MMR from SMR trials is clearly not useful. revisit in combintation with lab SMR below 









#Lab SMR 
SMR.lab <- readRDS("data/smr.lab.rds") %>% as.data.frame()
head(SMR.lab)

ggplot(SMR.lab, aes(Rest.Duration, O2.Kg.L.h, col=temp))+geom_point()+facet_wrap(~SS.ID, ncol=7)+
  theme_bw()+scale_color_viridis_c(option="plasma")+labs(x="Hour", y="MO2")
#appendix fig



#Some 24 and 48 hour trials. Compare to see if this is an issue
head(SMR.lab)
SMR.lab.48 <- SMR.lab %>% filter(Rest.Duration>40)
`%ni%` <- Negate(`%in%`)
SMR.lab.48 <- SMR.lab %>% filter(SS.ID %in% SMR.lab.48$SS.ID)
SMR.lab.24 <- SMR.lab %>% filter(SS.ID %ni% SMR.lab.48$SS.ID)




#48 hour 10 lowest after 5 removed 
SMR.lab.48$nrow <- 1:length(SMR.lab.48$SS.ID)
SMR.new.5low <- SMR.lab.48 %>% 
  filter(!is.na(O2.Kg.L.h) & O2.Kg.L.h>0) %>% group_by(SS.ID) %>%  slice_min(order_by=O2.Kg.L.h, n = 5)
SMR.lab.482 <- SMR.lab.48 %>% filter(nrow %ni% SMR.new.5low$nrow)

#10 lowest SMR
SMR.48.sum <- SMR.lab.482 %>% group_by(SS.ID) %>%  slice_min(order_by =O2.Kg.L.h, n = 10) %>%
  summarize(temp=mean(temp),SMR = mean(O2.Kg.L.h), weight.g=mean(weight)) %>% filter(!is.na(SMR)) %>% as.data.frame()
head(SMR.48.sum)


#24 hour 10 lowest after 5 removed 
SMR.lab.24$nrow <- 1:length(SMR.lab.24$SS.ID)
SMR.new.5low <- SMR.lab.24 %>% 
  filter(!is.na(O2.Kg.L.h) & O2.Kg.L.h>0) %>% group_by(SS.ID) %>%  slice_min(order_by=O2.Kg.L.h, n = 5)
SMR.lab.242 <- SMR.lab.24 %>% filter(nrow %ni% SMR.new.5low$nrow)

#10 lowest SMR
SMR.24.sum <- SMR.lab.242 %>% group_by(SS.ID) %>%  slice_min(order_by =O2.Kg.L.h, n = 10) %>%
  summarize(temp=mean(temp),SMR = mean(O2.Kg.L.h), weight.g=mean(weight)) %>% filter(!is.na(SMR)) %>% as.data.frame()
head(SMR.24.sum)

ggplot(SMR.48.sum, aes(temp, SMR))+geom_point()+
  geom_point(data=SMR.24.sum, aes(temp, SMR), col="blue")


#looks the same. compare just at warmer temps where we have both. 
SMR.48.comp <- rbind(SMR.24.sum %>% mutate(type="24")%>% filter(temp>15 & temp<24),
                     SMR.48.sum %>% mutate(type="48") %>% filter(temp>15 & temp<24))
ggplot(SMR.48.comp, aes(temp, SMR, col=type))+geom_point()

head(SMR.48.comp)
SMR.48.comp$logSMR <- log(SMR.48.comp$SMR)
SMR.48.comp$logweight <- log(SMR.48.comp$weight.g)

#model
SMR48.lm <- lm(logSMR ~ type, data=SMR.48.comp)
drop1(SMR48.lm,test="Chisq")
summary(SMR48.lm)
#report





#summarize all lab MR data
head(SMR.lab)
#remove 5 lowest values 
SMR.lab$nrow <- 1:length(SMR.lab$SS.ID)
SMR.lab.5low <- SMR.lab %>% 
  filter(!is.na(O2.Kg.L.h) & O2.Kg.L.h>0) %>% group_by(SS.ID) %>%  slice_min(order_by=O2.Kg.L.h, n = 5)

SMR.lab2 <- SMR.lab %>% filter(nrow %ni% SMR.lab.5low$nrow)

#10 lowest SMR
MR.lab.sum <- SMR.lab2 %>% group_by(SS.ID) %>%  slice_min(order_by =O2.Kg.L.h, n = 10) %>%
  summarize(temp=mean(temp),SMR = mean(O2.Kg.L.h), weight.g=mean(weight)) %>% filter(!is.na(SMR)) %>% as.data.frame()
head(MR.lab.sum)

MR.lab.MMR <- SMR.lab %>% group_by(SS.ID) %>%  summarize(RMR=mean(O2.Kg.L.h), MMR = max(O2.Kg.L.h))
MR.lab.sum$MMR <- MR.lab.MMR$MMR[match(MR.lab.sum$SS.ID, MR.lab.MMR$SS.ID)]
MR.lab.sum$RMR <- MR.lab.MMR$RMR[match(MR.lab.sum$SS.ID, MR.lab.MMR$SS.ID)]
MR.lab.sum$chase.type <- SMR.lab$chase.type[match(MR.lab.sum$SS.ID, SMR.lab$SS.ID)]

head(MR.lab.sum)
ggplot(MR.lab.sum, aes(temp, SMR), col="black")+geom_point()+geom_smooth(method="lm",col="black")+
  geom_point(data=MR.lab.sum, aes(temp, RMR), col="blue")+
  geom_smooth(data=MR.lab.sum, aes(temp, RMR), col="blue", method="lm")+
  geom_point(data=MR.lab.sum, aes(temp, MMR), col="red")+
  geom_smooth(data=MR.lab.sum, aes(temp, MMR), col="red", method="lm")





#combine lab and field SMR data
MR.field.sum$data.type <- "field"
MR.field.sum$chase.type <- "none"
MR.lab.sum$data.type <- "lab"
SMR.comb <- rbind(MR.lab.sum %>% select(-RMR, -MMR), 
                  MR.field.sum %>% select(-RMR, -MMR)) #**problem here with MMR having values added to it???? 
head(SMR.comb)


#plot and model combined SMR 
SMR.comb$logSMR <- log(SMR.comb$SMR)
SMR.comb$logweight <- log(SMR.comb$weight.g)

ggplot(SMR.comb, aes(temp, logSMR, col=data.type))+geom_point()+geom_smooth(method="lm")+
  theme_bw()+scale_y_continuous(limits=c(2,6))+
  scale_color_viridis_d(end=0.8, name="Data type", option="plasma")+
  labs(x="Temperature (ºC)", y=bquote(~ln-SMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))



#Model to assess differences
SMR.LM <-lm(logSMR~ temp + logweight + temp*logweight + data.type, data = SMR.comb)
drop1(SMR.LM, test="Chisq")
summary(SMR.LM)
SMR.LM <-lm(logSMR~ temp + logweight, data = SMR.comb)
drop1(SMR.LM, test="Chisq")
summary(SMR.LM)



#predict
smr.preds <- merge(data.frame(weight=seq(1, 10, 0.1)), data.frame(temp= seq(1, 26, 0.1)))
smr.preds$logweight <- log(smr.preds$weight)
smr.preds$smr.pred <- exp(predict(SMR.LM, smr.preds))

#plot
ggplot(SMR.comb, aes(temp, logSMR, col=data.type))+geom_point()+geom_smooth(method="lm")+
  theme_bw()+scale_y_continuous(limits=c(2,6))+
  scale_color_viridis_d(end=0.8, name="Data type")+
  labs(x="Temperature (ºC)", y=bquote(~log-SMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))

ggplot(SMR.comb, aes(temp, SMR))+geom_point()+
  geom_smooth(data=smr.preds %>% filter(weight==3 | weight==6 | weight==9),
              aes(temp, smr.pred, col=as.factor(weight)), inherit.aes = F)+
  theme_bw()+scale_color_viridis_d(name="Weight (g)", end=0.9)+
  labs(x="Temperature", y=bquote(~SMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))
#














#MMR lab and field ----
head(MR.lab.sum)
head(MMR.field)

#errors with chase type assignment here, fix:
MR.lab.sum$chase.type <- ifelse(MR.lab.sum$chase.type=="1.min", "none", "1.min")
MR.lab.sum$chase.type <- ifelse(MR.lab.sum$temp>23, "none", MR.lab.sum$chase.type)

ggplot(MMR.field, aes(TempMean, O2.Kg.L.h))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2), col="black")+
  geom_point(data=MR.lab.sum, aes(temp, MMR, col=chase.type), inherit.aes = F)+
  geom_smooth(data=MR.lab.sum, aes(temp, MMR), col="blue", method="lm")+
  theme_bw()+scale_x_continuous(limits=c(0,25))+scale_y_continuous(limits=c(0,700))+
  scale_shape_discrete(name="Protocol")+
  #scale_color_viridis_c(name="Weight (g)")+xlab("Temperature (ºC)")+
  ylab(bquote(~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'))


#compare MMR at 18 between chase and non 
mmr18 <- MR.lab.sum %>% filter(temp>12 & temp<20)
ggplot(data=mmr18, aes(temp, MMR, col=chase.type), inherit.aes = F)+geom_point()
summary(lm(MMR~chase.type, data=mmr18))
#sig difference. Drop non chase here at 18 moving forward



mmr.f <- MR.lab.sum %>% filter(temp<12 | temp>20 | chase.type=="1.min")
ggplot(data=mmr.f, aes(temp, MMR, col=chase.type), inherit.aes = F)+geom_point()


#model
mmr.f$lMMR <- log(mmr.f$MMR)
mmr.f$logweight <- log(mmr.f$weight.g)

MMR.LM <-lm(MMR~ temp + logweight + temp*logweight, data = mmr.f)
drop1(MMR.LM, test="Chisq")

MMR.LM <-lm(MMR~ temp + logweight, data = mmr.f)
drop1(MMR.LM, test="Chisq")
summary(MMR.LM)


#preds
head(smr.preds)
smr.preds$mmr.pred <- predict(MMR.LM, smr.preds)


ggplot(mmr.f, aes(temp, log(MMR)))+geom_point()+
  geom_smooth(data=smr.preds %>% filter(weight==3 | weight==6 | weight==9),
              aes(temp, log(mmr.pred), col=as.factor(weight)), inherit.aes = F)+
  theme_bw()+scale_color_viridis_d(name="Weight (g)", end=0.9)+
  labs(x="Temperature", y=bquote(~ln~MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), title="B) MMR - fitted")+
  scale_y_continuous(limits=c(0,7))




#Aerobic scope ----
head(smr.preds)
smr.preds$AS <- smr.preds$mmr.pred-smr.preds$smr.pred
smr.preds$AS <- ifelse(smr.preds$AS<0,0, smr.preds$AS)
ggplot(smr.preds %>% filter(weight==3 | weight==6 | weight==9),
       aes(temp, AS, col=as.factor(weight)))+geom_smooth()+
  theme_bw()+scale_y_continuous(limits=c(0,500))+
  scale_color_viridis_d(end=0.8, name="Weight (g)")+
  labs(x="Temperature (ºC)", y=bquote(~AS* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'),title="E) AS - Fish Weight")+
  
  ggplot(smr.preds, aes(temp, weight, fill=AS))+geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(option="magma",name="AS")+
  labs(x="Temperature (ºC)", y="Fish Weight (g)",title="F) AS - Temperature & Weight")


#SMR, MMR, AS on same plots 
head(smr.preds)
smr.preds %>% filter(AS==max(AS))
ggplot(smr.preds %>% filter(weight==5), aes(temp, smr.pred), col="blue")+geom_smooth()+
  geom_smooth(data=smr.preds %>% filter(weight==5), aes(temp, mmr.pred), col="red")+
  geom_smooth(data=smr.preds %>% filter(weight==5), aes(temp, AS, col=as.factor(weight)), linetype=2, col="green", fill="NA")+
  scale_color_viridis_d(end=0.8, name="Weight (g)")+theme_bw()+
  labs(x="Temperature (ºC)", y=bquote(~MR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'),title="G) Metabolic Metrics (5 g fish)")+
  scale_x_continuous(limits=c(0,30))+scale_y_continuous(limits=c(0,650))







#Data type plots
ggplot(SMR.comb, aes(temp, logSMR, col=data.type))+geom_point()+geom_smooth(method="lm")+
  theme_bw()+scale_y_continuous(limits=c(0,7))+
  scale_colour_manual(values = c("black", "blue"), name="Data type")+
  labs(x="Temperature (ºC)", y=bquote(~ln-SMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'),
       title="Field & Lab Data", subtitle="A) Standard Metabolic Rate")+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(fill = "white", colour = "black"))+
  
  
  ggplot(MMR.field, aes(TempMean, log(O2.Kg.L.h)))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2), col="black")+
  geom_point(data=MR.lab.sum, aes(temp, log(MMR), pch=chase.type), col="blue")+
  geom_smooth(data=MR.lab.sum, aes(temp, log(MMR)), col="blue", method="lm")+
  annotate("label", x=23, y=6.9, label="Lab", col="blue")+
  annotate("label", x=21.2, y=5.2, label="Field", col="black")+
  theme_bw()+scale_x_continuous(limits=c(0,25))+scale_y_continuous(limits=c(4,7))+
  scale_shape_manual(name="Lab Protocol", values=c(17,19))+
  scale_color_viridis_c(name="Weight (g)")+
  labs(x="Temperature (ºC)", y=(bquote(~ln - MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')')), 
       subtitle="B) Maximum Metabolic Rate")+
  theme(legend.position = c(0.85, 0.15),legend.background = element_rect(fill = "white", colour = "black"))+
  
  plot_layout(ncol=1)
#paper Fig


#smr
ggplot(SMR.comb, aes(temp, log(SMR)))+geom_point()+
  geom_smooth(data=smr.preds %>% filter(weight==3 | weight==6 | weight==9),
              aes(temp, log(smr.pred), col=as.factor(weight)), inherit.aes = F, fill=NA)+
  theme_bw()+scale_color_viridis_d(name="Weight (g)", end=0.9)+
  scale_y_continuous(limits=c(0,8))+
  labs(x="Temperature", y=bquote(~ln-SMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), title="A) SMR - fitted")+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(fill = "white", colour = "black"))+
  
  ggplot(smr.preds, aes(temp, weight, fill=smr.pred))+geom_raster()+
  scale_fill_viridis_c(option="magma", name="SMR", limits=c(0,750))+
  theme_classic()+labs(x="Temperature (ºC)", y="Weight (g)", title='B) SMR - Temperature & Weight')+
  
  #mmr
  ggplot(mmr.f, aes(temp, log(MMR)))+geom_point()+
  geom_smooth(data=smr.preds %>% filter(weight==3 | weight==6 | weight==9),
              aes(temp, log(mmr.pred), col=as.factor(weight)), inherit.aes = F, fill=NA)+
  theme_bw()+scale_color_viridis_d(name="Weight (g)", end=0.9)+
  scale_y_continuous(limits=c(3,8))+
  labs(x="Temperature", y=bquote(~ln-MMR* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'), title="C) MMR - fitted")+
  theme(legend.position = c(0.85, 0.2),legend.background = element_rect(fill = "white", colour = "black"))+
  
  ggplot(smr.preds, aes(temp, weight, fill=mmr.pred))+geom_raster()+
  scale_fill_viridis_c(option="magma", name="MMR", limits=c(0,750))+
  theme_classic()+labs(x="Temperature (ºC)", y="Weight (g)", title='D) MMR - Temperature & Weight')+
  
  
  #AS
  ggplot(smr.preds %>% filter(weight==5), aes(temp, smr.pred))+geom_smooth(col="black")+
  geom_smooth(data=smr.preds %>% filter(weight==5), aes(temp, mmr.pred), col="red")+
  geom_smooth(data=smr.preds %>% filter(weight==5), aes(temp, AS, col=as.factor(weight)), linetype=2, col="green", fill="NA")+
  scale_color_viridis_d(end=0.8, name="Weight (g)")+theme_bw()+
  annotate("label", x=10, y=460, label="MMR", col="red")+
  annotate("label", x=20, y=100, label="SMR", col="black")+
  annotate("label", x=20, y=440, label="AS", col="green")+
  scale_y_continuous(limits=c(0,700))+
  labs(x="Temperature (ºC)", y=bquote(~MO[2]* ' (mg' ~O[2] ~kg^'–1'~hour^'–1'*')'),title="E) Metabolic Metrics (5 g fish)")+
  
  ggplot(smr.preds, aes(temp, weight, fill=AS))+geom_raster()+
  theme_classic()+
  scale_fill_viridis_c(option="magma",name="AS", limits=c(0,450))+ 
  labs(x="Temperature (ºC)", y="Weight (g)",title="F) AS - Temperature & Weight") +
  
  plot_layout(ncol=2)

#paper Fig


















#model g O2 g day for bioenergetics model ----

#SMR
head(SMR.comb)
SMR.comb$SMR.g <- SMR.comb$SMR/1000/1000*24
ggplot(SMR.comb, aes(temp, SMR.g, col=weight.g))+geom_point()
SMR.comb$logSMR.g <- log(SMR.comb$SMR.g)
ggplot(SMR.comb, aes(temp, logSMR.g, col=weight.g))+geom_point()

SMR.gO2.LM <-lm(logSMR.g~ temp + logweight, data = SMR.comb)
drop1(SMR.gO2.LM, test="Chisq")
summary(SMR.gO2.LM)
exp(-6.542053)

#SMR = RA × W^RB × e^RQ×T
#SMR = 0.00137 x W ^ -0.277  x e^(0.080 x T)

RA <- 0.00144
RB <- -0.284
RQ <- 0.080



#convert to fish bioenergetics 4 terms 
fb <- merge(data.frame(temp=seq(1,33,1)), data.frame(weight=seq(1,10,1)))
ft <- exp(RQ*fb$temp)
Rmax <- RA * fb$weight ^ RB
fb$Met.J <- Rmax * ft * fb$weight * 13560
head(fb)
ggplot(fb, aes(temp, Met.J, col=weight))+geom_point()
#looks good. 

#activity scalar, fish we clearly observed to increase activity majorly with 
#temp in the lab from 6-24/26ºC
fb4 <- read.csv("data/Parameters_official.csv")
hist(fb4$BACT)

#approximate an activity scalar of 1 at 0C, 1.5 at 32C (ie 50% increase from SMR)
0.5 / 32
ACT <- 0
BACT <- 0.016

fb$act <- ACT + BACT*fb$temp
ggplot(fb, aes(temp, act))+geom_point()

fb$Met.act <- fb$Met.J+(fb$Met.J*fb$act)
ggplot(fb, aes(temp, Met.act, col=weight))+geom_point()
ggplot(fb, aes(temp, Met.J))+geom_point()+
  geom_point(data=fb, aes(temp, Met.act), col='red')





#CTmax ----

#CTmax
CTmax <- readRDS("~/github/Silver_Shiner_Bioenergetics/data/CTmax.RDS")
head(as.data.frame(CTmax))


CTmeans <- CTmax %>% group_by(Type, StartTemp) %>% summarise(mean=mean(Temp))
CTmeans


ggplot(CTmax, aes(as.factor(round(StartTemp,0)), Temp, col=Type))+
  stat_summary(fun = mean, 
               geom = "point") + 
  stat_summary(fun.data = mean_cl_normal,  
               geom = "errorbar", width=0.1) +
  geom_label(data=CTmeans, aes(as.factor(round(StartTemp,0)), mean, label=round(mean,1)), position=position_dodge(width=1))+
  geom_point(aes("6",6), pch=15, size=2, col="black")+
  geom_point(aes("11",11), pch=15, size=2, col="black")+
  geom_point(aes("18",18), pch=15, size=2, col="black")+
  geom_point(aes("20",20), pch=15, size=2, col="black")+
  xlab("Acclimation Temperature (ºC)")+ylab("Temperature (ºC)")+scale_color_manual(values=c("orange","red"))+
  scale_y_continuous(limits=c(0,35))+theme_bw()+labs(title="A) CTmax - Temperature")+
  ggplot(CTmax, aes(weight.g, Temp, col=Type))+geom_point()+geom_smooth(method="lm")+
  xlab("Weight (g)")+ylab("Temperature (ºC)")+scale_color_manual(values=c("orange","red"))+
  scale_y_continuous(limits=c(0,35))+
  plot_layout(ncol=1)+theme_bw()+labs(title="B) CTmax - Fish weight")
#paper Fig



#models
CTagg.lm <- lm(Temp ~ as.factor(StartTemp) + weight.g, data=CTmax %>% filter(Type=="Agitation"))
summary(CTagg.lm)
anova(CTagg.lm)
drop1(CTagg.lm, test="Chisq")

CTmax.lm <- lm(Temp ~ as.factor(StartTemp) + weight.g, data=CTmax %>% filter(Type=="CTmax"))
summary(CTmax.lm)
anova(CTmax.lm)
drop1(CTmax.lm, test="Chisq")
#report 











#growth and consumption -----

#lab growth
growth <- readRDS("data/growth.records.RDS")
head(growth) #recorded growth data amongst experiments 

ggplot(growth %>% filter(temp!=21), aes(as.factor(temp), growth.g))+geom_boxplot()+theme_bw()+
  labs(x="Temperature (ºC)", y="Growth (g)")

ggplot(growth %>% filter(temp!=21), aes(temp, growth.g))+geom_point()+theme_bw()+geom_smooth(method="lm", formula=y~poly(x,2))+
  labs(x="Temperature (ºC)", y="Growth (g)")

ggplot(growth %>% filter(temp!=21), aes(s.mass, growth.g, col=as.factor(temp)))+geom_point()+geom_smooth(method="lm")+theme_bw()+
  scale_color_viridis_d(option="plasma", end=0.9, name="ºC")+
  labs(x="Starting Weight (g)", y="Growth (g)")

#growth interpolated by day during experiments
growth.days <- readRDS("data/growth.by.day.RDS")
ggplot(growth.days, aes(day, weight, col=SS.ID))+geom_line()+facet_wrap(~temp, scales='free')
growth.days %>% group_by(temp) %>% summarise(min(day),max(day))


#look at growth rates in group level temperature experiment
growth.temp <- growth %>% filter(temp!=21)


growth.lm <-lm(growth.g ~ s.mass + temp + s.mass*temp, data = growth.temp %>% filter(!is.na(growth.g)))
drop1(growth.lm, test="Chisq")
summary(growth.lm)
#report, sig interaction between temperature and size. 









#consumption data
cons <- readRDS("~/github/Silver_Shiner_Bioenergetics/data/consumption.RDS")
head(cons) #feeding, waste, consumption rates by tank and day for both individual and temperature experiments 

#temperature effects plots 

#calc means per temp
cons.mean <- cons %>% group_by(temp, date) %>% summarise(cons.p=mean(cons.p), waste.g=mean(waste.g), cons.g=mean(cons.g))

#removing 21C here, which was individual feeding, plotting 26 seperate because of different timing 
ggplot(cons %>% filter(temp!=21 & temp!=26), 
       aes(date, cons.p, col=as.factor(temp), group=tank))+geom_point(alpha=0.5)+geom_line(alpha=0.5)+
  geom_point(data=cons.mean %>% filter(temp!=21 & temp!=26), 
             aes(date, cons.p, col=as.factor(temp)), inherit.aes = F, size=3)+
  geom_line(data=cons.mean %>% filter(temp!=21 & temp!=26), 
            aes(date, cons.p, col=as.factor(temp)), inherit.aes = F, size=1)+
  theme_bw()+labs(x="Date", y="Consumption (% body weight)")+
  scale_colour_viridis_d(option="plasma", end=0.8, name="Temp (ºC)")+
  scale_y_continuous(limits=c(0,0.25))+
  theme(legend.position = c(0.15, 0.8), legend.box.background = element_rect(color = "black"))+
  
  ggplot(cons %>% filter(temp==26), 
         aes(date, cons.p, col=as.factor(temp), group=tank))+geom_point(alpha=0.5)+geom_line(alpha=0.5)+
  geom_point(data=cons.mean %>% filter(temp!=21 & temp==26), 
             aes(date, cons.p, col=as.factor(temp)), inherit.aes = F, size=3)+
  geom_line(data=cons.mean %>% filter(temp!=21 & temp==26), 
            aes(date, cons.p, col=as.factor(temp)), inherit.aes = F, size=1)+
  theme_bw()+labs(x="Date", y="Consumption (% body weight)")+
  scale_colour_viridis_d(option="plasma", begin=0.9, name="Temp (ºC)")+
  scale_y_continuous(limits=c(0,0.25))+
  theme(legend.position = c(0.15, 0.8), legend.box.background = element_rect(color = "black"))



#filter to just max C periods
cons.og <- cons %>% filter(temp!=21 & temp!=26 & date>="2023-12-09") %>% group_by(temp, tank) %>%
  summarise(cons.p=mean(cons.p)) %>% as.data.frame()
cons.og

cons.26 <- cons %>% filter(temp==26 & date>="2024-02-20") %>% group_by(temp, tank) %>%
  summarise(cons.p=mean(cons.p)) %>% as.data.frame()
cons.26

cons.temp <- rbind(cons.og, cons.26)
cons.temp %>% group_by(temp) %>% summarise(mean(cons.p))

ggplot(cons.temp, aes(as.factor(temp), cons.p))+geom_boxplot()+theme_bw()+
  labs(x="Temperature (ºC)", y="Maximum Consumption (% BW)")+scale_y_continuous(limits=c(0,0.25))

ggplot(cons.temp, aes(temp, cons.p))+geom_point()+theme_bw()+geom_smooth()+
  labs(x="Temperature (ºC)", y="Maximum Consumption (% BW)")



#model temp effects 
cons.gam <- gam(cons.p ~ s(temp, k=5), data=cons.temp)
summary(cons.gam)

#predict
temps.pred <- data.frame(temp=seq(1, 30, 0.1))
temps.pred$cmax.pred <- predict(cons.gam, temps.pred, type="response")

ggplot(cons.temp, aes (temp, cons.p))+geom_point()+
  geom_line(data=temps.pred, aes(temp, cmax.pred), col="blue")+
  theme_bw()+labs(x="Temperature (ºC)", y="Maximum Consumption (% body weight)")
#peak should be reasonably reliable for CTO estimation (max predicted value)

# Find the maximum value
temps.pred %>% filter(cmax.pred == max(cmax.pred))
#23C is CTO 


#calc Q10 for 6 to 24C to get CQ
head(cons.temp)
cons.temp %>% group_by(temp) %>% summarise(mean(cons.p))
CQ <- (0.154 / 0.0310)^(10/(24-6))
CQ 
#good. 










#individual feeding experiment 
cons.ind <- cons %>% filter(temp==21)
head(cons.ind)

#plot consumption rates
ggplot(cons.ind, aes(date, cons.p, col=weight, group=tank))+geom_point()+geom_line()+scale_color_viridis_c(name="Weight (g)")+
  theme_bw()+labs(x="Date", y="Consumption (% body weight)")+
  theme(legend.position = c(0.1,0.8), legend.box.background = element_rect(color = "black"))

#averaging Feb 21, 22, was noted the 23 fish were acting shy and not feeding for some reason. 
cons.ind.max <- cons.ind %>% filter(date=="2024-02-21" | date=="2024-02-22") %>% 
  group_by(tank.temp) %>% summarise(cons=mean(cons.p), weight=mean(weight))
cons.ind.max 

ggplot(cons.ind.max, aes(log(weight), log(cons)))+geom_point()+geom_smooth(method="lm")


#convert these to highest value at 23..
temps.pred %>% filter(temp==21 | temp==23)
temps.pred$cmax.pred[2]/temps.pred$cmax.pred[1]
#not necessary 

#use this to correct values to 23 
cons.ind.max$cons.corr <- cons.ind.max$cons
ggplot(cons.ind.max, aes(weight,cons.corr))+geom_point()+geom_smooth(method='lm')

#just log weight:
cons.ind.max$log.weight <- log(cons.ind.max$weight)
ind.lm.corr <- lm(cons.corr ~ log.weight, data=cons.ind.max)
summary(ind.lm.corr)

#typical structure, log both
ind.lm.corr2 <- lm(log(cons.corr) ~ log(weight), data=cons.ind.max)
summary(ind.lm.corr2)
exp(-0.52) #these values are abnormally high for FB4.. may be an issue with
#larger fish not feeding to their max extent in the lab. 

#fitted values from log weight model for plotting
cons.preds <- data.frame(weight=seq(1,10, 0.1))
cons.preds$log.weight <- log(cons.preds$weight)

cons.preds <- cbind(cons.preds, predict(ind.lm.corr, cons.preds, interval = "confidence", level = 0.95))
head(cons.preds)

cons.preds$c.pred <- predict(ind.lm.corr, cons.preds, type="response")
head(cons.preds)

ggplot(cons.ind.max, aes(weight, cons.corr))+geom_point()+
  geom_ribbon(data=cons.preds, aes(weight, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.5, inherit.aes = F)+
  geom_line(data=cons.preds, aes(weight, fit), color = "blue")+
  theme_bw()


#pulling the intercept from the log weight model here, seems reasonable, and using the
#slope from redbelly dace from FB4, a similar species in terms of size and ecology 
#0.305 is the intercept. 
dace <- fb4 %>% filter(Species=="Dace (adult & juvenile)")
dace #-0.310 is dace slope

#plug into FB4 metrics 
CA <- 0.305
CB <- -0.310

#from temp model
CTM <- 33 #CTmax
CTO <- 23 #from gam above
CQ <- 2.436 #from Q10 above

#standard scaling equation 2 from fB4:
CY <- log(CQ) * (CTM - CTO + 2)
CZ <- log(CQ) * (CTM - CTO)
CX <- (CZ^2 * (1+(1+40/CY)^0.5)^2)/400

fb <- merge(data.frame(temp=seq(1,33,1)), data.frame(weight=seq(1,10,1)))
head(fb)

fb$V <- (CTM - fb$temp) / (CTM - CTO)
fb$ft <- ifelse(fb$temp < CTM, fb$V^CX * exp(CX * (1 - fb$V)), 0)
fb$Cmax <-  CA + log(fb$weight)*CB
fb$Cmax <- CA*fb$weight^CB
fb$Cons.p <- 1 #max consumption at 1
fb$C <- fb$Cmax * fb$Cons.p * fb$ft
fb$Cons.g <- fb$C*fb$weight
fb$Cons.J <- fb$Cons.g*2460

ggplot(fb, aes(temp, ft, col=weight))+geom_point()+
ggplot(fb, aes(temp, Cmax, col=weight))+geom_point()+
ggplot(fb, aes(temp, C, col=weight))+geom_point()+
ggplot(fb, aes(temp, Cons.J, col=weight))+geom_point()
# 






#all final consumption and feeding plots 

#calc wastes from 6 to 24 as example plot
waste.avg <- cons %>% group_by(temp, date) %>% summarise(waste.g=mean(waste.g))

#plots
ggplot(growth %>% filter(temp!=21), aes(as.factor(temp), growth.g))+geom_boxplot()+theme_bw()+
  # scale_y_continuous(limits=c(-2.5,1))+
  labs(x="Temperature (ºC)", y="Growth (g)", title="Growth", subtitle="A) Growth - Temperature")+
  
  ggplot(growth %>% filter(temp!=21), aes(s.mass, growth.g, col=as.factor(temp)))+geom_point()+geom_smooth(method="lm")+theme_bw()+
  scale_color_viridis_d(option="plasma", end=0.9, name="ºC")+
  labs(x="Starting Weight (g)", y="Growth (g)", subtitle="B) Growth - Fish Weight")+
  theme(legend.position = c(0.15,0.31), legend.box.background = element_rect(color = "black"),
        legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3,"cm"))+
  scale_y_continuous(limits=c(-3.4,1))+
  
  ggplot(cons %>% filter(temp!=21 & temp!=26), 
         aes(date, cons.p, col=as.factor(temp), group=tank))+geom_point(alpha=0.5)+geom_line(alpha=0.5)+
  geom_point(data=cons.mean %>% filter(temp!=21 & temp!=26), 
             aes(date, cons.p, col=as.factor(temp)), inherit.aes = F, size=3)+
  geom_line(data=cons.mean %>% filter(temp!=21 & temp!=26), 
            aes(date, cons.p, col=as.factor(temp)), inherit.aes = F, size=1)+
  theme_bw()+labs(x="Date", y="Consumption (g food / g fish)",title="Consumption", subtitle="C) Consumption - Temperature (6-24ºC)")+
  scale_colour_viridis_d(option="plasma", end=0.8, name="Temp (ºC)")+
  scale_y_continuous(limits=c(0,0.25))+
  theme(legend.position = c(0.15, 0.7), legend.box.background = element_rect(color = "black"),
        legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3,"cm"))+
  
  
  ggplot(cons %>% filter(temp==26), 
         aes(date, cons.p, col=as.factor(temp), group=tank))+geom_point(alpha=0.5)+geom_line(alpha=0.5)+
  geom_point(data=cons.mean %>% filter(temp!=21 & temp==26), 
             aes(date, cons.p, col=as.factor(temp)), inherit.aes = F, size=3)+
  geom_line(data=cons.mean %>% filter(temp!=21 & temp==26), 
            aes(date, cons.p, col=as.factor(temp)), inherit.aes = F, size=1)+
  theme_bw()+labs(x="Date", y="Consumption (g food / g fish)", subtitle="D) Consumption - Temperature (26ºC)")+
  scale_colour_viridis_d(option="plasma", begin=0.9, name="Temp (ºC)")+
  scale_y_continuous(limits=c(0,0.25))+
  theme(legend.position = c(0.15, 0.8), legend.box.background = element_rect(color = "black"),
        legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3,"cm"))+
  
  ggplot(cons %>% filter(temp!=21 & temp!=26), 
         aes(date, waste.g, col=as.factor(temp), group=tank))+geom_point(alpha=0.5)+geom_line(alpha=0.5)+
  geom_point(data=waste.avg%>% filter(temp!=21 & temp!=26), aes(date, waste.g, col=as.factor(temp)), inherit.aes = F, size=3)+
  geom_line(data=waste.avg%>% filter(temp!=21 & temp!=26), aes(date, waste.g, col=as.factor(temp)), inherit.aes = F, size=1)+
  theme_bw()+labs(x="Date", y="Waste Feed (g)", subtitle="E) Waste Feed - Temperature (6-24ºC)")+
  scale_colour_viridis_d(option="plasma", end=0.8, name="Temp (ºC)")+
  scale_y_continuous(limits=c(0,1.2))+
  theme(legend.position = c(0.15, 0.71), legend.box.background = element_rect(color = "black"),
        legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3,"cm"))+
  
  ggplot(cons.ind, aes(date, cons.p, col=weight, group=tank))+geom_point()+geom_line()+scale_color_viridis_c(name="Weight (g)")+
  theme_bw()+labs(x="Date", y="Consumption (g food / g fish)", subtitle="F) Consumption - Fish Weight (21ºC)")+
  theme(legend.position = c(0.15,0.65), legend.box.background = element_rect(color = "black"))+
  scale_y_continuous(limits=c(0,0.4))+theme(legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3,"cm"))+
  
  ggplot(cons.temp, aes (as.factor(temp), cons.p))+geom_boxplot()+
  theme_bw()+labs(x="Temperature (ºC)", y="Maximum Consumption (g food / g fish)", subtitle="G) Maximum Consumption - Temperature")+
  scale_y_continuous(limits=c(0,0.25))+
  
  ggplot(cons.ind.max, aes(weight, cons.corr))+geom_point()+
  geom_ribbon(data=cons.preds, aes(weight, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.5, inherit.aes = F)+
  geom_line(data=cons.preds, aes(weight, fit), color = "blue")+
  theme_bw()+
  labs(x="Weight (g)", y="Maximum Consumption (g food / g fish)", subtitle="H) Maximum Consumption - Fish Weight (21ºC)")+
  
  
  plot_layout(ncol=2)
#paper Fig

















#shiner life history data ----

#from:
#https://onlinelibrary.wiley.com/doi/full/10.1111/eff.12598?casa_token=-QknFT340tMAAAAA%3AMCIPdA5LBfRDRKqcECeq0SxJQeQgVlXiRBbXSTi7FIcHFL9a-Y_ck6oHQ0pLJTT-ITSCeDqS2ZDlH1VhkA
#L age=128.70 (1- e(−0.758(age year+0.764))

#for reference in energetics modeling 
SSage <- data.frame(age=c(0:4))
SSage$length.mm <- 128.70*(1-exp(-0.758*(SSage$age+0.764)))

#length to weight 
ss.lw <- read.csv("data/SS.LW.csv")
ggplot(ss.lw, aes(length.mm, weight.g))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2))
ggplot(ss.lw, aes(log(length.mm), log(weight.g)))+geom_point()+geom_smooth(method="lm")

ss.lw$logweight <- log(ss.lw$weight.g)
ss.lw$loglength <- log(ss.lw$length.mm)

lw.lm <- lm(logweight ~ loglength, data=ss.lw %>% filter(!is.na(weight.g)))
summary(lw.lm)

SSage$loglength <- log(SSage$length.mm)
SSage$weight.g <- exp(predict(lw.lm, SSage))
ggplot(SSage, aes(age, length.mm))+geom_path()+
  ggplot(SSage, aes(age, weight.g))+geom_path()


#so these are mean sizes at age. For now just assume that's the growth rate
SSage$Growth.g <- SSage$weight.g-lag(SSage$weight.g)
SSage$Growth.g[1] <- SSage$weight.g[1]
head(SSage)
#will use to parameterize annual energetics model 
















#bioenergetics modeling -----
summary(fb4$FA)
summary(fb4$UA) 

#set up metadata
shiner.meta <- data.frame(oxycal=13560, #fb4
                          ED=5522, #dace fb4
                          EDP=2460, #Meyer et al. 2016:
                          #https://meridian.allenpress.com/jfwm/article/7/2/388/209673/Growth-Food-Consumption-and-Energy-Status-of
                          CA=0.305, #intercept individual model
                          CB=-0.310, #slope individual model
                          CQ=2.436, #Q10 6-24 C
                          CTO=23, #from gam of cons temp
                          CTM=33, #CTmax
                          RA=0.00137, #from g2 model
                          RB=-0.277, #from g2 model
                          RQ=0.080, #from g2 model
                          ACT=0.0, #assume equals smr at 0C
                          BACT=0.016,#from scalar 0 to 50% at 33C
                          FA=0.16, #median FB4
                          UA=0.07, #median FB4
                          SDA=0.15) #dace
shiner.meta



#Consumption EQ2 
CY <- log(shiner.meta$CQ) * (shiner.meta$CTM - shiner.meta$CTO + 2)
CZ <- log(shiner.meta$CQ) * (shiner.meta$CTM - shiner.meta$CTO)
CX <- (CZ^2 * (1+(1+40/CY)^0.5)^2)/400





#basic model amongst sizes
fb.max <- merge(data.frame(temp=seq(1,33,1)), data.frame(weight=seq(1,10,1)))
head(fb.max)

fb.max$V <- (shiner.meta$CTM - fb.max$temp) / (shiner.meta$CTM - shiner.meta$CTO)
fb.max$ft <- ifelse(fb.max$temp < shiner.meta$CTM, fb.max$V^CX * exp(CX * (1 - fb.max$V)), 0)
fb.max$Cmax <- shiner.meta$CA * fb.max$weight ^ shiner.meta$CB
fb.max$Cons.p <- 1
fb.max$C <- fb.max$Cmax * fb.max$Cons.p * fb.max$ft
fb.max$Cons.g <- fb.max$C*fb.max$weight
fb.max$Cons.J <- fb.max$Cons.g*shiner.meta$EDP 

#wastes
fb.max$Eg = shiner.meta$FA*fb.max$Cons.J 
fb.max$Ex = shiner.meta$UA*(fb.max$Cons.J-fb.max$Eg)

#SDA
fb.max$SDA <- shiner.meta$SDA *(fb.max$Cons.J-fb.max$Eg)  

#shiner.metabolism
fb.max$R.ft <- exp(shiner.meta$RQ*fb.max$temp)
fb.max$Rmax <- shiner.meta$RA * fb.max$weight ^ shiner.meta$RB
fb.max$Met.J <- fb.max$Rmax * fb.max$R.ft * fb.max$weight * shiner.meta$oxycal
fb.max$Met.act <- fb.max$Met.J+(fb.max$Met.J*(shiner.meta$BACT*fb.max$temp))

#growth
fb.max$Growth.J <- fb.max$Cons.J - (fb.max$Met.act+fb.max$SDA+fb.max$Eg+fb.max$Ex)
fb.max$Growth.g <- fb.max$Growth.J/shiner.meta$ED 
fb.max$Growth.f <- as.factor(ifelse(fb.max$Growth.g<0,"negative","positive"))

head(fb.max)
fb.max %>% filter(Growth.g==max(Growth.g))
fb.max %>% filter(Cons.g==max(Cons.g))


#plots
ggplot(fb.max, aes(temp, Cons.J, col=weight))+geom_point()+scale_color_viridis_c(name="Weight (g)")+theme_bw()+
  labs(title="A) Consumption", x="Temperature (ºC)", y="Consumption (Joules)")+theme(legend.position = "none")+
  
  plot_spacer()+
  
  ggplot(fb.max, aes(temp, Met.J, col=weight))+geom_point(pch=4)+scale_color_viridis_c()+theme_bw()+theme(legend.position = "none")+
  labs(title="B) Metabolism", x="Temperature (ºC)", y="Metabolism (Joules)")+scale_y_continuous(limits=c(0,1500))+
  ggplot(fb.max, aes(temp, SDA, col=weight))+geom_point(pch=4)+scale_color_viridis_c()+theme_bw()+theme(legend.position = "none")+
  labs(x="Temperature (ºC)", y="SDA (Joules)", title="C) SDA")+scale_y_continuous(limits=c(0,1500))+
  ggplot(fb.max, aes(temp, Eg, col=weight))+geom_point(pch=4)+scale_color_viridis_c()+theme_bw()+theme(legend.position = "none")+
  labs(x="Temperature (ºC)", y="Egestion (Joules)", title="D) Egestion")+scale_y_continuous(limits=c(0,1500))+
  ggplot(fb.max, aes(temp, Ex, col=weight))+geom_point(pch=4)+scale_color_viridis_c()+theme_bw()+theme(legend.position = "none")+
  labs(x="Temperature (ºC)", y="Excretion (Joules)", title="E) Excretion")+scale_y_continuous(limits=c(0,1500))+
  ggplot(fb.max, aes(temp, Growth.J, col=weight, pch=Growth.f))+geom_point()+scale_color_viridis_c()+theme_bw()+
  scale_shape_manual(values = c(4, 19))+labs(title="Growth")+theme(legend.position = "none")+
  labs(x="Temperature (ºC)", y="Growth (Joules)", title="F) Growth (Joules)")+
  ggplot(fb.max, aes(temp, Growth.g, col=weight, pch=Growth.f))+geom_point()+scale_color_viridis_c()+theme_bw()+
  scale_shape_manual(values = c(4, 19))+theme(legend.position = "none")+
  labs(x="Temperature (ºC)", y="Growth (g)", title="G) Growth (grams)")+
  plot_layout(ncol=2)
#paper Fig






#fit across year 
source("worksheets/energetics.functions.R")

#set up grow object with daily creek temps 
creek.temps <- readRDS("data/creek.temps.RDS")
head(creek.temps)
mean(creek.temps$temp)
max(creek.temps$temp)

#plot for appendix:
ggplot(creek.temps, aes(ord.date, temp, col=temp))+geom_point()+theme_bw()+scale_color_viridis_c(option="magma", name="ºC")+
  labs(x="Ordinal Date (2017)", y="Temperature (ªC)")
#appendix Fig



#grow object
SS.grow <- creek.temps %>% group_by(ord.date) %>% summarise(temp=mean(temp)) %>%  select(day=ord.date, temp) %>% as.data.frame()
SS.grow$weight <- NA
SS.grow$weight[1] <- 4.7
SS.grow$E <- NA
SS.grow$E[1] <- SS.grow$weight[1]*as.numeric(shiner.meta$ED)
head(SS.grow)


#fit.p - determine consumption p value to acheive set growth rate
head(SSage)

#fit a 1 to 2 year old growth 
p <- fit.p(IW = 4.7, 
           FW = 7.4, 
           p = 0.5,  #starting value
           W.tol = 0.0001, #tolerance for difference from value
           max.iter = 25) #max iterations

p 



#fit this with grow:
outputs <- grow(data=SS.grow, meta=shiner.meta, ndays=nrow(SS.grow)-1,  p=0.3364)
energetics <- as.data.frame(outputs[[1]])
FW <- outputs[[2]]
head(energetics)
FW
energetics %>% filter(temp<10 & Growth.g<0) %>% summarise(max(temp)) 
energetics %>% filter(temp>10 & Growth.g<0) %>% summarise(min(temp)) 

ggplot(energetics, aes(day, Cons.J, col=temp))+geom_point()+scale_color_viridis_c(option="plasma", name="(ºC)")+theme_bw()+
  labs(x="Day of Year", y="Consumption (Joules)", title="Field Bioenergetics (Age 1 Fish)", subtitle="A) Consumption")+theme(legend.position="none")+
  ggplot(energetics, aes(day, Met.J, col=temp))+geom_point()+scale_color_viridis_c(option="plasma", name="(ºC)")+theme_bw()+
  labs(x="Day of Year", y="Metabolism (Joules)", subtitle="B) Metabolism")+
  ggplot(energetics, aes(day, Growth.g, col=temp))+geom_point()+scale_color_viridis_c(option="plasma", name="(ºC)")+theme_bw()+
  labs(x="Day of Year", y="Growth (grams)", subtitle="C) Growth")+theme(legend.position="none")+
  ggplot(energetics, aes(day, weight, col=temp))+geom_point()+scale_color_viridis_c(option="plasma", name="(ºC)")+theme_bw()+
  labs(x="Day of Year", y="Weight (grams)", subtitle="D) Weight")+theme(legend.position="none")
#paper Fig

#end energetics 













#metabolism metrics comparison ----

#make predictions for metabolism for 5g fish
met.preds <- merge(data.frame(weight=5), data.frame(temp= seq(1, 33, 0.1)))
met.preds$logweight <- log(met.preds$weight)
met.preds$smr.pred <- exp(predict(SMR.LM, met.preds))
met.preds$mmr.pred <- predict(MMR.LM, met.preds)
met.preds$AS <- met.preds$mmr.pred-met.preds$smr.pred
met.preds$AS <- ifelse(met.preds$AS<0,0, met.preds$AS)
met.preds$FAS <- met.preds$mmr.pred/met.preds$smr.pred
met.preds$AD <- met.preds$AS-met.preds$smr.pred


#grab 5gram fish bioenergetics
fb.max5 <- fb.max %>% filter(weight==5)


#feild metabolism
head(energetics)
ggplot(energetics, aes(temp, Growth.g))+geom_smooth(fill="NA")+theme_bw()


#lab growth
ggplot(growth %>% filter(temp!=21), aes(temp, growth.g))+theme_bw()+geom_smooth(method="lm", formula=y~poly(x,2), fill="NA")+
  labs(x="Temperature (ºC)", y="Growth (g)")

#model these and combine with met.preds
field.gam <- gam(Growth.g ~ s(temp), data=energetics)
met.preds$growth.field <- predict(field.gam, met.preds)
head(met.preds)
max(energetics$temp)
met.preds$growth.field <- ifelse(met.preds$temp>25, NA,as.numeric(met.preds$growth.field))
ggplot(met.preds, aes(temp, growth.field))+geom_line()

head(growth)
lab.gam <-  gam(growth.g ~ s(temp, k=4), data=growth %>% filter(temp!=21))
met.preds$growth.lab <- predict(lab.gam, met.preds)
head(met.preds)
met.preds$growth.lab <- ifelse(met.preds$temp>26 | met.preds$temp<6, NA,as.numeric(met.preds$growth.lab))
ggplot(growth %>% filter(temp!=21), aes(temp, growth.g))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2), fill="NA")+
  geom_line(data=met.preds, aes(temp, growth.lab))




#plots

ggplot(met.preds, aes(temp, AS))+geom_line(col="blue")+theme_bw()+
    labs(x="Temperature (ºC)", y="Aerobic Scope") +
ggplot(data=met.preds, aes(temp, FAS))+geom_line(col="blue")+theme_bw()+
    labs(x="Temperature (ºC)", y="Factorial Aerobic Scope") +
ggplot(data=met.preds, aes(temp, AD))+geom_line(col="blue")+theme_bw()+
     labs(x="Temperature (ºC)", y="Aerobic Differential")+
   
ggplot(data=fb.max5, aes(temp, Growth.g))+geom_line(col="#0CB702")+theme_bw()+
     labs(x="Temperature (ºC)", y="Growth Potential (g)")+
ggplot(met.preds, aes(temp, growth.field))+geom_line(col="#0CB702")+theme_bw()+
      scale_x_continuous(limits=c(0,33))+labs(x="Temperature (ºC)", y="Field growth (g)")+
ggplot(met.preds, aes(temp, growth.lab))+theme_bw()+geom_line(col="#0CB702")+
       scale_x_continuous(limits=c(0,33))+labs(x="Temperature (ºC)", y="Lab growth (g)")



#combined plot 

#standardize variables
z_score_standardize <- function(x) {
  (x - mean(x)) / sd(x)
}

met.preds$AS.z <- z_score_standardize(met.preds$AS)
met.preds$FAS.z <- z_score_standardize(met.preds$FAS)
met.preds$AD.z <- z_score_standardize(met.preds$AD)

met.preds.f <- met.preds %>% filter(!is.na(growth.field))
met.preds.f$G.field.z <- z_score_standardize(met.preds.f$growth.field)
met.preds.f %>% filter(growth.field==max(growth.field))

met.preds.f2 <- met.preds %>% filter(!is.na(growth.lab))
met.preds.f2$G.lab.z <- z_score_standardize(met.preds.f2$growth.lab)
fb.max5$Growth.z <- z_score_standardize(fb.max5$Growth.g)
met.preds.f2 %>% filter(growth.lab==max(growth.lab))

#combine into one long dataframe to set colour scheme with legend. 
met.comb <- rbind(met.preds %>% select(temp, AS.z) %>% mutate(metric="AS") %>% rename(value=AS.z),
                  met.preds %>% select(temp, FAS.z) %>% mutate(metric="FAS") %>% rename(value=FAS.z),
                  met.preds %>% select(temp, AD.z) %>% mutate(metric="AD") %>% rename(value=AD.z),
                  met.preds.f %>% select(temp, G.field.z) %>% mutate(metric="Growth.Field") %>% rename(value=G.field.z),
                  met.preds.f2 %>% select(temp, G.lab.z) %>% mutate(metric="Growth.Lab") %>% rename(value=G.lab.z),
                  fb.max5 %>% select(temp, Growth.z) %>% mutate(metric="Growth.Potential") %>% rename(value=Growth.z))
head(met.comb)
met.comb$metric <- factor(met.comb$metric, levels=c("AS","FAS","AD","Growth.Field","Growth.Lab","Growth.Potential"))

ggplot(data=met.comb, aes(temp, value, col=metric))+geom_line(size=1)+
  geom_vline(aes(xintercept=25), linetype=2, col="#E68613")+
  geom_vline(aes(xintercept=33), linetype=2, col="red")+
  annotate("label", x=25, y=1.8, label="CTag", col="#E68613")+
  annotate("label", x=33, y=1.8, label="CTmax", col="red")+
  labs(x="Temperature (ºC)", y="Z score")+theme_classic()+
  scale_colour_viridis_d(option="turbo", name="Metric")




#add combined plot to panels

upperplots <- (ggplot(met.preds, aes(temp, AS))+geom_line(col="#00A9FF",size=1)+theme_bw()+
                 labs(x=NULL, y=bquote(~'mg' ~O[2] ~kg^'–1'~hour^'–1'), title="Metabolic and Bioenergetic Metrics",
                      subtitle="A) Aerobic Scope") +
  ggplot(data=met.preds, aes(temp, FAS))+geom_line(col="#A494FF", size=1)+theme_bw()+
    labs(x=NULL, y="MMR/SMR", subtitle="B) Factorial Aerobic Scope"))/
  
  (ggplot(data=met.preds, aes(temp, AD))+geom_line(col="#FF61CC", size=1)+theme_bw()+
     labs(x=NULL, y="AS - SMR", subtitle="C) Aerobic Differential")+
     
  ggplot(data=fb.max5, aes(temp, Growth.g))+geom_line(col="#0CB702", size=1, linetype=2)+theme_bw()+
  labs(x=NULL, y=bquote("g" ~day^'–1'), subtitle="D) Growth Potential"))/
  
  (ggplot(met.preds, aes(temp, growth.field))+geom_line(col="#ABA300", size=1, linetype=2)+theme_bw()+
  scale_x_continuous(limits=c(0,33))+labs(x="Temperature (ºC)", y=bquote("g" ~day^'–1'), subtitle="E) Growth (Field)")+
    
  ggplot(met.preds, aes(temp, growth.lab))+theme_bw()+geom_line(col="#00C19A", size=1, linetype=2)+
  scale_x_continuous(limits=c(0,33))+labs(x="Temperature (ºC)", y="g",subtitle="F) Growth (Lab)"))

bottomplot <- ggplot(data=met.comb, aes(temp, value, col=metric, linetype=metric))+
  geom_vline(aes(xintercept=25), linetype=2, col="black")+
  geom_vline(aes(xintercept=33), linetype=2, col="black")+
  annotate("label", x=25, y=1.9, label="CTagitation", col="black")+
  annotate("label", x=33, y=1.9, label="CTmax", col="black")+
  geom_line(size=1)+
  scale_color_manual(values=c("#00A9FF", "#A494FF","#FF61CC","#ABA300","#00C19A","#0CB702"),name="Metric")+
  scale_linetype_manual(values=c(1,1,1,2,2,2))+guides(linetype='none')+
  labs(x="Temperature (ºC)", y="Z score", title="G) Standardized Metrics")+theme_classic()

upperplots + bottomplot + plot_layout(heights=c(1,1,1,2))
#paper Fig

