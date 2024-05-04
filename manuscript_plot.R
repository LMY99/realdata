# Make plots for manuscript
library(tidyverse)
library(splines2)
library(foreach)
rm(list=ls()); gc(verbose=FALSE)
load("biocard_result_group16nonzeros_ungrouped.RData")
# Make Point-CI plot for covariate effect ----
indice <- (R/2+1):R
summary_fixed <- apply(coefs[2:4,,indice],c(1,2),
                       function(x) c(mean=mean(x),HDInterval::hdi(x)))
dimnames(summary_fixed)[[2]] <- c("APOE4 vs Non-APOE4","Female vs Male","Education")
dimnames(summary_fixed)[[1]] <- c("mean","lower","upper")
#dimnames(summary_fixed)[[3]] <- colnames(Y)
sf <- t(Reduce(cbind, list(
  summary_fixed[,,1],
  summary_fixed[,,2],
  summary_fixed[,,3],
  summary_fixed[,,4],
  summary_fixed[,,5],
  summary_fixed[,,6],
  summary_fixed[,,7],
  summary_fixed[,,8],
  summary_fixed[,,9],
  summary_fixed[,,10],
  summary_fixed[,,11]
)))
library(ggplot2)

bionames <- c("MMSE","LogMem","DSST","ENT-THICK","HIPPO","ENT-VOL",
           "MTL","SPARE-AD","t-tau","p-tau181","AB42/AB40 Ratio")

sf <- data.frame(sf, variable=rownames(sf), biomarker=rep(bionames,each=3),
                 x=rep(c(-0.1,0,0.1),11)+rep(1:11,each=3),
                 row.names = NULL)
sink("table.txt")
sf %>% mutate(text=sprintf("% 5.2f(% 5.2f,% 5.2f)",mean,lower,upper)) %>%
  select(biomarker, text, variable) %>%
  pivot_wider(names_from = biomarker, values_from = text) %>%
  t() %>%
  janitor::row_to_names(1) %>%
  xtable::xtable() %>%
  print()
sink()
pdf("Fixed_effect.pdf",height=10)
print(
  ggplot(sf) + 
    geom_point(aes(x=x,y=mean,color=variable,shape=variable)) + 
    geom_linerange(aes(x=x,ymin=lower,ymax=upper,color=variable)) +
    scale_x_continuous(breaks=1:11,labels=bionames,minor_breaks=NULL) +
    coord_flip() + geom_hline(yintercept=0,linetype='solid',color='gray') +
    theme_bw() + theme(panel.grid.major.x=element_blank(),
                       panel.grid.minor.x=element_blank(),
                       panel.grid.major.y=element_blank()) +
    ylab('Value') + xlab('Biomarker')
)
dev.off()
# Make Goodness-of-fit Plot for select variables
fitted_values <- array(0, dim(Y))
colnames(fitted_values) <- bionames
colnames(Y) <- bionames
for(i in 1:K){
  fit <- covar.list[[i]] %*% coefs[,i,indice] + REs[df$ID,i,indice]
  fitted_values[,i] <- apply(fit,1,mean)
  fitted_values[is.na(Y[,i]),i] <- NA 
}
fitted_values <- as.data.frame(fitted_values)
Y$type <- "Data"
Y$id <- df$ID
Y$age <- df$ageori
fitted_values$type <- "Fitted"
fitted_values$id <- df$ID
fitted_values$age <- df$ageori
select_biom <- c(-1,-2,-4,-5,-7,-8,-10,-11)
gof_data <- rbind(Y, fitted_values)[,select_biom]
pdf("Goodness_of_fit_Grid3.pdf")
print(
gof_data %>% gather("Biomarker","Value",DSST,"ENT-VOL","t-tau") %>%
  ggplot() + geom_line(aes(x=age,y=Value,group=interaction(id,type),color=type),
                       alpha=0.3) +
  facet_wrap(vars(Biomarker),ncol=1,scales='free_y') +
  theme_classic() + theme(strip.background = element_blank(),
                          panel.grid.major = element_line(colour = "grey75"),
                          panel.grid.minor = element_line(colour = "grey90")) +
  scale_color_manual(values=c(Data='black',Fitted='darkred')) +
  scale_x_continuous(limits=c(50,100))
)
dev.off()
# Make CI plot for Unit-Scaled Biomarker Trajectory
agex <- seq(from=min(boundary.knot),to=max(boundary.knot),by=0.1)
basis <- ibs(agex, knots=knot.list[[1]], Boundary.knots = boundary.knot, 
            degree=2, intercept=TRUE)
basis <- basis[,3:(ncol(basis)-2)]
spline_std <- 
foreach(i=1:K,.combine=rbind) %do% {
  curves <- basis %*% coefs[-(1:4),i,indice]
  curves <- apply(curves, 2, function(x) x/max(x))
  frame <- apply(curves, 1, function(x) c(value=mean(x), HDInterval::hdi(x))) %>% t() %>%
    as.data.frame() %>% mutate(age=agex, Biomarker=bionames[i])
  frame
} %>% mutate(Biomarker=factor(Biomarker, levels=bionames), type='Curve')

spline_std <- rbind(spline_std, 
foreach(i=1:K,.combine=rbind) %do% {
  curves <- basis %*% coefs[-(1:4),i,indice]
  curves <- apply(curves, 2, function(x) x/max(x))
  curve_mean <- apply(curves, 1, mean)
  frame <- apply(curves, 2, function(x) agex[max(which(diff(x,differences=2)>=0))+1]) %>% 
    mean() %>%
    as.data.frame() %>% mutate(Biomarker=bionames[i]) %>% rename(age=1)
  frame$value <- approx(agex, curve_mean, frame$age)$y
  frame$lower <- NA; frame$upper <- NA
  frame
} %>% mutate(Biomarker=factor(Biomarker, levels=bionames), type='Mark'),
fill=TRUE
) %>% filter(!is.na(Biomarker))

colors0 <- c("darkred","red","orange","yellow","green","darkgreen","turquoise4",
             "royalblue","blue","purple","black")
names(colors0) <- bionames

inflect_order <- spline_std %>% filter(type=='Mark') %>% 
  select(age) %>% as.matrix() %>% as.vector() %>%
  order(decreasing = FALSE)
bionames <- bionames[inflect_order]
colors0 <- colors0[inflect_order]
spline_std <- mutate(spline_std, Biomarker=factor(Biomarker, levels=bionames))

p1 <- 
  ggplot() + 
  geom_line(aes(x=age,y=value,group=Biomarker,color=Biomarker), data=subset(spline_std, type=='Curve')) +
  geom_ribbon(aes(x=age,ymin=lower,ymax=upper,group=Biomarker,fill=Biomarker),alpha=0.1, 
              data=subset(spline_std, type=='Curve')) +
  geom_point(aes(x=age,y=value,color=Biomarker),data=subset(spline_std, type=='Mark'),shape='I',size=4) +
  scale_color_manual(values=colors0) +
  scale_fill_manual(values=colors0) +
  scale_x_continuous(limits=c(50,100)) +
  scale_y_continuous(limits=c(0,1)) +
  ylab("Abnormality") + theme(legend.justification = c(0,0.5),
                              axis.text.y = element_text())
pdf(width=14); print(p1);dev.off()

load('AD_Diagnose_Age.rda')
load('DECAGE.rda')
p2 <- 
rbind(cbind(age=AD_Diagnose_Age,status='Dementia'), cbind(age=decs,status='Onset')) %>%
  as.data.frame() %>%
  filter(!is.na(age)) %>%
  mutate(age=as.double(age)) %>%
  ggplot() + geom_density(aes(x=age,group=status,fill=status),bounds=c(50,100),alpha=0.3) +
  scale_x_continuous(limits=c(50,100)) +
  scale_fill_manual(values=c('Dementia'='black','Onset'='grey'),
                    labels=c('Dementia(n=28)','Symptom Onset(n=83)')) +
  theme_bw() + theme(axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     legend.justification = c(0,0.5))
pdf("SplineStd_noCI.pdf",width=14)
print(
cowplot::plot_grid(p1,p2,ncol=1,align='v',rel_heights=c(7,3))
)
dev.off()
