rm(list=ls())
par(mfrow=c(1,1))

setwd('/home/boucherf/Documents/Flo_BACKUPS/Travail/Articles/Alooids_height/data')

########################################################################
########################################################################
######################## PART I: Load datasets #########################
########################################################################
########################################################################

# Georeferenced data from SANBI: GBIF has too many gaps outside of South Africa
# environmental variables as follows:
# leaf Area Index from ftp://neoftp.sci.gsfc.nasa.gov/geotiff/MOD15A2_M_LAI/ and NDVI from same source
# climate from WORLDCLIM sampled at 10 min (because our coordinate data are not more accurate)
QDS=read.table('QDS_binom_with_coord_and_clim.txt',h=T)
head(QDS)
dim(QDS)
sp.SANBI=as.character(sort(unique(QDS$Taxon)))
sp.SANBI

height=read.table(file='data_height_Alooideae_name_match2.txt',h=T)
library(ape)
tree=read.tree('MCC_dated_sp_level.tree')
plot(tree)
sp.phylo=as.character(sort(tree$tip.label))
sp.height=as.character(sort(unique(height$Taxon_phylo)))

sp.phylo[-which(sp.phylo%in%sp.height)] # 1 species included in phylo lacks height data
sp.height[-which(sp.height%in%sp.phylo)] # 6 species not in phylo
sp.phylo[-which(sp.phylo%in%sp.SANBI)] # 66 species in phylo with no GPS data
sp.SANBI[-which(sp.SANBI%in%sp.phylo)] # 82 spp. in SANBI records but not in phylo

######################################################
# Distribution of height across Alooids
library(moments)
hist(log10(height$Mean_vegetative_height_m),breaks=50)
agostino.test(log10(na.exclude(height$Mean_vegetative_height_m)))
hist(log10(height$Mean_flowering_height_m),breaks=50)
agostino.test(log10(na.exclude(height$Mean_flowering_height_m)))

cor(log10(height$Mean_flowering_height_m),log10(height$Mean_vegetative_height_m),use="na.or.complete",method='pearson')

########################################################################
########################################################################
############### PART II.1: DIVERSIFICATION: STRAPP #####################
########################################################################
########################################################################

# Analyse BAMM output
library(BAMMtools)
library(coda)

bamm_chain=read.table('mcmc_out.txt',sep=',',header=T)
head(bamm_chain) ; dim(bamm_chain)
apply(bamm_chain[-c(1:800),],2,effectiveSize) # 20% burnin

edata=getEventData(phy=tree,eventdata = 'event_data.txt',burnin=0.2)
summary(edata)

plot.bammdata(edata, lwd=2,legend=T)
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css)

par(mfrow=c(1,1))
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
plot.bammdata(best, lwd = 2,legend=T)
addBAMMshifts(best, cex=2.5)

# STRAPP test
edata$tip.label[-which(edata$tip.labe%in%height$Taxon_phylo)]
height$Taxon_phylo[-which(height$Taxon_phylo%in%edata$tip.label)]

LVH_bamm=rep(NA,length(edata$tip.label)) ; names(LVH_bamm)=edata$tip.label
for (i in 1:length(LVH_bamm)){
  if (names(LVH_bamm)[i]%in%height$Taxon_phylo){
    LVH_bamm[i]=log10(height$Mean_vegetative_height_m[which(height$Taxon_phylo==names(LVH_bamm)[i])]) 
  }
}

# And now test
strapp_div=traitDependentBAMM(ephy=edata,traits=LVH_bamm,reps=999,rate='net diversification')
strapp_spec=traitDependentBAMM(ephy=edata,traits=LVH_bamm,reps=999,rate='speciation')
strapp_div #r=-0.13, p=0.69
strapp_spec #r=-0.04, p=0.78


########################################################################
########################################################################
############### PART II.2: DIVERSIFICATION: MuSSE ######################
########################################################################
########################################################################

# determine clusters of height
clus3=kmeans(log10(height$Mean_vegetative_height_m),centers=3) # 3 clusters
10^(clus3$centers) # the true centers of clusters (in meters, not log-scale)
cat3=clus3$cluster ; names(cat3)=height$Taxon_phylo
# reorder the clusters by increasing height (if necessary)
cat3[which(cat3==1)]='a'
cat3[which(cat3==2)]='b'
cat3[which(cat3==3)]='c'
cat3[which(cat3=='a')]=rank(clus3$centers)[1]
cat3[which(cat3=='b')]=rank(clus3$centers)[2]
cat3[which(cat3=='c')]=rank(clus3$centers)[3]
cat3=as.numeric(cat3) ; names(cat3)=height$Taxon_phylo
table(cat3) # number of species in each cluster

if ((cat3['Haworthia_chloracantha']!=1)|(cat3['Aloidendron_dichotomum']!=3)){
  stop('!!!!! There is probably an issue with the ordering of clusters !!!!!')
}

# library(diversitree)
# sub_tree=drop.tip(tree,tip=tree$tip.label[-which(tree$tip.label%in%names(cat3))])
# # full model
# MuSSE3.0=make.musse(tree=sub_tree,states=cat3,k=3,sampling.f = 0.4)
# # We first forbid direct transition between very short and tall
# MuSSE3=constrain(MuSSE3.0,q13~0,q31~0) 
# # Model (i): BD model with different rates across height classes
# fit_full3=find.mle(func=MuSSE3,x.init=starting.point.musse(sub_tree, k=3, q.div=5, yule=FALSE)[argnames(MuSSE3)])
# # Model (ii): BD with same rates across classes
# MuSSE3.4=constrain(MuSSE3,mu2~mu1,mu3~mu1,lambda2~lambda1,lambda3~lambda1)
# fit_full3.4=find.mle(func=MuSSE3.4,x.init=starting.point.musse(sub_tree, k=3, q.div=5, yule=FALSE)[argnames(MuSSE3.4)])
# # Model (iii): PB model with different rates across height classes
# MuSSE3.1=constrain(MuSSE3,mu1~0,mu2~0,mu3~0)
# fit_full3.1=find.mle(func=MuSSE3.1,x.init=starting.point.musse(sub_tree, k=3, q.div=5, yule=FALSE)[argnames(MuSSE3.1)])
# # Model (iv):  PB model with same rate across classes
# MuSSE3.2=constrain(MuSSE3.1,lambda2~lambda1,lambda3~lambda1)
# fit_full3.2=find.mle(func=MuSSE3.2,x.init=starting.point.musse(sub_tree, k=3, q.div=5, yule=FALSE)[argnames(MuSSE3.2)])
# 
# # calculate AIC weights
# aics=rep(NA,4) ; names(aics)=c('BD','BD_1_r','PB','PB_1_lambda')
# aics[1]=2*(length(fit_full3$par)-fit_full3$lnLik)
# aics[2]=2*(length(fit_full3.4$par)-fit_full3.4$lnLik)
# aics[3]=2*(length(fit_full3.1$par)-fit_full3.1$lnLik)
# aics[4]=2*(length(fit_full3.2$par)-fit_full3.2$lnLik)
# Daics=aics-min(aics)
# aicw_MuSSE=exp(-0.5*Daics)/sum(exp(-0.5*Daics))
# aicw_MuSSE
# 
# pars=rbind(fit_full3$par.full,fit_full3.4$par.full,fit_full3.1$par.full,fit_full3.2$par.full)
# pars=cbind(aicw_MuSSE,pars)
# write.table(pars,file='All_pars_four_models_Musse.txt',row.names=F,quote=F,sep='\t')

# Now run the best model on 100 trees randomly selected from the posterior
trees=read.tree('trees100_MrBayes_dated_chronos_ingroup_only_sp_level.trees')
trees

# pars_BD_diff_rates_100trees=as.data.frame(matrix(NA,100,12))
# colnames(pars_BD_diff_rates_100trees)=names(fit_full3$par.full)
# 
# for (i in 1:length(trees)){
#   rm(fit_full3)
#   sub_tree=drop.tip(trees[[i]],tip=trees[[i]]$tip.label[-which(trees[[i]]$tip.label%in%names(cat3))])
#   # full model
#   MuSSE3.0=make.musse(tree=sub_tree,states=cat3,k=3,sampling.f = 0.4)
#   # We first forbid direct transition between very short and tall
#   MuSSE3=constrain(MuSSE3.0,q13~0,q31~0) 
#   fit_full3=find.mle(func=MuSSE3,x.init=starting.point.musse(sub_tree, k=3, q.div=5, yule=FALSE)[argnames(MuSSE3)])
#   pars_BD_diff_rates_100trees[i,]=fit_full3$par.full
#   print(i)
# }
# 
# write.table(pars_BD_diff_rates_100trees,file='All_pars_four_models_Musse.txt',row.names=F,quote=F,sep='\t')
# 
# pars_BD_diff_rates_100trees=cbind(pars_BD_diff_rates_100trees$lambda1-pars_BD_diff_rates_100trees$mu1,pars_BD_diff_rates_100trees$lambda2-pars_BD_diff_rates_100trees$mu2,pars_BD_diff_rates_100trees$lambda3-pars_BD_diff_rates_100trees$mu3,pars_BD_diff_rates_100trees)
# colnames(pars_BD_diff_rates_100trees)[1:3]=c("NetDiv1","NetDiv2","NetDiv3")
# write.table(pars_BD_diff_rates_100trees,file='All_pars_four_models_Musse_100trees_posterior.txt',row.names=F,quote=F,sep='\t')
# 
# apply(pars_BD_diff_rates_100trees,2,mean)
# apply(pars_BD_diff_rates_100trees,2,sd)


########################################################################
########################################################################
################## PART III: EVOLUTION OF PLANT HEIGHT #################
########################################################################
########################################################################

########################################################################
##################### PART III.1: Maximum-likelihood ###################
########################################################################

###############################
# load functions needed to fit the FPK model
library(BBMV)

height_binom=cbind(height[,1],sapply(as.character(height$Taxon_phylo),function(x){paste(strsplit(x,split='_')[[1]][1],'_',strsplit(x,split='_')[[1]][2],sep='')}),height[,c(3:8)])
colnames(height_binom)=colnames(height)

# Create a table with log10 height measurement per species, matching names in the phylogeny
sp_list=sort(tree$tip.label)
sp_intersect2=sp_list[which(sp_list%in%height_binom$Taxon_phylo)]
tab2=cbind(sp_intersect2,rep(NA,length(sp_intersect2)),rep(NA,length(sp_intersect2))) ; colnames(tab2)=c('Species','Veg_height','Flower_height')
row.names(tab2)=NULL
for (i in 1:dim(tab2)[1]){
  sp_match=which(height_binom$Taxon_phylo==sp_intersect2[i])
  if (length(sp_match)==1){
    tab2[i,2:3]=as.numeric(height_binom[sp_match,c(4,5)])
  }
  else{
    veg=mean(as.numeric(height_binom[sp_match,c(4)]),na.rm=T)
    flo=mean(as.numeric(height_binom[sp_match,c(5)]),na.rm=T)
    tab2[i,2:3]=c(veg,flo)
  }
}
tab2=as.data.frame(tab2)
for (i in 2:3){tab2[,i]=as.numeric(as.character(tab2[,i]))}
head(tab2)
LVH=log10(tab2$Veg_height) ; names(LVH)=tab2$Species
LVH=na.exclude(LVH)
sub_tree=drop.tip(tree,tip=tree$tip.label[-which(tree$tip.label%in%names(LVH))])

###########################################
# Fit different evolutionary models using Maximum-likelihood
par(mfrow=c(1,1))
# Complete FPK model
ll_FPK4=lnL_FPK(tree=sub_tree,trait=LVH,Npts=100,a=NULL,b=NULL,c=NULL) # the full model
fit4=find.mle_FPK(model=ll_FPK4,init.optim=c(-5,5,-5,0))
get.landscape.FPK(fit=fit4)

# BBMV model: Create four different likelihood functions, veg height: min=1cm max=50m
ll_BBMV4=lnL_BBMV(sub_tree,trait=LVH,Npts=100,bounds=c(min(LVH),log10(50)),a=NULL,b=NULL,c=NULL)
ll_BBMV2=lnL_BBMV(sub_tree,trait=LVH,Npts=100,bounds=c(min(LVH),log10(50)),a=0,b=NULL,c=NULL)
ll_BBMV1=lnL_BBMV(sub_tree,trait=LVH,Npts=100,bounds=c(min(LVH),log10(50)),a=0,b=0,c=NULL)
ll_BBMV0=lnL_BBMV(sub_tree,trait=LVH,Npts=100,bounds=c(min(LVH),log10(50)),a=0,b=0,c=0) # this is the BBM model

# fit the four BBMV models and plot the macroevolutionary landscapes
par(mfrow=c(1,1))
fit4b=find.mle_FPK(model=ll_BBMV4)
get.landscape.FPK(fit=fit4b)
fit2b=find.mle_FPK(model=ll_BBMV2)
get.landscape.FPK(fit=fit2b)
fit1b=find.mle_FPK(model=ll_BBMV1)
get.landscape.FPK(fit=fit1b)
fit0b=find.mle_FPK(model=ll_BBMV0)
get.landscape.FPK(fit=fit0b)

# fit BM and OU models
library(geiger)
BM=fitContinuous(sub_tree,LVH,model='BM')
OU=fitContinuous(sub_tree,LVH,model='OU')

# compare the fits using AIC
fit4$aic 
OU$opt$aic 
BM$opt$aic 
fit4b$aic 
fit2b$aic 
fit1b$aic 
fit0b$aic 

# AIC weights
aics=c(fit4$aic,OU$opt$aic,BM$opt$aic,fit4b$aic,fit2b$aic,fit1b$aic,fit0b$aic)
Daics=aics-min(aics)
aicw_evol=exp(-0.5*Daics)/sum(exp(-0.5*Daics))
names(aicw_evol)=c('FPK','OU','BM','BBMV4','BBMV2','BBMV1','BBM')
aicw_evol

write.table(aicw_evol,file='FPK_BBMV_results_Sept2018/AICs_evol_models_MCC_tree.txt',row.names=F,quote=F,sep='\t')
save(fit4b,file='FPK_BBMV_results_Sept2018/fit4b_MCC_tree.Rdata')

# measure characteristic time of the best model and compare it to tree depth
charac_time(fit=fit4b)
max(branching.times(sub_tree))

# Uncertainty and distribution of trait at root
fit4b$par
Uncert4=Uncertainty_FPK(fit=fit4b,tree=sub_tree,trait=LVH,Npts=100,effort_uncertainty= 50,scope_a=c(-2,2),scope_b=c(-2,2),scope_c=c(-5,5))
Uncert4
10^(fit4b$root[which(fit4b$root[,2]==max(fit4b$root[,2])),1]) # ML root value
10^(Uncert4$CI95_root) # 95% CI
save(Uncert4,file='FPK_BBMV_results_Sept2018/Uncert4b_MCC_tree.Rdata')


########################################################################
######################## PART III.2: MCMC runs #########################
########################################################################

chain1=MH_MCMC_FPK(tree=sub_tree,trait=LVH,bounds=c(min(LVH),log10(50)),Nsteps=200000,record_every=100,plot_every=100,Npts=50,pars_init=c(1,4,0,5,10),prob_update=c(0.05,0.1,0.12,0.3,0.43),verbose=TRUE,plot=TRUE,save_to='FPK_BBMV_results_Sept2018/MCMC_LVH_BBMV4_50pts.Rdata',save_every=100,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.5,0.5,0.5,0.5,1),prior.only=F)

chain2=MH_MCMC_FPK(tree=sub_tree,trait=LVH,bounds=c(min(LVH),log10(50)),Nsteps=200000,record_every=100,plot_every=100,Npts=50,pars_init=c(-2,-5,1,5,1),prob_update=c(0.05,0.1,0.12,0.3,0.43),verbose=TRUE,plot=TRUE,save_to='FPK_BBMV_results_Sept2018/MCMC_LVH_BBMV4_50pts_run2.Rdata',save_every=100,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.5,0.5,0.5,0.5,1),prior.only=F)

   # load the two chains
load('FPK_BBMV_results_Sept2018/MCMC_LVH_BBMV4_50pts.Rdata') ; chain1=chain
load('FPK_BBMV_results_Sept2018/MCMC_LVH_BBMV4_50pts_run2.Rdata') ; chain2=chain

# Measure effective sample size and convergence of the two chains
library(coda)
chain1=as.data.frame(chain1)
chain2=as.data.frame(chain2)
apply(chain1[-c(1:floor(dim(chain1)[1]/10)),],2,effectiveSize)
apply(chain2[-c(1:floor(dim(chain2)[1]/10)),],2,effectiveSize)

c1=as.mcmc(chain1,start=floor(dim(chain1)[1]/10),end=dim(chain1)[1])
c2=as.mcmc(chain2,start=floor(dim(chain2)[1]/10),end=dim(chain2)[1])
gelman.plot(x=list(c1[,c(2:8)],c2[,c(2:8)]))
gelman.diag(x=list(c1[,c(2:8)],c2[,c(2:8)]))

# Plot the traces
par(mfrow=c(3,3)) 
for (i in 2:9){
  plot(chain1[-c(1:floor(dim(chain1)[1]/10)),i],type='l')
  lines(chain2[-c(1:floor(dim(chain2)[1]/10)),i],type='l',col=2)
}

# Now plot landscape estimated: this will produce Fig. 1 of the article
# ML estimate
load('fit4b_MCC_tree.Rdata')
# combined MCMC chains
chain_comb=rbind(chain1[-c(1:floor(dim(chain1)[1]/10)),],chain2[-c(1:floor(dim(chain2)[1]/10)),])
dim(chain_comb)
#get.landscape.FPK.MCMC(chain=chain_comb,bounds=fit4b$par_fixed$bounds,Npts=100,burnin=0.0,probs.CI=c(0.05,0.95),COLOR_MEDIAN='darkred',COLOR_FILL='darkred',transparency=0.3,main='',ylab='Probability density',xlab='Plant height (log-10 scale)',xlim=NULL,ylim=c(0,1.5))
get.landscape.FPK.MCMC(chain=chain_comb,bounds=fit4b$par_fixed$bounds,Npts=100,burnin=0.0,probs.CI=c(0.05,0.95),COLOR_MEDIAN='BLACK',COLOR_FILL='grey55',transparency=0.3,main='',ylab='Probability density',xlab='Plant height (log-10 scale)',xlim=NULL,ylim=c(0,1.25))

hist(LVH,add=T,breaks=25,freq=F,col=adjustcolor(col='grey10',alpha.f=0.2))

chain=chain_comb;bounds=fit4b$par_fixed$bounds;Npts=100;burnin=0.0;probs.CI=c(0.25,0.75);COLOR_MEDIAN='BLACK';COLOR_FILL='grey5';transparency=0.3
chain2=chain[-c(1:floor(dim(chain)[1]*burnin)),]  # remove burnin
all_V=as.data.frame(matrix(NA,dim(chain2)[1],Npts))
step=(bounds[2]-bounds[1])/(Npts-1)
SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
for (i in 1:dim(chain2)[1]){
  temp=chain2[i,'a']*SEQ^4+chain2[i,'b']*SEQ^2+chain2[i,'c']*SEQ #potential
  all_V[i,]=exp(-temp)/sum(exp(-temp)*step)
}
lines(apply(all_V,2,function(x){quantile(x,probs=c(0.5))})~seq(from=bounds[1],to=bounds[2],length.out=Npts),col=COLOR_MEDIAN,lwd=3,type='l')
polygon(x=c(seq(from=bounds[1],to=bounds[2],length.out=100),seq(from=bounds[2],to=bounds[1],length.out=100)),y=c(apply(all_V,2,function(x){quantile(x,probs=probs.CI[1])}),rev(apply(all_V,2,function(x){quantile(x,probs=probs.CI[2])}))),col=adjustcolor(col=COLOR_FILL,alpha.f=transparency))

V=fit4b$par$a*SEQ^4+fit4b$par$b*SEQ^2+fit4b$par$c*SEQ #potential
lines((exp(-V)/sum(exp(-V)*step))~seq(from=bounds[1],to=bounds[2],length.out=Npts),lty='dashed',col=1,lwd=3)

#legend(x='topright',legend=c('Maximum-Likelihood estimate','Median of the posterior','50% credible interval','95% credible interval','Observed distribution'),fill=c(1,2,adjustcolor(col='red',alpha.f=transparency),adjustcolor(col='darkred',alpha.f=transparency),col=adjustcolor(col='blue',alpha.f=0.2)))


########################################################################
################## PART III.3: ML fit on 100 trees #####################
########################################################################

fit4b_ML_100trees=as.data.frame(matrix(NA,100,5)) ; colnames(fit4b_ML_100trees)=c("sigsq","a","b","c","root")
root_100_trees=as.data.frame(matrix(NA,100,101)) ; root_100_trees[,1]=fit4b$root$root
for (i in 1:100){
  sub_tree=drop.tip(trees[[i]],tip=trees[[i]]$tip.label[-which(trees[[i]]$tip.label%in%names(LVH))])
  ll_BBMV4=lnL_BBMV(sub_tree,trait=LVH,Npts=100,bounds=c(min(LVH),log10(50)),a=NULL,b=NULL,c=NULL)
  fit4b_temp=find.mle_FPK(model=ll_BBMV4)
  fit4b_ML_100trees[i,]=c(fit4b_temp$par,10^(fit4b_temp$root[which(fit4b_temp$root[,2]==max(fit4b_temp$root[,2])),1]))
  root_100_trees[,(i+1)]=fit4b_temp$root$density
  print(i)
}

save(root_100_trees,file='root_ACE_100_trees.Rdata')
overall_dens=apply(root_100_trees[,-c(1)],1,mean)
plot(overall_dens~as.numeric(10^(root_100_trees[,1])),type='l',ylim=c(0,0.03))
cum_dens=rep(NA,length(overall_dens))
for (i in length(cum_dens):1){cum_dens[i]=sum(overall_dens[c(i:length(overall_dens))])}
cum_dens
10^root_100_trees[which(cum_dens<0.95)[1],1] # 4 cm
10^root_100_trees[which(cum_dens<0.9)[1],1] # 13 cm
10^root_100_trees[which(cum_dens<0.5)[1],1] # 5.8 m

  
Tc=rep(NA,100)
SEQ = seq(from = -1.5, to = 1.5, length.out = 100)
bounds=c(min(LVH),log10(50))
for (i in 1:100){
  Vq = fit4b_ML_100trees[i,'a'] * SEQ^4 + fit4b_ML_100trees[i,'b'] * SEQ^2 + fit4b_ML_100trees[i,'c'] * SEQ
  Mat = DiffMat_forward(Vq)
  vp = diag(Mat$diag)
  Tc[i] = (2 * (bounds[2] - bounds[1])^2/fit4b_ML_100trees[i,'sigsq'])/(100-1)^2/abs(sort(Re(vp), decreasing = T)[2])
}
fit4b_ML_100trees=cbind(fit4b_ML_100trees,Tc)

save(fit4b_ML_100trees,file='fit4b_100_posterior_tree.Rdata')
apply(fit4b_ML_100trees,2,summary)
boxplot(fit4b_ML_100trees[,c('a','b','c')])

all_V_ML=matrix(NA,100,100)
Npts=100
step=(fit4b$par_fixed$bounds[2]-fit4b$par_fixed$bounds[1])/(Npts-1)
SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  for (i in 1:dim(fit4b_ML_100trees)[1]){
    temp=fit4b_ML_100trees[i,'a']*SEQ^4+fit4b_ML_100trees[i,'b']*SEQ^2+fit4b_ML_100trees[i,'c']*SEQ #potential
    all_V_ML[i,]=exp(-temp)/sum(exp(-temp)*step)
  }
xlim=fit4b$par_fixed$bounds
ylim=c(0,1.2)
par(mfrow=c(1,1),mar=c(5,4,1,1))
bounds=xlim
plot(apply(all_V_ML,2,function(x){quantile(x,probs=c(0.5))})~seq(from=bounds[1],to=bounds[2],length.out=Npts),col='red',lwd=3,type='l',main=NULL,ylab='N.exp(-V)',xlab='Plant height (log-10 scale)',xlim=xlim,ylim=ylim)
hist(LVH,add=T,breaks=25,freq=F,col=adjustcolor(col='grey10',alpha.f=0.2))
probs.CI=c(0.1,0.9)
polygon(x=c(seq(from=bounds[1],to=bounds[2],length.out=Npts),seq(from=bounds[2],to=bounds[1],length.out=Npts)),y=c(apply(all_V_ML,2,function(x){quantile(x,probs=probs.CI[1])}),rev(apply(all_V_ML,2,function(x){quantile(x,probs=probs.CI[2])}))),col=adjustcolor(col='darkred',alpha.f=0.3))
probs.CI=c(0.25,0.75)
polygon(x=c(seq(from=bounds[1],to=bounds[2],length.out=Npts),seq(from=bounds[2],to=bounds[1],length.out=Npts)),y=c(apply(all_V_ML,2,function(x){quantile(x,probs=probs.CI[1])}),rev(apply(all_V_ML,2,function(x){quantile(x,probs=probs.CI[2])}))),col=adjustcolor(col='red',alpha.f=0.3))

add.ML.landscape.FPK(fit=fit4b,Npts=100,COLOR=1,LTY='dashed')

########################################################################
########################################################################
################ PART IV: CONSEQUENCES OF PLANT HEIGHT #################
########################################################################
########################################################################

########################################################################
################ PART IV.1: Diversification in space ###################
########################################################################

######################################################
# Genus level analysis
genus=sapply(as.character(QDS$Taxon),function(x){strsplit(x,split='_')[[1]][1]})
table(genus) # 11 genera
Genus_QDS=unique(cbind(genus,as.character(QDS$QDS)))
Genus_QDS=na.exclude(Genus_QDS)
head(Genus_QDS)
Area_genus=table(Genus_QDS[,1])
Area_genus=Area_genus*625 # each QDS is c. 25x25 km
Rich_genus=c(400,7,6,1,6,23,3,42,18,2,4) # From Manning 2014
names(Rich_genus)=names(Area_genus)
plot(Rich_genus~as.numeric(Area_genus),log='xy',pch=19)
Rich_genus['Aristaloe']=NA # if we want to remove Aristaloe (only one species)

names(Rich_genus)
crownAge_genus=c(20.8,14.6,19.5,NA,6.8,11.5,NA,12.8,12.5,NA,6.8) # median from MCC tree, NA for singletons
stemAge_genus=c(22.5,19.3,25.8,NA,9.0,14.2,5.5,19.3,15.6,NA,9.0) # median from MCC tree
divMOM=log(Rich_genus)/stemAge_genus

median_height_genus=rep(NA,length(Rich_genus)) ; names(median_height_genus)=names(Rich_genus)
for (i in 1:length(median_height_genus)){
  median_height_genus[i]=median(height$Mean_vegetative_height_m[which(height$Genus==names(median_height_genus)[i])],na.rm=T)
}
median_height_genus

# How does median height in a genus influence diversification rates per time and per unit area?
# prepare genus skeleton tree
genus_tree=drop.tip(tree,tip=tree$tip.label[-which(tree$tip.label%in%c('Aloe_ferox','Aloidendron_dichotomum','Tulista_minima','Haworthiopsis_coarctata','Haworthia_pubescens','Kumara_plicatilis','Gonialoe_variegata','Gasteria_glauca','Astroloba_corrugata','Aristaloe_aristata','Aloiampelos_commixta'))])
genus_tree$tip.label=sapply(genus_tree$tip.label,FUN=function(x){return(strsplit(x,split='_')[[1]][1])})
plot(genus_tree)

##################################################################@
# Relationship between diversification rate vs. species packing at the genus level
genus2=sapply(unique(as.character(QDS$Taxon)),function(x){strsplit(x,split='_')[[1]][1]})
table(genus2) # 11 genera
Rich_genus_SA=table(genus2)[names(Rich_genus)]
Rich_genus_SA['Aristaloe']=1
Rich_genus_SA['Gonialoe']=3
dat_genus2=as.data.frame(cbind(names(median_height_genus),median_height_genus,Area_genus,Rich_genus,stemAge_genus,crownAge_genus,Rich_genus_SA))
colnames(dat_genus2)[1]='genus' ; row.names(dat_genus2)=NULL
for (i in 2:7){dat_genus2[,i]=as.numeric(as.character(dat_genus2[,i]))}

library(caper)
compdatgenus2=comparative.data(phy=genus_tree,data=dat_genus2,names.col = 'genus')
compdatgenus3=comparative.data(phy=genus_tree,data=dat_genus2[,c(1,2,3,7)],names.col = 'genus')
summary(pgls(I(log(Rich_genus_SA)/log(Area_genus))~I(log(median_height_genus)),data=compdatgenus3,lambda='ML'))


# plot genus level
names(stemAge_genus)=names(Rich_genus)
divdens=I(log(Rich_genus_SA)/log(Area_genus)) #; divdens['Aristaloe']='NA'
par(mar=c(4,4,2,2))

plot(as.numeric(divdens)~as.numeric(log(median_height_genus)),pch=19,xlim=c(-3.5,2.5),ylim=c(0.,0.4),xlab='Median height genus (log-scale)',ylab='Species packing - log(n)/log(Area)')
text(y=as.numeric(divdens)+0.012,x=as.numeric(log(median_height_genus)),label=names(median_height_genus),col='grey50',font=)


######################################################
# Analysis on non-named clades from the MCC and on 100 trees from the posterior

library(phytools)


div_Area_varying_t=as.data.frame(matrix(NA,11,1)) ; colnames(div_Area_varying_t)=c('time')
div_Area_varying_t$time=c(12,11,10,9,8,7,6,5,4,3,2)

cut_and_fit_div_area_and_rate_100trees=function(timeslice,tree){
  slices=treeSlice(tree,slice=(max(branching.times(tree))-timeslice))
  # create skeleton tree for pgls
  skel_list=c()
  for (i in 1:length(slices)){skel_list=c(skel_list,slices[[i]]$tip.label[1])}
  skeleton=drop.tip(tree,tip=tree$tip.label[-which(tree$tip.label%in%skel_list)])
  # need to rename these tips
  clades=rep(NA,length(slices))
  for (i in 1:length(slices)){clades[i]=paste('Clade_',i,sep='')}
  for (tip in 1:length(skeleton$tip.label)){
    num=which(skel_list==skeleton$tip.label[tip])
    skeleton$tip.label[tip]=clades[num]
  }
  Rich_subtree=Crown_age_subtree=Area_subtree=Median_height_subtree=rep(NA,length(slices))
  for (i in 1:length(slices)){
    subtree=slices[[i]]
    tips=subtree$tip.label
    Rich_subtree[i]=length(tips)
    Crown_age_subtree[i]=max(branching.times(subtree))
    Median_height_subtree[i]=median(height_binom$Mean_vegetative_height_m[which(height_binom$Taxon_phylo%in%tips)])
    Area_subtree[i]=length(unique(QDS$QDS[which(QDS$Taxon%in%tips)]))*625
  }
  Area_subtree[which(Area_subtree==0)]=NA
  Div_Area=log(Rich_subtree)/log(Area_subtree)
  Div_Rate=log(Rich_subtree)/timeslice
  data_sub=as.data.frame(cbind(clades,Rich_subtree,Crown_age_subtree,Area_subtree,Median_height_subtree,Div_Area,Div_Rate))
  for (i in 2:7){data_sub[,i]=as.numeric(as.character(data_sub[,i]))}
  comp.dat=comparative.data(phy=skeleton,data=data_sub,names.col=clades)
  p=pgls(Div_Area~log10(Median_height_subtree),data=comp.dat,lambda='ML')
  r2=summary(p)$adj.r.squared
  coef=summary(p)$coefficients[2,1]
  pval=summary(p)$coefficients[2,4]
  p2=pgls(Div_Rate~log10(Median_height_subtree),data=comp.dat,lambda='ML')
  r22=summary(p2)$adj.r.squared
  coef2=summary(p2)$coefficients[2,1]
  pval2=summary(p2)$coefficients[2,4]
  res_divdens=c(timeslice,coef,pval,r2,length(slices))
  res_divrate=c(timeslice,coef2,pval2,r22,length(slices))
  return(list(Div_Area=Div_Area,Div_Rate=Div_Rate,res_divdens=res_divdens,res_divrate=res_divrate,COR=cor(na.exclude(cbind(Div_Area,Div_Rate)))[1,2]))
}
cut_and_fit_div_area_and_rate_100trees(10,tree)

trees=read.tree('trees100_MrBayes_dated_chronos_ingroup_only_sp_level.trees')

# Get results 100 trees for how plant height influences diversity density
NCLADES=as.data.frame(matrix(NA,100,11))
SLOPE=PVAL=R2=NCLADES

for (i in 1:100){
  for(j in 1:length(div_Area_varying_t$time)){
    temp=try(cut_and_fit_div_area_and_rate_100trees(timeslice=div_Area_varying_t$time[j],tree=trees[[i]]))
    if (class(temp)!='try-error'){
      SLOPE[i,j]=temp$res_divdens[2] ; PVAL[i,j]=temp$res_divdens[3]; R2[i,j]=temp$res_divdens[4] ; NCLADES[i,j]=temp$res_divdens[5] 
    }
  }
}

par(mfrow=c(2,2))
boxplot(SLOPE,col='grey55',names=as.character(I(-div_Area_varying_t$time)),main='SLOPE divdens/H')
abline(h=0,col=2)
boxplot(PVAL,col='grey55',names=as.character(I(-div_Area_varying_t$time)),main='PVAL divdens/H')
abline(h=0.05,col=2)
boxplot(R2,col='grey55',names=as.character(I(-div_Area_varying_t$time)),main='R2 divdens/H')
abline(h=0,col=2)
boxplot(NCLADES,col='grey55',names=as.character(I(-div_Area_varying_t$time)),main='NCLADES divdens/H')

R2_sp_packing=R2
apply(R2_sp_packing,2,summary)

PVAL_sp_packing=PVAL
apply(PVAL_sp_packing,2,summary)

slopes=cbind(SLOPE,rep('Diversity_density',100))

# Get results 100 trees for how plant height influences diversification rate
NCLADES=as.data.frame(matrix(NA,100,11))
COR=SLOPE=PVAL=R2=NCLADES

for (i in 1:100){
  for(j in 1:length(div_Area_varying_t$time)){
    temp=try(cut_and_fit_div_area_and_rate_100trees(timeslice=div_Area_varying_t$time[j],tree=trees[[i]]))
    if (class(temp)!='try-error'){
      SLOPE[i,j]=temp$res_divrate[2] ; PVAL[i,j]=temp$res_divrate[3]; R2[i,j]=temp$res_divrate[4] ; NCLADES[i,j]=temp$res_divrate[5] ; COR[i,j]=temp$COR
    }
  }
}

par(mfrow=c(2,2))
boxplot(SLOPE,col='grey55',names=as.character(I(-div_Area_varying_t$time)),main='SLOPE divrate/H')
boxplot(PVAL,col='grey55',names=as.character(I(-div_Area_varying_t$time)),main='PVAL divrate/H')
abline(h=0.05,col=2)
boxplot(R2,col='grey55',names=as.character(I(-div_Area_varying_t$time)),main='R2 divrate/H')
boxplot(NCLADES,col='grey55',names=as.character(I(-div_Area_varying_t$time)),main='NCLADES divrate/H')

slopes=as.data.frame(slopes)
All_slopes_df=as.data.frame(matrix(NA,dim(SLOPE)[1]*dim(SLOPE)[2],3))
colnames(All_slopes_df)=c('slope','time','measure')
All_slopes_df$measure=rep('Diversity_density',dim(All_slopes_df)[1])
All_slopes_df$time=rep(-div_Area_varying_t$time,100)
All_slopes_df$slope=as.vector(t(slopes[,c(1:11)]))

slopes2=cbind(SLOPE,rep('Diversification_rate',100))
All_slopes_df2=as.data.frame(matrix(NA,dim(SLOPE)[1]*dim(SLOPE)[2],3))
colnames(All_slopes_df2)=c('slope','time','measure')
All_slopes_df2$measure=rep('Diversification_rate',dim(All_slopes_df2)[1])
All_slopes_df2$time=rep(-div_Area_varying_t$time,100)
All_slopes_df2$slope=as.vector(t(slopes2[,c(1:11)]))
All_slopes_df_tot=rbind(All_slopes_df,All_slopes_df2)
All_slopes_df_tot$slope=as.numeric(as.character(All_slopes_df_tot$slope))
All_slopes_df_tot$time=factor(x=as.character(All_slopes_df_tot$time),levels=c(-div_Area_varying_t$time))
All_slopes_df_tot$measure=as.factor(as.character(All_slopes_df_tot$measure))


par( mfrow=c(1,2))
boxplot(All_slopes_df_tot$slope[which(All_slopes_df_tot$measure=='Diversification_rate')]~All_slopes_df_tot$time[which(All_slopes_df_tot$measure=='Diversification_rate')],notch=T,col='grey55',xlab='Time slice (Ma)',ylab='Slope of the relationship: diversification rate vs. plant height')
abline(h=0,col=2,lty='dashed',lwd=2)
text(labels=c('A.'),x=0.6,y=0.25,cex=1.5)
boxplot(All_slopes_df_tot$slope[which(All_slopes_df_tot$measure=='Diversification_rate')]~All_slopes_df_tot$time[which(All_slopes_df_tot$measure=='Diversification_rate')],notch=T,col='grey55',add=T)

boxplot(All_slopes_df_tot$slope[which(All_slopes_df_tot$measure=='Diversity_density')]~All_slopes_df_tot$time[which(All_slopes_df_tot$measure=='Diversity_density')],notch=T,col='grey55',xlab='Time slice (Ma)',ylab='Slope of the relationship: species packing vs. plant height')
abline(h=0,col=2,lty='dashed',lwd=2)
text(labels=c('B.'),x=0.6,y=0.038,cex=1.5)


# CORRELATION BETWEEN DIVERSIFICATION RATE AND SPECIES PACKING AMONG UNNAMED CLADES FROM THE PHYLOGENETIC POSTERIOR
#apply(COR,2,summary)

