source('r/function_soil_hydro.R')
# 
# calculate max E reduction ####
soil.par.df = data.frame(
        type = c('sand','loamy sand',
                 'loam','clay loam','light clay'),
        b = c(2.79,4.26,5.25,8.17,11.55),
        psi.e = -c(0.68,0.36,3.48,2.58,4.58),
        k.sat = c(264.3,79.8,19.1,13.9,5.5)
)
swc.vec = seq(0.01,0.38,by=0.005)

e.frac.ls = list()

for (i in seq_along(soil.par.df$type)) {
        
        temp.df = soil.par.df[soil.par.df$type == soil.par.df$type[i],]
        
        
        e.frac.ls[[i]]= e.frac.func(swc.vec,
                                    psi.e = temp.df$psi.e * 10^-3,
                                    b = temp.df$b,
                                    k.sat = temp.df$k.sat,
                                    swc.sat = 0.38,
                                    lai = 0.5)
}
library(mgcv)
# get john's ros data####
ros.pd.df=read.csv('ROS_MD_PM_LWP-GX-TDR_20121031-20130516_L1.csv')
# normalise to each species
get.norm.func = function(df.in){
        df.in = df.in[order(df.in$TDR),]
        vec.in = df.in$Cond
        quant.5.95 = quantile(vec.in,c(.05,.95),na.rm=T)
        df.in$gs.norm = (vec.in - quant.5.95[[1]]) / ( quant.5.95[[2]] -  quant.5.95[[1]])
        fit.gam = gam(gs.norm~s(TDR,k=6),data = df.in)
        df.in$gs.norm.gam = fit.gam$fitted.values
        # plot(gs.norm.gam~TDR,data = df.in)
        return(df.in)
}

tmp.ls = split(ros.pd.df,ros.pd.df$Species)

tmp.ls= lapply(tmp.ls,get.norm.func)

ros.pd.df.norm = do.call(rbind,tmp.ls)

# we could also fit for g1 but there's dry and wet treatment
# library(plantecophys)
# 
# ros.pd.df.norm$group = paste0(ros.pd.df.norm$gxDate,'_',ros.pd.df.norm$Species)
# 
# g1.vec=fitBBs(ros.pd.df.norm,'group',gsmodel='BBOpti',
#       varnames = list(ALEAF = "Photo", GS = "Cond", VPD = "VpdL",
#                       Ca = "CO2S"))
# g1.fit.df = coef(g1.vec)
# 
# ros.pd.df.norm.g1 = merge(ros.pd.df.norm,g1.fit.df[,c('group','g1')],all.x =T)

pdf('reduction of gs over swc.pdf',width = 8,height = 8*0.618)
# makes plots####
with(ros.pd.df.norm,plot(gs.norm~TDR,col=Species,pch=16))
# plot gam fits to each species of the raw data
for (i in  seq_along(unique(ros.pd.df.norm$Species))) {
        points(gs.norm.gam~TDR,
               data = ros.pd.df.norm[ros.pd.df.norm$Species == unique(ros.pd.df.norm$Species)[i],],
               type='l',lty='dashed',col=i)
}

# plot the max E
for (i in  seq_along(soil.par.df$type)) {
        points(e.frac.ls[[i]]~c(swc.vec * 100),
               type='l',col = 'grey')
}
dev.off()

