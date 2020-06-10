# functions for soil hydrau####
##### # calculate soil psi
psi.s.func <- function(swc,
                       psi.e = -0.74e-3,#KPa
                       b = 6.4, 
                       swc.sat = 0.38){
        psi.e * (swc / swc.sat)^-b
}
# k.soil.func(0.2,psi.e = -1)
# k.soil.func(0.2,psi.e = -.1)

# calculate soil conductivity
k.soil.func <- function(swc,
                        swc.sat = 0.13,
                        psi.e = -0.74e-3,#KPa
                        b = 6.4, 
                        k.sat = 11.4){
        
        # if(swc>swc.sat){
        #   warning(paste0('Given SWC ',swc,' is larger than satureate ',swc.sat,
        #                  '/n A max of ',k.sat,' is given'))
        #   k.s = k.sat
        # }else{
        psi.soil <- psi.s.func(swc = swc,
                               psi.e = psi.e,#KPa
                               b = b, 
                               swc.sat = swc.sat)
        
        k.s <- k.sat * (psi.e / psi.soil)^(2+3/b)
        
        k.s[which(k.s> k.sat)] <- k.sat
        return(k.s)
}

# calculate soil conductance
k.s.p.func <- function(swc,
                       r.l=0.48,
                       lai=8,
                       l.v=1490,
                       r.root=0.0005,
                       swc.sat = 0.38,
                       psi.e = -0.74e-3,#KPa
                       b = 6.4, 
                       k.sat = 11.4){
        ks <- r.l*2*pi*
                k.soil.func(swc,swc.sat,psi.e,b,k.sat)/
                lai/log10(1/sqrt(pi * l.v)/r.root)
        
        return(ks)
}

# calculate Emax 
e.frac.func <- function(swc,
                        k.plant = 0.02,
                        swc.sat = 0.38,
                        psi.e = -0.74e-3,#KPa
                        b = 6.4, 
                        k.sat = 11.4,
                        psi.min = -2,
                        r.l=0.48,
                        lai=8,
                        l.v=1490,
                        r.root=0.0005){
        
        psi.soil <- psi.s.func(swc = swc,
                               psi.e = psi.e,#KPa
                               b = b, 
                               swc.sat = swc.sat)
        
        psi.soil <- max(psi.soil , psi.min)
        
        R.p = 1 / k.plant
        R.s = 1 / k.s.p.func(swc,swc.sat=swc.sat,psi.e=psi.e,b=b,k.sat=k.sat,
                             r.l=r.l,
                             lai=lai,
                             l.v=l.v,
                             r.root=r.root)
        # -R.p * (psi.soil - psi.min)/ psi.min / (R.p + R.s)
        k.l <- 1 / (R.p + R.s)
        e.max <- k.l*(psi.soil - psi.min)
        # 
        R.s.sat <- 1 / k.s.p.func(swc.sat,swc.sat=swc.sat,psi.e=psi.e,b=b,k.sat=k.sat,
                                  r.l=r.l,
                                  lai=lai,
                                  l.v=l.v,
                                  r.root=r.root)
        
        k.l.sat <- 1/(R.p + R.s.sat)
        
        e.max.sat <- k.l.sat*(psi.e - psi.min)
        
        return(e.max/e.max.sat)
}
