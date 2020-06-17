# functions for soil hydrau####
##### # calculate soil psi
psi.s.func <- function(swc,
                       psi.e = -0.74e-3,#KPa
                       b = 6.4, 
                       swc.sat = 0.38){
        psi.e * (swc / swc.sat)^-b
}

psi.s.func(0.1)
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
e.max.func <- function(swc,
                       k.plant = 0.02,
                       swc.sat = 0.38,
                       psi.e = -0.74e-3,#KPa
                       b = 6.4, 
                       k.sat = 11.4,
                       psi.min = -2,
                       r.l=0.48,
                       lai=8,
                       l.v=1490,
                       r.root=0.0005,
                       p50 = -3,
                       plc.slope = 50
){
        
        psi.soil <- psi.s.func(swc = swc,
                               psi.e = psi.e,#KPa
                               b = b, 
                               swc.sat = swc.sat)
        
        # psi.soil <- max(psi.soil , psi.min)
        
        # calculate max e by plant
        e.plant.func <- function(psi,psi.s,p50,plc.slope){
                KPfnc(psi,plc.slope,p50) * (psi.s - psi)
        }
        
        # get optimal p.psi that max transpiration
        psi.p.opt = optimize(e.plant.func,interval=c(-15,0),maximum=TRUE,
                             psi.s = psi.soil,p50=p50,plc.slope=plc.slope)$maximum
        
        
        # get the K plant at the optimal psi
        k.plant.opt <- k.plant * KPfnc(psi.p.opt,plc.slope,p50)
        e.plant.opt = k.plant.opt * (psi.soil - psi.p.opt)
        # calculate resistent
        R.p = 1 / k.plant.opt
        R.s = 1 / k.s.p.func(swc,swc.sat=swc.sat,psi.e=psi.e,b=b,k.sat=k.sat,
                             r.l=r.l,
                             lai=lai,
                             l.v=l.v,
                             r.root=r.root)
        # -R.p * (psi.soil - psi.min)/ psi.min / (R.p + R.s)
        # get the K plant at the optimal psi
        k.l <- 1 / (R.p + R.s)
        e.max <- k.l*(psi.soil - psi.p.opt)
        
        e.max <- max(0,e.max)
        e.max <- max(e.plant.opt,e.max)
        return(e.max)
}
# 
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
                        r.root=0.0005,
                        p50 = -3,
                        plc.slope = 50
){
        e.psi <- e.max.func(swc = swc,
                            k.plant = k.plant,
                            swc.sat = swc.sat,
                            psi.e = psi.e,#KPa
                            b = b, 
                            k.sat = k.sat,
                            psi.min = psi.min,
                            r.l=r.l,
                            lai=lai,
                            l.v=l.v,
                            r.root=r.root,
                            p50 = p50,
                            plc.slope = plc.slope
        )
        
        e.sat <- e.max.func(swc = swc.sat,
                            k.plant = k.plant,
                            swc.sat = swc.sat,
                            psi.e = psi.e,#KPa
                            b = b, 
                            k.sat = k.sat,
                            psi.min = psi.min,
                            r.l=r.l,
                            lai=lai,
                            l.v=l.v,
                            r.root=r.root,
                            p50 = p50,
                            plc.slope = plc.slope
        )
        
      return(e.psi/e.sat)
}

# curve(e.frac.func(x),from=.01,to=.38)
# e.frac.func <- function(swc,
#                         k.plant = 0.02,
#                         swc.sat = 0.38,
#                         psi.e = -0.74e-3,#KPa
#                         b = 6.4, 
#                         k.sat = 11.4,
#                         psi.min = -2,
#                         r.l=0.48,
#                         lai=8,
#                         l.v=1490,
#                         r.root=0.0005){
#         
#         psi.soil <- psi.s.func(swc = swc,
#                                psi.e = psi.e,#KPa
#                                b = b, 
#                                swc.sat = swc.sat)
#         
#         psi.soil <- max(psi.soil , psi.min)
#         
#         R.p = 1 / k.plant
#         R.s = 1 / k.s.p.func(swc,swc.sat=swc.sat,psi.e=psi.e,b=b,k.sat=k.sat,
#                              r.l=r.l,
#                              lai=lai,
#                              l.v=l.v,
#                              r.root=r.root)
#         # -R.p * (psi.soil - psi.min)/ psi.min / (R.p + R.s)
#         k.l <- 1 / (R.p + R.s)
#         e.max <- k.l*(psi.soil - psi.min)
#         # 
#         R.s.sat <- 1 / k.s.p.func(swc.sat,swc.sat=swc.sat,psi.e=psi.e,b=b,k.sat=k.sat,
#                                   r.l=r.l,
#                                   lai=lai,
#                                   l.v=l.v,
#                                   r.root=r.root)
#         
#         k.l.sat <- 1/(R.p + R.s.sat)
#         
#         e.max.sat <- k.l.sat*(psi.e - psi.min)
#         
#         return(e.max/e.max.sat)
# }
