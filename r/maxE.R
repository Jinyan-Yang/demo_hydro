#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plantecophys)
server <- function(input, output, session) {
    
    
    maxEPlot_fun <- function(input){
            
            # plc 
            # p = plant psi
            # sx is the slope of plc
            # px= p50
            
            KPfnc <- function(P, SX, PX,psi.cruit=-50,crit.psi.test = FALSE) {
                    if(crit.psi.test == TRUE){
                            if(P > psi.cruit){
                                    X <- 50
                                    V <- (X - 100) * log(1 - X/100)
                                    p <- (P/PX)^((PX * SX)/V)
                                    relk <- (1 - X/100)^p
                                    
                            }else{
                                    0.1+0.0001*P
                            }
                    }else{
                            X <- 50
                            V <- (X - 100) * log(1 - X/100)
                            p <- (P/PX)^((PX * SX)/V)
                            relk <- (1 - X/100)^p
                            
                    }
                    return(1-relk)
            }
            
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
                                swc.sat = 0.38,
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
        # plot####
        
        par(mfrow=c(2,2), mar=c(5,5,1,1))
        # plot psi
        curve((psi.s.func(swc=x,b = input$b,
                          swc.sat = 0.38,
                          psi.e = -0.74e-3)),
              from = .01,to=.38,
              xlab="VWC",
              ylab=expression(psi[soil]~('MPa')),
              lwd=2,
              ylim=c(-10,0))
        # plot soil conductivity
        curve((k.soil.func(swc=x,b = input$b,
                          swc.sat = 0.38,
                          psi.e = -0.74e-3,
                          k.sat = input$k.sat)),
              from = .01,to=.38,
              xlab="VWC",
              ylab=expression(k[soil]~('mol MPa'^-1~s^-1^m-2)),
              lwd=2)
        # plot fractional reduction of conductance
        curve((k.s.p.func(swc=x,
                          r.l=input$r.l,
                          lai=input$lai,
                          l.v=input$l.v,
                          r.root=input$r.root,
                          swc.sat = 0.38,
                          psi.e = -0.74e-3,#KPa
                          b = input$b, 
                          k.sat = input$k.sat) /
                   k.s.p.func(0.38,
                              r.l=input$r.l,
                              lai=input$lai,
                              l.v=input$l.v,
                              r.root=input$r.root,
                              swc.sat = 0.38,
                              psi.e = -0.74e-3,#KPa
                              b = input$b, 
                              k.sat = input$k.sat)), from = .01,to=.38,
              xlab="VWC",
              ylab=expression(K[soil]*'/'*K[soil.sat]),
              lwd=2)
        
        # plot max E refduction
        
        swc.vec = seq(.01,.338,by=0.001)
        frac.vec <- c()
        for (i in seq_along(swc.vec)) {
                frac.vec[i] <- e.frac.func(swc=swc.vec[i],
                                           r.l=input$r.l,
                                           lai=input$lai,
                                           l.v=input$l.v,
                                           r.root=input$r.root,
                                           k.plant = input$k.plant,
                                           swc.sat = 0.38,
                                           psi.e = -0.74e-3,#KPa
                                           b = input$b, 
                                           k.sat = input$k.sat,
                                           psi.min = input$psi.min,
                                           p50 = input$p50,
                                           plc.slope = input$plc.slope)
                
        }
        
        plot(frac.vec~swc.vec,type='l',lwd=2,
             xlab="VWC",ylab=expression(E*'/'*E[max])
             )
        # curve(e.frac.func(swc=x,
        #                    r.l=input$r.l,
        #                    lai=input$lai,
        #                    l.v=input$l.v,
        #                    r.root=input$r.root,
        #                    k.plant = input$k.plant,
        #                    swc.sat = 0.38,
        #                    psi.e = -0.74e-3,#KPa
        #                    b = input$b, 
        #                    k.sat = input$k.sat,
        #                    psi.min = input$psi.min,
        #                    p50 = input$p50,
        #                    plc.slope = input$plc.slope), from = .01,to=.38,
        #       xlab="VWC",
        #       ylab=expression(E*'/'*E[max]),
        #       lwd=2)
        
        
    }
    
    output$maxEPlot <- renderPlot(
        maxEPlot_fun(input),width = 700, height = 700*0.618
    )
}


ui <- fluidPage(
    titlePanel("Using soil hydraulics to predict plant water aviability"),
    
    
    fluidRow(
        column(5,
               sliderInput("lai", "Leaf area index (m2 m-2):",
                           min=0.5, max=8, value=0.5, step=0.5),
               sliderInput("l.v", "Root length density (m m-3):",
                           min=500, max=5000, value=5000, step=100),
               sliderInput("r.root", "Mean radius of water absorbing roots (m):",
                           min=0.0001, max=0.001, value=0.001, step=0.0001),
               sliderInput("k.plant", "plant hydraulic conductance (mol m–2 MPa-1 s–1):",
                           min=0.02, max=0.1, value=0.02, step=0.01),
               sliderInput("r.l", "Root radius (m):",
                           min=0.1, max=2, value=0.5, step=2),
               sliderInput("psi.min", "Plant critical water potential (MPa):",
                           min=-8, max=-1, value=-8, step=0.5)
               
        ),
        column(5,
               sliderInput("b", "soil retention parameter:",
                           min=2, max=12, value=2, step=.5),
               sliderInput("k.sat", "soil saturated hydraulic conductivity (mol m–1 s–1 MPa–1):",
                           min=5, max=250, value=250, step=5),
               
               sliderInput("p50", "Plant PSI at 50% PLC:",
                           min=-10, max=-.5, value=-3, step=.5),
               sliderInput("plc.slope", "Slope of PLC curve:",
                           min=1, max=100, value=50, step=5)
        )
    ),
    plotOutput("maxEPlot")
    
)

library(shiny)
x <- shinyApp(ui, server)


