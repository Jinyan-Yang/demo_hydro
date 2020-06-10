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
    
    
    tuzetplot_fun <- function(input){
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
        # ####
        
        par(mfrow=c(1,2), mar=c(5,5,1,1))

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
                              k.sat = input$k.sat)), from = .1,to=.38,
              xlab="VWC",
              ylab=expression(K[soil]*'/'*K[soil.sat]),
              lwd=2)
        
        
        curve((e.frac.func(swc=x,
                           r.l=input$r.l,
                           lai=input$lai,
                           l.v=input$l.v,
                           r.root=input$r.root,
                           k.plant = input$k.plant,
                           swc.sat = 0.38,
                           psi.e = -0.74e-3,#KPa
                           b = input$b, 
                           k.sat = input$k.sat,
                           psi.min = input$psi.min)), from = .1,to=.38,
              xlab="VWC",
              ylab=expression(E*'/'*E[max]),
              lwd=2)
        
        
    }
    
    output$tuzetplot <- renderPlot(
        tuzetplot_fun(input),width = 400*2, height = 400*0.618,
    )
}


ui <- fluidPage(
    titlePanel("Using soil hydraulics to predict plant water aviability"),
    
    plotOutput("tuzetplot"),
    
    fluidRow(
        column(5,
               sliderInput("lai", "Leaf area index (m2 m-2):",
                           min=0.5, max=8, value=4, step=0.5),
               sliderInput("l.v", "Root length density (m m-3):",
                           min=500, max=5000, value=1500, step=100),
               sliderInput("r.root", "Mean radius of water absorbing roots (m):",
                           min=0.0001, max=0.001, value=0.0005, step=0.0001),
               sliderInput("k.plant", "plant hydraulic conductance (mol m–2 MPa-1 s–1):",
                           min=0.02, max=0.1, value=0.02, step=0.01),
               sliderInput("r.l", "Root radius (m):",
                           min=0.1, max=2, value=0.5, step=0.1),
               sliderInput("psi.min", "Plant critical water potential (MPa):",
                           min=-8, max=-1, value=-2, step=0.5)
               
        ),
        column(5,
               sliderInput("b", "soil retention parameter:",
                           min=2, max=12, value=6, step=1),
               sliderInput("k.sat", "soil saturated hydraulic conductivity (mol m–1 s–1 MPa–1):",
                           min=5, max=250, value=10, step=5)
        )
    )
    
)

library(shiny)
x <- shinyApp(ui, server)


