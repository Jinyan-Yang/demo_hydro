##########################################################################################
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

# Goal: The goal here to to implement the Sperry model at leaf scale to 
# allow simple fitting and calculation with gas exchange data. 


library(shiny)
library(plantecophys)

server <- function(input, output, session) {
  
  # main function to calculate net gain and make plot
  plotSperryt_fun <- function(input){
    
    # ###### example for debuging#####
    # input <- data.frame(kl=4 ,
    #                     psis=-0.01,
    #                     VPD=1.5,
    #                     Jmax=160,
    #                     Vcmax=100,
    #                     sf=10,
    #                     psif=-5,
    #                     Ca = 400)
    
    # plc curve Weilbull
    KPfnc <- function(P, SX, PX,X=50){
      
      P = -P
      PX = -PX
      
      X <- X[1] # when fitting; vector may be passed but X cannot actually vary.
      V <- (X-100)*log(1-X/100)
      p <- (P/PX)^((PX*SX)/V)
      relk <- (1-X/100)^p
      
      return(relk)
    }
    
    # # plant water suuply
    e.supply.old <- function(psil, kl, psis,sf, psif){
      # change k as a function of psi
        kl.use = kl * KPfnc(psil,sf, psif)

      return(kl.use * (psis - psil))
    }
    # # 
    # # plant water suuply in integral way
    e.supply <- function(psil, kl, psis,sf, psif){
    
      x <-  integrate(KPfnc,psil,psis,SX=sf,PX=psif)
       
      y <- kl * x$value
      
      return(y)
    }
    
    # test difference between the two functions
    # psil.vec <- seq(-3,-1,by=0.1)
    # 
    # loss.old.vec <- loss.vec <- c()
    # 
    # for(i in seq_along(psil.vec)){
    #   loss.old.vec[i] <- e.supply.old(psil=psil.vec[i], kl = 1, psis=-0.1,sf=10, psif=-2)
    #   loss.vec[i] <- e.supply(psil=psil.vec[i], kl = 1, psis=-0.1,sf=10, psif=-2)
    # }
    # 
    # 
    # plot(loss.vec~psil.vec,col='red')
    # points(loss.old.vec~psil.vec)
    
    # get gs based on transpiration
    cal.gs.func <- function(psil,kl, psis,VPD,sf, psif,Patm = 101){
      
      gs = e.supply(psil, kl, psis,sf, psif) / VPD * Patm
      
      mmol2mol = 1e-3
      
      return(gs * mmol2mol)
    }
   
    
    # profit based on photosynthesis
    profit.func <- function(psil,kl, psis,VPD,Jmax,Vcmax, sf, psif,Ca,psi.crit){
      
      # carbon gain at given psil
      photo.given.df <- Photosyn(GS = cal.gs.func(psil,kl, psis,VPD, sf, psif),
                                 Jmax = Jmax,
                                 Vcmax = Vcmax,
                                 Ca = Ca, VPD= VPD)
      # carbon gain at critical psil
      photo.max.df <- Photosyn(GS = cal.gs.func(psi.crit,kl, psis,VPD, sf, psif),
                               Jmax = Jmax,
                               Vcmax = Vcmax,
                               Ca = Ca, VPD= VPD)
      
      return(photo.given.df$ALEAF / photo.max.df$ALEAF)
    }
    
    # hydraulic cost
    cost.func <- function(psil,sf, psif,psi.crit){
      # cost at given psil
      tmp1 = 1 - KPfnc(psil,sf, psif)
      # cost at cirtical psil
      tmp2 = 1 - KPfnc(psi.crit,sf, psif)
      
      return(tmp1 / tmp2)
    }
    
    # net gain for optimasation
    opt.func <- function(psil,kl, psis,VPD,Jmax,Vcmax, sf, psif,Ca,psi.crit){
      opt.target <- profit.func(psil,kl, psis,VPD,Jmax,Vcmax, sf, psif,Ca,psi.crit) - cost.func(psil, sf, psif,psi.crit)
      return(opt.target)
    }

    # calculate the critical PsiL as the PSil that gives max transpiration
    pt.value <- optimize(e.supply,c(-7,0),
                         kl=input$kl, psis=input$psis,sf=input$sf,psif=input$psif,
                         maximum = T)$maximum
    
    
    # calculate the oprimal Psil that max net gain
    opti.results <- optim(par = list(psil = -2),
                          fn = opt.func,
                          
                          control = list(fnscale = -1),method= 'Brent',lower = -2.999,upper = 0,
                          
                          kl=input$kl, 
                          psis=input$psis,
                          sf=input$sf,
                          psif=input$psif,
                          Ca = input$Ca,
                          VPD=input$VPD,
                          Jmax=input$Jmax,
                          Vcmax=input$Vcmax,
                          
                          psi.crit = pt.value)
    
    
    
    # plots
    par(mfrow=c(1,1), mar=c(5,5,1,1))
    # plot profit
    
    psil.vec <- seq(pt.value,input$psis,by=0.1)
    
    profit.vec <- c()
    cost.vec <- c()
    for (i in seq_along(psil.vec)){
      profit.vec[i] <- profit.func(psil = psil.vec[i], 
                                   kl=input$kl, 
                                   psis=input$psis,
                                   sf=input$sf,
                                   psif=input$psif,
                                   Ca = input$Ca,
                                   VPD=input$VPD,
                                   Jmax=input$Jmax,
                                   Vcmax=input$Vcmax,
                                   psi.crit = pt.value)
      
      cost.vec[i] <- cost.func(psil.vec[i],sf=input$sf,psif=input$psif,psi.crit = pt.value)
      
    }
    
    
    plot(profit.vec~psil.vec,xlab="LWP (MPa)",
         ylab="Relative gain/loss",
         lwd=2,
         ylim=c(0,1),type='l')
    
    points(cost.vec~psil.vec,lwd=2, lty='dotted',type='l')
    
    points(c(profit.vec-cost.vec)~psil.vec,lwd=2, lty=1,type='l',col='red')
    
    abline(v = opti.results$par,lty='dashed',col='grey',lwd=2)
    
    legend('bottom',legend=c('Carbon gain','Hydraulic cost','Net gain'),lty=c('solid','dotted','solid'),bty='n',col=c('black','black','red'))
    
    # curve(profit.func(psil = x, 
    #                   kl=input$kl, 
    #                   psis=input$psis,
    #                   sf=input$sf,
    #                   psif=input$psif,
    #                   Ca = input$Ca,
    #                   VPD=input$VPD,
    #                   Jmax=input$Jmax,
    #                   Vcmax=input$Vcmax,
    #                   psi.crit = pt.value), 
    #       xlab="LWP (MPa)",
    #       ylab="Relative gain/loss",
    #       from=pt.value, to=-0.05, 
    #       lwd=2,
    #       ylim=c(0,1))
    # # plot cost
    # curve(cost.func(x,sf=input$sf,psif=input$psif,psi.crit = pt.value),
    #       add=TRUE, lwd=1, lty='dotted',ylim=c(0,1))
    # 
    # # add optima
    # abline(v = opti.results$par,lty='dashed',col='grey',lwd=2)
    # 
    # legend('bottomleft',legend=c('Carbon gain','Hydraulic cost'),lty=c('solid','dotted'),bty='n')

  }
  
  output$Sperryplot <- renderPlot(
    plotSperryt_fun(input)
  )
}


ui <- fluidPage(
  titlePanel("Leaf Level Sperry Model "),
  
  plotOutput("Sperryplot"),
  

  fluidRow(
    # environmental drivers
    column(3,
           sliderInput("psis", "Soil water potential (MPa):",
                       min=-3, max=0, value=-0.01, step=0.05),
           sliderInput("VPD", "Vapour pressure deficit (kPa):",
                       min=0.5, max=5, value=1.5, step=0.05),
           sliderInput("Ca", "Ca (ppm):",
                       min=50, max=600, value=400, step=50)
           
           
    ),
    # hydraulics
    column(3,
           sliderInput("sf", "PLC Slope:",
                       min=0.5, max=50, value=10, step=0.5),
           sliderInput("psif", "PLC P50 (MPa):",
                       min=-8, max=-0.5, value=-5, step=0.25),
           sliderInput("kl", "Max hydraulic conductance (mmol per m2 per s):",
                       min=0.5, max=5, value=4, step=0.05)
    ),
    # photosynthtic
    column(3,
           sliderInput("Jmax", "Jmax:",
                       min=50, max=250, value=160, step=10),
           sliderInput("Vcmax", "Vcmax:",
                       min=25, max=150, value=100, step=10)
    )
  )
)
  


library(shiny)
x <- shinyApp(ui, server)


