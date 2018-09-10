# import libraries
library(shiny)
library(shinythemes)
library(stringi)
library(stringr)
library(gets)
library(tsoutliers)
library(stsm)
source("datagen.outlier.R")
##### simulation part

##########
# SERVER #
##########

#generic line initiating the SERVER 


##################
# User Interface #
##################

#generic line initiating the UI
ui <- shinyUI(fluidPage(theme = shinytheme("superhero"),
  
  #Add a title
  titlePanel("Outlier generation"),
  
  #This creates a layout with a left sidebar and main section
 sidebarLayout(
   sidebarPanel(numericInput('var_param_eps',"Choose the variance for epsilon",
                             value=10,step=0.001),
                sliderInput(inputId='outlier_eps',label="Additive outlier index (positive)",
                            min=0,max=100,
                            value=15),
                sliderInput(inputId='outlier_eps2',label="Additive outlier index (negative)",
                            min=0,max=100,
                            value=40),
                numericInput('var_param_eta',"Choose the variance for eta",
                                          value=6,step=0.001),
                sliderInput("outlier_eta", label = "Structural (level) break index (positive)", min = 0, 
                            max = 100, value =50),
                numericInput('var_param_zeta',"Choose the variance for zeta",
                              value=0.0001,step=0.001),
                sliderInput("outlier_zeta", label = "Trend change index (negative)", min = 0, 
                            max = 100, value =10),
              numericInput('var_param_omega',"Choose the variance for omega",
                       value=8,step=0.001),
              sliderInput("outlier_omega", label = "Seasonality change index (negative)", min = 0, 
                          max = 100, value =70)
              
   
               ),
                #numericInput(inputId='pvalue.input',label="Choose significance threshold (p.value)",
               #              value = 1/length(input$indicator.input),step=0.01)
                
   
   
    #beginning of sidebar section
   #this is for the selection of which time series we want to analyse
    #for which time rnage we want to analyse
    
    #usually includes inputs
    #sidebarPanel(),
    
    #beginning of main section
    mainPanel(
      plotOutput("plotsim")
    )
    #plot output
)
    #tableOutput()
  
  
  #Close the UI definition
))

server <- shinyServer(function(input, output) {
  output$plotsim <- renderPlot({
    pars <- c(var1 =input$var_param_eps, var2 = input$var_param_eta ,var3 =input$var_param_zeta, var4 = input$var_param_omega, a01 = 25)
    m <- stsm::stsm.model(model = "BSM", y = ts(seq(100),freq=7), #here a multinomial frequency in ts would be nice to have
                          pars = pars, nopars = NULL)
    ss <- char2numeric(m)
    if(input$outlier_eps==0){
      AOplusmag=0
    }else{
      AOplusmag=6
    }
    if(input$outlier_eps2==0){
      AOminusmag=0
    }else{
      AOminusmag=10
    }
    if(input$outlier_eta==0){
      LSplusmag=0
    }else{
      LSplusmag=8
    }
    if(input$outlier_zeta==0){
      TCminusmag=0
    }else{
      TCminusmag=8
    }
    if(input$outlier_omega==0){
      SAOminusmag=0
    }else{
      SAOminusmag=7
    }
    set.seed(123)
    y_m=datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                    n0 = 70, old.version = TRUE,
                    AOplus=c(input$outlier_eps),AOplusmag=AOplusmag, AOminus=c(input$outlier_eps2),AOminusmag=AOminusmag,
                    LSplus=c(input$outlier_eta),LSplusmag=LSplusmag,
                    TCminus=c(input$outlier_zeta), TCminusmag=TCminusmag,SAOminus=c(input$outlier_omega), SAOminusmag=SAOminusmag
    )$data
    set.seed(123)
    y= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                       n0 = 70, old.version = TRUE,
                       AOplus=c(input$outlier_eps),AOplusmag=0, AOminus=c(input$outlier_eps2),AOminusmag=0,
                       LSplus=c(input$outlier_eta),LSplusmag=0,
                       TCminus=c(input$outlier_zeta), TCminusmag=0,SAOminus=c(input$outlier_omega), SAOminusmag=0
    )$data
    plot(y_m,col="blue",ylab="",ylim=c(min(y,y_m)-2,max(y,y_m)))
    lines(y, pch=22, lty=2)
    legend(1,max(y,y_m), legend=c("BSM", "BSM with outliers"),
           col=c("blue","black"), lty=1:2, cex=0.8)
  })
  
  
  #Close the server definition
})

##############
# Launch App #
##############

#generic line that launches the app
shinyApp(ui = ui, server = server)