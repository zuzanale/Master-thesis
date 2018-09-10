# import libraries
library(shiny)
library(stringi)
library(stringr)
library(gets) 
source("plot.my.R")
##### simulation part

load("shiny.Rdata")
lts2.zoo=log(ts2.zoo)
indicator1 = ts(lts2.zoo,frequency =7) 
lts2.11.zoo=log(ts2.11.zoo)
indicator1.11 = ts(lts2.11.zoo,frequency =7) 
lts2.21.zoo=log(ts2.21.zoo)
indicator1.21 = ts(lts2.21.zoo,frequency =7) 
lts2.22.zoo=log(ts2.22.zoo)
indicator1.22 = ts(lts2.22.zoo,frequency =7) 
lts3.zoo=log(ts3.zoo)
indicator2 = ts(lts3.zoo,frequency =7) 
lts3.11.zoo=log(ts3.11.zoo)
lts3.11.zoo[lts3.11.zoo==-Inf]=0 #replace -inf by 0
indicator2.11 = ts(lts3.11.zoo,frequency =7) 
lts3.12.zoo=log(ts3.12.zoo)
lts3.12.zoo[lts3.12.zoo==-Inf]=0 #replace -inf by 0
indicator2.12 = ts(lts3.12.zoo,frequency =7) 
lts3.13.zoo=log(ts3.13.zoo)
lts3.13.zoo[lts3.13.zoo==-Inf]=0 #replace -inf by 0
indicator2.13 = ts(lts3.13.zoo,frequency =7) 
lts4.zoo=log(ts4.zoo)
lts4.zoo[lts4.zoo==-Inf]=0 #replace -inf by 0
indicator4 = ts(lts4.zoo,frequency =7) 
lts5.zoo=log(ts5.zoo)
lts5.zoo[lts5.zoo==-Inf]=0 #replace -inf by 0
indicator5 = ts(lts5.zoo,frequency =7) 

load("indicator_results.RData")

####corrections of different ascii
# Praparations (could be put into global.R) ------------
choices <-c("indicator1","indicator1.11","indicator1.21",
                     "indicator1.22",
                     "indicator2","indicator2.11","indicator2.12",
                     "indicator2.13","indicator4",
                     "indicator5")

# Replace name place holder with the actual one
names(choices)[1] <- paste0(
  "indicator1",
  stri_dup(intToUtf8(160), 1)) # Replace 6 with the desired number
  #"I want multiple spaces here")
names(choices)[2] <- paste0(
  "indicator1.11",
  stri_dup(intToUtf8(160), 1))
names(choices)[3] <- paste0(
  "indicator1.21",
  stri_dup(intToUtf8(160), 1))
names(choices)[4] <- paste0(
  "indicator1.22",
  stri_dup(intToUtf8(160), 1))
names(choices)[5] <- paste0(
  "indicator2",
  stri_dup(intToUtf8(160), 1))
names(choices)[6] <- paste0(
  "indicator2.11",
  stri_dup(intToUtf8(160), 1))
names(choices)[7] <- paste0(
  "indicator2.12",
  stri_dup(intToUtf8(160), 1))
names(choices)[8] <- paste0(
  "indicator2.13",
  stri_dup(intToUtf8(160), 1))
names(choices)[9] <- paste0(
  "indicator4",
  stri_dup(intToUtf8(160), 1))
names(choices)[10] <- paste0(
  "indicator5",
  stri_dup(intToUtf8(160), 1))

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
  titlePanel("Development of scores"),
  
  #This creates a layout with a left sidebar and main section
 sidebarLayout(
   sidebarPanel(
              
                                                       
     selectInput(inputId ="indicator", "Choose an indicator or its breakdown",
                            choices=choices,multiple=FALSE),
                
                sliderInput(inputId='timerange',label="Choose time period",
                            min=0,max=275,
                            value=c(0,275)),
                selectInput(inputId ="trans", "Transformation to apply",
                            choices=c("level","logarithm"),multiple=FALSE),
                selectInput(inputId ="NAs", "Missing values handling (before the transformation)",
                            choices=c("no handling","interpolation","mean",
                                      "median",
                                      "zero"),multiple=FALSE),
               conditionalPanel(
                 condition = "input.NAs == 'interpolation'",
                selectInput("interpolation.var", "Interpolation:",
                              choices=c("linear","spline", "stine"),multiple=FALSE)
               )),
                #numericInput(inputId='pvalue.input',label="Choose significance threshold (p.value)",
               #              value = 1/length(input$indicator.input),step=0.01)
                
   
   
    #beginning of sidebar section
   #this is for the selection of which time series we want to analyse
    #for which time rnage we want to analyse
    
    #usually includes inputs
    #sidebarPanel(),
    
    #beginning of main section
  #  mainPanel(
     
   #   plotOutput("plot1"),
   #   br(), br(),
  #    plotOutput("plot2"),
  #    br(), br(),
  #    plotOutput("plot3"),
  #    br(), br(),
  #    plotOutput("plot4"),
  #  )
   
   mainPanel(
     #plotOutput("plot1"),
     #   br(), br(),
     fluidRow(
       column(width =6,
              h2("Parent Indicator"),
              plotOutput("plot1")),
       column(width =6,
              h2("Breakdown 1"),
              plotOutput("plot2")),
       column(width =6,
              h2("Breakdown 2"),
              plotOutput("plot3")),
       column(width =6,
              h2("Breakdown 3"),
              plotOutput("plot4"))
     ))
)
    #tableOutput()
  
  
  #Close the UI definition
))

server <- shinyServer(function(input, output) {
  
  output$plot1 <- renderPlot({
    range=input$timerange[1]:input$timerange[2]
   # browser()
    time_series=eval(parse(text = input$indicator))
    time_series= time_series[range]
    time_series.copy=time_series
    #if interpolation is selected for the missing data
    require(imputeTS)
    #browser()
    if(input$NAs!="no handling"){
      time_series.copy=time_series
      if(input$NAs=="mean"){
        time_series[is.na(time_series)]=mean(time_series,na.rm=T)
      }else if(input$NAs=="median"){
        time_series[is.na(time_series)]=median(time_series,na.rm=T)
      }else if(input$NAs=="zero"){
        time_series[is.na(time_series)]=0
      }else if(input$NAs=="interpolation"){
        time_series=na.interpolation(time_series, option = input$interpolation.var)
      }
      time_series.copy=time_series
    }

    if (input$trans=="logarithm"){
      ltime_series=log(time_series)
    # ltime_series[is.na(ltime_series)]=
      
      if(input$NAs!="no handling"){
        plot(as.ts(ltime_series),ylab="",col='red',lty=2)
        lines(log(time_series.copy),col='black')
      }
      if(input$NAs=="mean"){
        ltime_series[is.na(ltime_series)]=mean(ltime_series,na.rm=T)
      }else if(input$NAs=="median"){
        ltime_series[is.na(ltime_series)]=median(ltime_series,na.rm=T)
      }else if(input$NAs=="zero"){
        ltime_series[is.na(ltime_series)]=0
      }else if(input$NAs=="interpolation"){
        ltime_series=na.interpolation(ltime_series, option = input$interpolation.var)
      }
      time_series.copy=ltime_series
      plot(as.ts(ltime_series),ylab="",col='red',lty=2)
      lines((time_series.copy),col='black')
      #,ylim=c(min(input$indicator,na.rm=T),max(input$indicator,na.rm=T)))
    }else{
      if(input$NAs!="no handling"){
        plot(as.ts(time_series),ylab="",col='red',lty=2)
        lines((time_series.copy),col='black')
      }
      plot(as.ts(time_series),ylab="",col='red',lty=2)
      lines((time_series.copy),col='black')
    }
    
   
  })
  output$plot2 <- renderPlot({
    range=input$timerange[1]:input$timerange[2]
    # browser()
    time_series=eval(parse(text = input$indicator))
    time_series= time_series[range]
    time_series.copy=time_series
    
    if(time_series==indicator1){
      time_series=indicator1.11
      #browser()
    }else if(input$indicator=="indicator2"){
      time_series=indicator2.11
    }
    
    #if interpolation is selected for the missing data
    require(imputeTS)
    #browser()
    if(input$NAs!="no handling"){
      time_series.copy=time_series
      if(input$NAs=="mean"){
        time_series[is.na(time_series)]=mean(time_series,na.rm=T)
      }else if(input$NAs=="median"){
        time_series[is.na(time_series)]=median(time_series,na.rm=T)
      }else if(input$NAs=="zero"){
        time_series[is.na(time_series)]=0
      }else if(input$NAs=="interpolation"){
        time_series=na.interpolation(time_series, option = input$interpolation.var)
      }
      time_series.copy=time_series
    }
    
    if (input$trans=="logarithm"){
      ltime_series=log(time_series)
      # ltime_series[is.na(ltime_series)]=
      
      if(input$NAs!="no handling"){
        plot(as.ts(ltime_series),ylab="")
      #  lines(log(time_series.copy),col='black')
      }
      if(input$NAs=="mean"){
        ltime_series[is.na(ltime_series)]=mean(ltime_series,na.rm=T)
      }else if(input$NAs=="median"){
        ltime_series[is.na(ltime_series)]=median(ltime_series,na.rm=T)
      }else if(input$NAs=="zero"){
        ltime_series[is.na(ltime_series)]=0
      }else if(input$NAs=="interpolation"){
        ltime_series=na.interpolation(ltime_series, option = input$interpolation.var)
      }
      time_series.copy=ltime_series
      plot(as.ts(ltime_series),ylab="")
   #   lines((time_series.copy),col='black')
      #,ylim=c(min(input$indicator,na.rm=T),max(input$indicator,na.rm=T)))
    }else{
      if(input$NAs!="no handling"){
        plot(as.ts(time_series),ylab="")
       # lines((time_series.copy),col='black')
      }
      plot(as.ts(time_series),ylab="")
    #  lines((time_series.copy),col='black')
    }
    
    
  })
  
  output$plot3 <- renderPlot({
    range=input$timerange[1]:input$timerange[2]
    # browser()
    time_series=eval(parse(text = input$indicator))
    time_series= time_series[range]
    time_series.copy=time_series
    
    if(time_series==indicator1){
      time_series=indicator1.21
    }else if(input$indicator=="indicator2"){
      time_series=indicator2.12
    }
    
    #if interpolation is selected for the missing data
    require(imputeTS)
    #browser()
    if(input$NAs!="no handling"){
      time_series.copy=time_series
      if(input$NAs=="mean"){
        time_series[is.na(time_series)]=mean(time_series,na.rm=T)
      }else if(input$NAs=="median"){
        time_series[is.na(time_series)]=median(time_series,na.rm=T)
      }else if(input$NAs=="zero"){
        time_series[is.na(time_series)]=0
      }else if(input$NAs=="interpolation"){
        time_series=na.interpolation(time_series, option = input$interpolation.var)
      }
      time_series.copy=time_series
    }
    
    if (input$trans=="logarithm"){
      ltime_series=log(time_series)
      # ltime_series[is.na(ltime_series)]=
      
      if(input$NAs!="no handling"){
        plot(as.ts(ltime_series),ylab="")
        #  lines(log(time_series.copy),col='black')
      }
      if(input$NAs=="mean"){
        ltime_series[is.na(ltime_series)]=mean(ltime_series,na.rm=T)
      }else if(input$NAs=="median"){
        ltime_series[is.na(ltime_series)]=median(ltime_series,na.rm=T)
      }else if(input$NAs=="zero"){
        ltime_series[is.na(ltime_series)]=0
      }else if(input$NAs=="interpolation"){
        ltime_series=na.interpolation(ltime_series, option = input$interpolation.var)
      }
      time_series.copy=ltime_series
      plot(as.ts(ltime_series),ylab="")
      #   lines((time_series.copy),col='black')
      #,ylim=c(min(input$indicator,na.rm=T),max(input$indicator,na.rm=T)))
    }else{
      if(input$NAs!="no handling"){
        plot(as.ts(time_series),ylab="")
        # lines((time_series.copy),col='black')
      }
      plot(as.ts(time_series),ylab="")
      #  lines((time_series.copy),col='black')
    }
    
    
  })
  output$plot4 <- renderPlot({
    range=input$timerange[1]:input$timerange[2]
    # browser()
    time_series=eval(parse(text = input$indicator))
    time_series= time_series[range]
    time_series.copy=time_series
    
    if(time_series==indicator1){
      time_series=indicator1.22
    }else if(input$indicator=="indicator2"){
      time_series=indicator2.13
    }
    
    #if interpolation is selected for the missing data
    require(imputeTS)
    #browser()
    if(input$NAs!="no handling"){
      time_series.copy=time_series
      if(input$NAs=="mean"){
        time_series[is.na(time_series)]=mean(time_series,na.rm=T)
      }else if(input$NAs=="median"){
        time_series[is.na(time_series)]=median(time_series,na.rm=T)
      }else if(input$NAs=="zero"){
        time_series[is.na(time_series)]=0
      }else if(input$NAs=="interpolation"){
        time_series=na.interpolation(time_series, option = input$interpolation.var)
      }
      time_series.copy=time_series
    }
    
    if (input$trans=="logarithm"){
      ltime_series=log(time_series)
      # ltime_series[is.na(ltime_series)]=
      
      if(input$NAs!="no handling"){
        plot(as.ts(ltime_series),ylab="")
        #  lines(log(time_series.copy),col='black')
      }
      if(input$NAs=="mean"){
        ltime_series[is.na(ltime_series)]=mean(ltime_series,na.rm=T)
      }else if(input$NAs=="median"){
        ltime_series[is.na(ltime_series)]=median(ltime_series,na.rm=T)
      }else if(input$NAs=="zero"){
        ltime_series[is.na(ltime_series)]=0
      }else if(input$NAs=="interpolation"){
        ltime_series=na.interpolation(ltime_series, option = input$interpolation.var)
      }
      time_series.copy=ltime_series
      plot(as.ts(ltime_series),ylab="")
      #   lines((time_series.copy),col='black')
      #,ylim=c(min(input$indicator,na.rm=T),max(input$indicator,na.rm=T)))
    }else{
      if(input$NAs!="no handling"){
        plot(as.ts(time_series),ylab="")
        # lines((time_series.copy),col='black')
      }
      plot(as.ts(time_series),ylab="")
      #  lines((time_series.copy),col='black')
    }
    
    
  })
 
  #Close the server definition
})

##############
# Launch App #
##############

#generic line that launches the app
shinyApp(ui = ui, server = server)