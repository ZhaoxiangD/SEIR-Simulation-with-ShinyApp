library(deSolve)
library(tidyverse)
library(ggplot2)
library(shiny)
library(shinythemes)
library(tableone) # for prettier tables to display
library(kableExtra)
library(plotly)
library(shinydashboard)
library(dashboardthemes)
library('ggpubr') 
library('ggsci')
library(shinyjs)
library(Polychrome)

load("contacts_china.Rdata")
wuhanpop <- read.csv('wuhanpop.csv')

# calculate MSE fundtion
mse <- function(lambda, data, func, times, initial_state, params){
  #` Returns MSE of estimation given lambda compared with original data
  #` @param lambda, numeric value contact rate
  #` @param data_toy, numeric vector of original data
  #` @param func, function of model
  #` @param times, observing time
  #` @param initial state 
  #` @param params_toy, list of parameters, including:
  #`  m(number of age groups)
  #   epsilon(latency rate)
  #   lambda(contact rate)
  #   gamma(recovery rate)
  #` @returns MSE of estimation
  predict <- ode(initial_state, times, func, c(lambda=lambda, params))
  mse <- sum((data-predict[,-1])^2)
  return(mse)
}

# accumulate dataset by frame function
accumulate_by <- function(dat, var) {
  var <- lazyeval::f_eval(var, dat)
  lvls <- plotly:::getLevels(var)
  dats <- lapply(seq_along(lvls), function(x) {
    cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
  })
  dplyr::bind_rows(dats)
}

# SEIR model function ------------
SEIR <- function(time, Y, params){
  #` Returns states of each compartment(S,E,I,R)
  #` @param time, a numeric value, observing time
  #` @param Y, a numeric vector or matrix, containing initial states in each 
  #   compartments stratified by age
  #` @param params, list of parameters, including:
  #   m(number of age groups, set to 1 by default),
  #   epsilon(latency rate)
  #   lambda(contact rate)
  #   gamma(recovery rate)
  #   mu(ground death rate)
  #` @return a list that records states of change by time in compartments (S,E,I,R)
  
  # create a vector to store each compartments' changes by time
  
  dY <- numeric(length(Y))
  with(params,{
    if (!is.matrix(lambda)){
      lambda <- as.matrix(lambda)
    }
    for (i in 1:m){
      dY[i] <- mu[i] - lambda[,i] %*% Y[2*m + seq(1:m)] * Y[i] - mu[i] * Y[i]
      dY[m+i] <- lambda[,i] %*% Y[2*m + seq(1:m)] * Y[i] - (epsilon+mu[i]) * Y[m+i]
      dY[2*m+i] <- epsilon * Y[m+i] - (mu[i]+gamma) * Y[2*m+i]
      dY[3*m+i] <- gamma * Y[2*m+i] - mu[i] * Y[3*m+i]
    }
    return(list(dY))
  })
}


# Age-structured SEIR function----------
SEIR_age <- function(time, Y, params){
  #` Returns states of each compartment(S,E,I,R)
  #` @param time, a numeric value, observing time
  #` @param Y, a numeric vector or matrix, containing initial states in each 
  #   compartments stratified by age
  #` @param params, list of parameters, including:
  #   m(number of age groups, set to 1 by default),
  #   alpha(proportion of infectiouness)
  #   beta(transmission rate)
  #   C(contacts of age group j made by age group i)
  #   dL(average incubation period)
  #   dI(average duration of infection)
  #   kappa(daily probability of an exposed individual becoming infectious)
  #   rho(proportion that infected cases are clinical), 0.4 if age<20, 0.8 otherwise
  #   gamma(recovery rate)
  #` @return a list that records states of change by time in compartments (S,E,I,R)
  
  # create a vector to store each compartments' changes by time
  dY <- numeric(length(Y))
  with(params,{
    # convert lambda to a numeric matrix for linear operation 
    for (i in 1:m){
      # dS/dt
      dY[i] <- -beta * Y[i] * C[i,] %*% Y[2*m + seq(1:m)] - alpha * beta * Y[i] * C[i,] %*% Y[3*m + seq(1:m)]
      # dE/dt
      dY[m+i] <- -kappa * Y[m+i] + beta * Y[i] * C[i,] %*% Y[2*m + seq(1:m)] + alpha * beta * Y[i] * C[i,] %*% Y[3*m + seq(1:m)]
      # dIc/dt
      dY[2*m+i] <- rho[i] * kappa * Y[m+i] - gamma * Y[2*m+i]
      # dIsc/dt (Infectious individuals that are sub-clinical)
      dY[3*m+i] <- (1 - rho[i]) * kappa * Y[m+i] - gamma * Y[3*m+i]
      # dR/dt (infectious individuals that are clinical)
      dY[4*m+i] <- gamma * (Y[2*m+i] + Y[3*m+i])
    }
    return(list(dY))
  })
}


# ui-----------------------
### side panel---------
sidebar <- dashboardSidebar(
  # ui for model prediction plot
  conditionalPanel(condition="input.tabselected==1",
                   fileInput('dataset', 'Input dataset'),
                   #choose Model
                   selectInput('modelInput', 'Choose Model',
                               choices = c("SEIR", "SEIcIscR"),
                               selected = "SEIR"),
                   #Choose parameters 
                   selectInput('exampleInput', 'Choose Model parameters',
                               choices = c("Normal", "Example1", 'Example2', 'Example3'),
                               selected = 'Normal'),
                   br(),
                   # change parameters
                   uiOutput("epsilon"),
                   uiOutput('lambda'),
                   uiOutput('gamma'),
                   uiOutput('mu')
  ),
  # ui for phase plane plot
  conditionalPanel(condition="input.tabselected==2",
                   br(),
                   # slider heads
                   splitLayout(fluidRow(p(strong("Scenario 1")), align="center"),
                               fluidRow(p(strong("Scenario 2")), align="center")),
                   uiOutput("epsilon2"),
                   uiOutput('lambda2'),
                   uiOutput('gamma2'),
                   uiOutput('mu2')
  )
)

# app head
header <- dashboardHeader(title="SEIR Simulation App")

### main panel----------
body <- dashboardBody(
  useShinyjs(),
  ### changing theme
  mainPanel( width = '200%',
    tabsetPanel(
      # prediction plot
      tabPanel("Trajectories", value = 1,
               fluidRow(box(sliderInput("time", "Time", min = 0, max = 1000, value = 500), width = 12)),
               uiOutput('plotui'),
               uiOutput('standplotui')
               ),
      #ui output
      
      # phase plane plot
      tabPanel("Phase Planes and Equilibrium", value = 2,
               fluidRow(box(plotOutput('phaseplot')), box(plotOutput('phaseplotSR'))),
               uiOutput('tableui')
               ),
      tabPanel("Model and Theoretical Analysis", 
               uiOutput('pdfui')
               ),
      id = "tabselected")
  )
)

ui <- dashboardPage(header, sidebar, body, skin = 'black')

# server--------------------------
server <- function(input, output, session) {
  
  # UI for PDF output
  output$pdfui <- renderUI({
    if (input$modelInput == 'SEIcIscR'){
      tags$iframe(style="height:650px; width:100%; scrolling=yes", 
                  src="seiir_latex.pdf")
    } else {
      tags$iframe(style="height:650px; width:100%; scrolling=yes", 
                  src="seir_latex.pdf")
    }
  })
  

  # UI for plot change by model
 output$plotui <- renderUI({
   if (input$modelInput == 'SEIcIscR'){
     fluidRow(box(plotlyOutput('ageplot'), width = 12,height = 800))
   } else {
     fluidRow(box(plotlyOutput('trajplot'), width = 12))
     }
 })
 
 # plot2 in SEIr
 output$standplotui <- renderUI({
   if (input$modelInput == 'SEIcIscR'){
     div()
   } else {
     fluidRow(box(plotOutput('standplot'), width = 12))
   }
 })
 
 # eq table in tab2
 output$tableui <- renderUI({
   if (input$modelInput == 'SEIR'){
     fluidRow(box(htmlOutput('eqtable'), width = 12))
   } else {
     div()
   }
 })
  # Intial SEIR state
  initial_state <- c(S=0.999, E=0.001, I=0, R=0)
  m_g = 16
  prop_age = wuhanpop$propage
  E0 = Isc0 = R0 = rep(0, m_g)
  Ic0 = rep(0.0002, m_g)
  S0 = prop_age - (E0 + Isc0 + Ic0 + R0)
  initial_state_age <- c(S0,E0,Ic0,Isc0,R0) # intial state for SEIcIscR
  ### global parameters(intial states, time)---------------
  
  # Time 
  #times <- 0:1000
  times_reactive <- reactive({
    req(input$time)
    0:input$time
  })
  
  # augmented time
  times_reactive_stand <- reactive({
    req(input$time)
    0:(100*input$time)
  })
  
  ### Render UI for prediction-----------------
  output$epsilon <- renderUI(
    if (input$modelInput == 'SEIR') {
      sliderInput("epsilonInput", "Average Latency Period (Epsilon)", min = 0, max = 1,value = 1/7)}
    else {
      sliderInput("alphaInput", "Proportion of Infectiousness (Alpha)", min = 0, max = 1,value = 0.25)
    })
  
  output$lambda <- renderUI(
    if (input$modelInput == 'SEIR') {
      sliderInput("lambdaInput", "Contact Rate (Lambda)", min = 0, max = 1, value = 3/14)}
    else {
      sliderInput("betaInput", "Transmission Rate (Beta)", min = 0, max = 1,value = 0.2)
    })
  
  output$gamma <- renderUI(
    if (input$modelInput == 'SEIR') {
      sliderInput("gammaInput", "Average Infectious Period (Gamma)", min = 0, max = 1, value = 1/14)}
    else {
      sliderInput("dLInput", "Average Incubation Period (dL)", min = 0, max = 15,value = 6.4)
    })
  
  output$mu <- renderUI(
    if (input$modelInput == 'SEIR') {
      sliderInput("muInput", "Nature Birth/Death Rate (Mu)", min = 0, max = 0.0001, value = 1/365/76)}
    else {
      sliderInput("dIInput", "Average Duration of Infection (dI)", min = 0, max = 15,value = 7)
    })
  ###

  # Model reactive
  
  ### Reactive parameters for prediction----------------
  params_reactive <- reactive({
    if (input$modelInput == 'SEIR') { # for SEIR model
    req(input$epsilonInput)
    req(input$lambdaInput)
    req(input$gammaInput)
    req(input$muInput)
    list( m = 1,
      epsilon=input$epsilonInput, 
         lambda=input$lambdaInput, 
         gamma=input$gammaInput, 
         mu=input$muInput)
    } else { # for SEIcIscR model
      req(input$alphaInput)
      req(input$betaInput)
      req(input$dLInput)
      req(input$dIInput)
      C = contacts_china$all # contacts of age group j made by age group i
      list( m = 16,
            alpha=input$alphaInput, 
            beta=input$betaInput, 
            dL=input$dLInput, 
            dI=input$dIInput,
            C = contacts_china$all,
            kappa = 1-exp(-1/input$dLInput),
            gamma = 1-exp(-1/input$dIInput),
            rho = c(rep(0.4, 4), rep(0.8,12)))
    }
    })
  
  observe({
    if (input$modelInput == 'SEIR'){
      if (!is.null(input$dataset)) { # update parameters based on estimation of data
        file_path <- input$dataset
        data_toy <- read.csv(file_path$datapath)
        params <- params_reactive()
        opt <- optimize(f=mse, interval=c(0,2), 
                        data=data_toy, func=SEIR, times=times_reactive_stand(), 
                        initial_state=initial_state, params=params[-3])
        lambda_hat <- opt$minimum
        updateSliderInput(session,'lambdaInput',min = 0, max = 1,value = lambda_hat)
      } else{
        NULL
      }
    }else{
        NULL
      }
  })
  
  ### Tab 2 sidebar ui--------------------
  output$epsilon2 <- renderUI(
    if (input$modelInput == 'SEIcIscR') { # input for SEIcR in tab2
      splitLayout(sliderInput("alphaInput_2.1", "Alpha", min = 0, max = 1,value = input$alphaInput),
                  sliderInput("alphaInput_2.2", "Alpha", min = 0, max = 1,value = 1/4))
      }
    else { # for SEIR in tab2
      splitLayout(sliderInput("epsilonInput_2.1", "Epsilon", min = 0, max = 1,value = input$epsilonInput),
                  sliderInput("epsilonInput_2.2", "Epsilon", min = 0, max = 1,value = 1/4))
    })
  

  output$lambda2 <- renderUI(
    if (input$modelInput == 'SEIR') { #for SEIR in tab2
      splitLayout(sliderInput("lambdaInput_2.1", "Lambda", min = 0, max = 1,value = input$lambdaInput),
                  sliderInput("lambdaInput_2.2", "Lambda", min = 0, max = 1,value = 1/2))
    }
    else {
      splitLayout(sliderInput("betaInput_2.1", "Beta", min = 0, max = 1,value = input$betaInput),
                  sliderInput("betaInput_2.2", "Beta", min = 0, max = 1,value = 1/2))
    })
  
  output$gamma2 <- renderUI(
    if (input$modelInput == 'SEIR') {
      splitLayout(sliderInput("gammaInput_2.1", "Gamma", min = 0, max = 1,value = input$gammaInput),
                  sliderInput("gammaInput_2.2", "Gamma", min = 0, max = 1,value = 1/7))
    }
    else {
      splitLayout(sliderInput("dLInput_2.1", "dL", min = 0, max = 15,value = input$dLInput),
                  sliderInput("dLInput_2.2", "dL", min = 0, max = 15,value = 5))
    })
  
  output$mu2 <- renderUI(
    if (input$modelInput == 'SEIR') {
      splitLayout(sliderInput("muInput_2.1", "Mu", min = 0, max = 0.0001,value = input$muInput),
                  sliderInput("muInput_2.2", "Mu", min = 0, max = 0.0001,value = 0.00002))
    }
    else {
      splitLayout(sliderInput("dIInput_2.1", "dI", min = 0, max = 15,value = input$dIInput),
                  sliderInput("dIInput_2.2", "dI", min = 0, max = 15,value = 8))
    })
  ### Reactive parameters for phase plane---------------
  params_reactive_2.1 <- reactive({ # secinrao 1
    if (input$modelInput == 'SEIR') { # params for SEIR
      req(input$epsilonInput_2.1)
      req(input$lambdaInput_2.1)
      req(input$gammaInput_2.1)
      req(input$muInput_2.1)
      list( m = 1,
            epsilon=input$epsilonInput_2.1, 
            lambda=input$lambdaInput_2.1, 
            gamma=input$gammaInput_2.1, 
            mu=input$muInput_2.1)
  } else { #params for SEIcIscR
    req(input$alphaInput_2.1)
    req(input$betaInput_2.1)
    req(input$dLInput_2.1)
    req(input$dIInput_2.1)
    C = contacts_china$all # contacts of age group j made by age group i
    list( m = 16,
          alpha=input$alphaInput_2.1, 
          beta=input$betaInput_2.1, 
          dL=input$dLInput_2.1, 
          dI=input$dIInput_2.1,
          C = contacts_china$all,
          kappa = 1-exp(-1/input$dLInput_2.1),
          gamma = 1-exp(-1/input$dIInput_2.1),
          rho = c(rep(0.4, 4), rep(0.8,12)))
  } })
  
  params_reactive_2.2 <- reactive({ # secinrao 2
    if (input$modelInput == 'SEIR') {
      req(input$epsilonInput_2.2)
      req(input$lambdaInput_2.2)
      req(input$gammaInput_2.2)
      req(input$muInput_2.2)
      list( m = 1,
            epsilon=input$epsilonInput_2.2, 
            lambda=input$lambdaInput_2.2, 
            gamma=input$gammaInput_2.2, 
            mu=input$muInput_2.2)
  } else {
    req(input$alphaInput_2.2)
    req(input$betaInput_2.2)
    req(input$dLInput_2.2)
    req(input$dIInput_2.2)
    C = contacts_china$all # contacts of age group j made by age group i
    list( m = 16,
          alpha=input$alphaInput_2.2, 
          beta=input$betaInput_2.2, 
          dL=input$dLInput_2.2, 
          dI=input$dIInput_2.2,
          C = contacts_china$all,
          kappa = 1-exp(-1/input$dLInput_2.2),
          gamma = 1-exp(-1/input$dIInput_2.2),
          rho = c(rep(0.4, 4), rep(0.8,12)))
    
  } })
  
  
  # update parameters default------------
  observe({ if (input$exampleInput == 'Example2'){
    updateSliderInput(session,'epsilonInput',min = 0, max = 1,value = 1/7)
    updateSliderInput(session,'lambdaInput',min = 0, max = 1,value = 0.999)
    updateSliderInput(session,'gammaInput',min = 0, max = 1,value = 1/14)
    updateSliderInput(session,'muInput',min = 0, max = 0.0001,value = 1/365/76)
    
    
  } else if (input$exampleInput == 'Example1'){
    updateSliderInput(session,'epsilonInput',min = 0, max = 1,value = 1/7)
    updateSliderInput(session,'lambdaInput',min = 0, max = 1,value = 1/14)
    updateSliderInput(session,'gammaInput',min = 0, max = 1,value = 1/14)
    updateSliderInput(session,'muInput',min = 0, max = 0.0001,value = 1/365/76)
    
  } else if (input$exampleInput == 'Example3'){
    updateSliderInput(session,'epsilonInput',min = 0, max = 1,value = 1/7)
    updateSliderInput(session,'lambdaInput',min = 0, max = 1,value = 0.001)
    updateSliderInput(session,'gammaInput',min = 0, max = 1,value = 1/14)
    updateSliderInput(session,'muInput',min = 0, max = 0.0001,value = 1/365/76)
    
  } else {
    updateSliderInput(session,'epsilonInput',min = 0, max = 1,value = 1/7)
    updateSliderInput(session,'lambdaInput',min = 0, max = 1,value = 3/14)
    updateSliderInput(session,'gammaInput',min = 0, max = 1,value = 1/14)
    updateSliderInput(session,'muInput',min = 0, max = 0.0001,value = 1/365/76)
  }
  })

  ### Reactive result for prediction--------------
  result_reactive <- reactive({
    if (input$modelInput == 'SEIR') {
      ode(initial_state, times_reactive(), SEIR, params_reactive())
    } else if (input$modelInput == 'SEIcIscR') {
      ode(initial_state_age, times_reactive(), SEIR_age, params_reactive())
    } 
    })
  
  
  result_reactive_stand <- reactive({ # for non interactable plot
    if (input$modelInput == 'SEIR') {
      ode(initial_state, times_reactive_stand(), SEIR, params_reactive())
    } else if (input$modelInput == 'SEIcIscR') {
      ode(initial_state, times_reactive_stand(), SEIR, params_reactive())
    } else {
      ode(initial_state, times_reactive_stand(), SEIR, params_reactive())
    }
  })
  
  ### Reactive result for phase plane--------------
  result_reactive_2.1 <- reactive({ # # secinrao 1
    if (input$modelInput == 'SEIR') {
      ode(initial_state, times_reactive_stand(), SEIR, params_reactive_2.1())
    } else if (input$modelInput == 'SEIcIscR') {
      ode(initial_state_age, times_reactive(), SEIR_age, params_reactive_2.1())
    }
  })

  result_reactive_2.2 <- reactive({# secinrao 2
    if (input$modelInput == 'SEIR') {
      ode(initial_state, times_reactive_stand(), SEIR, params_reactive_2.2())
    } else if (input$modelInput == 'SEIcIscR') {
      ode(initial_state_age, times_reactive(), SEIR_age, params_reactive_2.2())
    }
  })
  
  ### phase plane plot------------
  # phase plot 1
  output$phaseplot <- renderPlot({
    if (input$modelInput == 'SEIR') {
      result1 <- as.data.frame(result_reactive_2.1())
      result1$group <- 1
      result2 <- as.data.frame(result_reactive_2.2())
      result2$group <- 2
      result <- rbind(result1, result2)
      ggplot(result, aes(x=S*100, y=R*100, color = as.factor(group))) + 
        geom_path(linewidth = 0.7) + 
        scale_color_discrete(labels = c('Scenario 1', 'Scenario 2')) +
        theme_classic2() + 
        labs(x = 'S (%)', y = 'R (%)')+
        theme(legend.title = element_blank(),
              legend.position = 'bottom',
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              axis.text = element_text(size = 10, face = 'bold'),
              axis.title = element_text(size = 12, face = 'bold'))
      } else {
        m = 16
        result1 <- result_reactive_2.1()
        if(m==1){
          result1_agg <- result1[,-1]
        } else {
          result1_agg <- sapply(0:4, function(x){
            apply(result1[,x*m+seq(1,m)+1],1,sum)
          })
        }
        colnames(result1_agg) <- c("S","E","Ic","Isc","R")
        result1_agg <- as.data.frame(result1_agg)
        result1_agg$group <- 1

        result2 <- result_reactive_2.2()
        if(m==1){
          result2_agg <- result2[,-1]
        } else {
          result2_agg <- sapply(0:4, function(x){
            apply(result2[,x*m+seq(1,m)+1],1,sum)
          })
        }
        colnames(result2_agg) <- c("S","E","Ic","Isc","R")
        result2_agg <- as.data.frame(result2_agg)
        result2_agg$group <- 2
        
        result <- rbind(result1_agg, result2_agg)
        
        ggplot(result, aes(x=S*100, y=Ic*100, color = as.factor(group))) + 
          geom_path(linewidth = 0.7) + 
          scale_color_discrete(labels = c('Scenario 1', 'Scenario 2')) +
          theme_classic2() + 
          labs(x = 'S (%)', y = 'Ic (%)')+
          theme(legend.title = element_blank(),
                legend.position = 'bottom',
                axis.line = element_line(linewidth = 0.5),
                axis.ticks = element_line(linewidth = 0.5),
                axis.text = element_text(size = 10, face = 'bold'),
                axis.title = element_text(size = 12, face = 'bold'))
      }
    })
  
  # phase plot 2
  output$phaseplotSR <- renderPlot({
    if (input$modelInput == 'SEIR'){
      result1 <- as.data.frame(result_reactive_2.1())
      result1$group <- 1
      result2 <- as.data.frame(result_reactive_2.2())
      result2$group <- 2
      result <- rbind(result1, result2)
      p <- ggplot(result, aes(x=E, y=I, color = as.factor(group))) + 
        geom_path(linewidth = 0.7) + 
        scale_color_discrete(labels = c('Scenario 1', 'Scenario 2')) +
        theme_classic2() + 
        labs(x = 'E (%)', y = 'I (%)')+
        theme(legend.title = element_blank(),
              legend.position = 'bottom',
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              axis.text = element_text(size = 10, face = 'bold'),
              axis.title = element_text(size = 12, face = 'bold'))
    } else {
      m = 16
      result1 <- result_reactive_2.1()
      result1_agg <- sapply(0:4, function(x){
          apply(result1[,x*m+seq(1,m)+1],1,sum)
        })
      colnames(result1_agg) <- c("S","E","Ic","Isc","R")
      result1_agg <- as.data.frame(result1_agg)
      result1_agg$group <- 1
      print(result1)
      
      result2 <- result_reactive_2.2()
      result2_agg <- sapply(0:4, function(x){
          apply(result2[,x*m+seq(1,m)+1],1,sum)
        })
      colnames(result2_agg) <- c("S","E","Ic","Isc","R")
      result2_agg <- as.data.frame(result2_agg)
      result2_agg$group <- 2
      
      result <- rbind(result1_agg, result2_agg)
      print(result)
      
      p <- ggplot(result, aes(x=E*100, y=R*100, color = as.factor(group))) + 
        geom_path(linewidth = 0.7) + 
        scale_color_discrete(labels = c('Scenario 1', 'Scenario 2')) +
        theme_classic2() + 
        labs(x = 'E (%)', y = 'R (%)')+
        theme(legend.title = element_blank(),
              legend.position = 'bottom',
              axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_line(linewidth = 0.5),
              axis.text = element_text(size = 10, face = 'bold'),
              axis.title = element_text(size = 12, face = 'bold'))
      
    }
    p
      
      })
  
  ### predition plot----------
  #Trajactory plot for SEIR
  output$trajplot <- renderPlotly({
    req(input$muInput)
    if (input$modelInput == 'SEIR'){
    result <- result_reactive()
    result <- pivot_longer(as.data.frame(result), cols = c(colnames(result)[-1]))
    result$name <- factor(result$name, levels = c('S', 'E', 'I', 'R'))
    result$time_chunk <- floor(result$time/10) # use time chunk to reduce the frame in the plotly
    result <- accumulate_by(result, ~ time_chunk)
    params <- params_reactive()
    r0 <- tryCatch((params$lambda * params$epsilon)/((params$epsilon + params$mu)*(params$gamma + params$mu)),
                   finally = 0)
    if (r0 > 1){
      ss <- 1/r0
      ii <- (params$mu*(1-ss))/(params$lambda * ss)
      ee <- ((params$gamma + params$mu)/params$epsilon)*ii
      rr <- 1-ss-ii-ee
    } else {
      ss <- 1
      ii <- 0
      ee <- 0
      rr <- 0
    }
      p <- ggplot(result, aes(frame = frame, x=time, y=value, color = name)) + 
        geom_hline(aes(yintercept = ss), 
                   color = pal_lancet(alpha = 0.8)(4)[1],
                   linetype = 'dashed') + 
        geom_hline(aes(yintercept = ee), 
                   color = pal_lancet(alpha = 0.8)(4)[2],
                   linetype = 'dashed') + 
        geom_hline(aes(yintercept = ii), 
                   color = pal_lancet(alpha = 0.8)(4)[3],
                   linetype = 'dashed') + 
        geom_hline(aes(yintercept = rr), 
                   color = pal_lancet(alpha = 0.8)(4)[4],
                   linetype = 'dashed') + 
        scale_color_lancet(name = '') + 
        labs(x = 'Time (days)', y = 'Proportion (%)') +
        scale_y_continuous(labels = c(0, 25,50,75,100),
                           breaks = c(0,0.25,0.5,0.75,1)) +
        theme_classic2() + 
        theme(axis.line = element_line(linewidth = 1.5),
              axis.ticks = element_line(linewidth = 1.5)) +
        geom_line(linewidth = 1)
    fig <- ggplotly(p)
    
    fig <- fig %>% animation_opts(
      frame = 100, 
      transition = 0, 
      redraw = FALSE
    )
    fig <- fig %>% layout(legend = list(orientation = 'h'))
    fig
    } else {
      div()      #fig <- fig %>% layout(legend = list(orientation = 'h'))
    }
  })
  outputOptions(output, "trajplot", suspendWhenHidden = FALSE)
  
  # trajectory plot for SEIcIscR
  output$ageplot <- renderPlotly({
    if (input$modelInput == 'SEIcIscR'){
      result <- result_reactive()
      s <- as.data.frame(result[,c(1,2:17)])
      s <- pivot_longer(s, cols = c(colnames(s)[-1]))
      s$group <- 'S'
      e <- as.data.frame(result[,c(1,18:33)])
      e <- pivot_longer(e, cols = c(colnames(e)[-1]))
      e$name <- as.numeric(e$name) - 16
      e$group <- 'E'
      ic <- as.data.frame(result[, c(1,34:49)])
      ic <- pivot_longer(ic, cols = c(colnames(ic)[-1]))
      ic$name <- as.numeric(ic$name) - 16 - 16
      ic$group <- 'Ic'
      isc <- as.data.frame(result[,c(1,50:65)])
      isc <- pivot_longer(isc, cols = c(colnames(isc)[-1]))
      isc$name <- as.numeric(isc$name) - 16 - 16 -16
      isc$group <- 'Isc'
      r <- as.data.frame(result[,c(1,66:81)])
      r <- pivot_longer(r, cols = c(colnames(r)[-1]))
      r$name <- as.numeric(r$name) - 16 - 16 - 16 - 16
      r$group <- 'R'
      
      result_4_plot <- rbind(s,e,ic,isc,r)
      print(result_4_plot)
      result_4_plot$time_chunk <- floor(result_4_plot$time/100) # use time chunk to reduce the frame in the plotly
      result_4_plot <- accumulate_by(result_4_plot, ~ time_chunk)
      result_4_plot$name <- factor(result_4_plot$name, levels = c(1:16))
      
      color <- createPalette(36,  c("#ff0000", "#00ff00", "#0000ff"))
      
      p <- ggplot(result_4_plot, aes(frame = frame, x=time, y=value*100, color = name)) +
        facet_wrap(~ group, ncol = 1, scales = 'free_y') +
        geom_line() +
        labs(x = 'Time (Days)', y = 'Proportion (%)') + 
        scale_color_manual(name = 'Age Group',
                           values = as.character(color))+
        theme_classic2()
      
      fig <- ggplotly(p, height=750)
      
      fig <- fig %>% animation_opts(
        frame = 100, 
        transition = 0, 
        redraw = FALSE
      )
      
      # Get the names of the legend entries
      df <- data.frame(id = seq_along(fig$x$data), legend_entries = unlist(lapply(fig$x$data, `[[`, "name")))
      # Extract the group identifier
      df$legend_group <- gsub("^\\((.*?),\\d+\\)", "\\1", df$legend_entries)
      # Add an indicator for the first entry per group
      df$is_first <- !duplicated(df$legend_group)
      
      for (i in df$id) {
        # Is the layer the first entry of the group?
        is_first <- df$is_first[[i]]
        # Assign the group identifier to the name and legendgroup arguments
        fig$x$data[[i]]$name <- df$legend_group[[i]]
        fig$x$data[[i]]$legendgroup <- fig$x$data[[i]]$name
        # Show the legend only for the first layer of the group 
        if (!is_first) fig$x$data[[i]]$showlegend <- FALSE
      }
    fig
    } else {
      div()
    }
  })
  
  outputOptions(output, "ageplot", suspendWhenHidden = FALSE)
  
  # non interactable trajactory table with augmented time
  output$standplot <- renderPlot({
    if (input$modelInput == 'SEIR'){
    result <- result_reactive_stand()
    result <- pivot_longer(as.data.frame(result), cols = c(colnames(result)[-1]))
    result$name <- factor(result$name, levels = c('S', 'E', 'I', 'R'))
    params <- params_reactive()
    r0 <- tryCatch((params$lambda * params$epsilon)/((params$epsilon + params$mu)*(params$gamma + params$mu)),
                   finally = 0)
    if (r0 > 1){
      ss <- 1/r0
      ii <- (params$mu*(1-ss))/(params$lambda * ss)
      ee <- ((params$gamma + params$mu)/params$epsilon)*ii
      rr <- 1-ss-ii-ee
    } else {
      ss <- 1
      ii <- 0
      ee <- 0
      rr <- 0
    }
    p <- ggplot(result, aes(x=time, y=value, color = name)) + 
      geom_hline(aes(yintercept = ss), 
                 color = pal_lancet(alpha = 0.8)(4)[1],
                 linetype = 'dashed') + 
      geom_hline(aes(yintercept = ee), 
                 color = pal_lancet(alpha = 0.8)(4)[2],
                 linetype = 'dashed') + 
      geom_hline(aes(yintercept = ii), 
                 color = pal_lancet(alpha = 0.8)(4)[3],
                 linetype = 'dashed') + 
      geom_hline(aes(yintercept = rr), 
                 color = pal_lancet(alpha = 0.8)(4)[4],
                 linetype = 'dashed') + 
      labs(x = 'Time (days)', y = 'Proportion (%)') +
      scale_y_continuous(labels = c(0, 25,50,75,100),
                         breaks = c(0,0.25,0.5,0.75,1)) +
      scale_color_lancet() + 
      theme_classic2() + 
      theme(axis.line = element_line(linewidth = 1),
            axis.ticks = element_line(linewidth = 1),
            legend.title = element_blank(),
            legend.position = 'bottom',
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, family = 'sans')) +
      geom_line(linewidth = 1.2)
    } else {
      model <- result_reactive()
      m <- 16
      if(m==1){
        model_agg <- model[,-1]
      } else {
        model_agg <- sapply(0:4, function(x){
          apply(model[,x*m+seq(1,m)+1],1,sum)
        })
      }
      colnames(model_agg) <- c("S","E","Ic","Isc","R")
      p <- matplot(model_agg, type="l", lty=1, main="SEIR model", xlab="Time")
    }
    p
  })
  
  # equilibrium table
  output$eqtable <- renderText({
    params <- params_reactive_2.1()
    r0_1 <- tryCatch((params$lambda * params$epsilon)/((params$epsilon + params$mu)*(params$gamma + params$mu)),
                   finally = 0)
    if (r0_1 > 1){
      ss_1 <- 1/r0_1
      ii_1 <- (params$mu*(1-ss_1))/(params$lambda * ss_1)
      ee_1 <- ((params$gamma + params$mu)/params$epsilon)*ii_1
      rr_1 <- 1-ss_1-ii_1-ee_1
    } else {
      ss_1 <- 1
      ii_1 <- 0
      ee_1 <- 0
      rr_1 <- 0
    }
    
    params <- params_reactive_2.2()
    r0_2 <- tryCatch((params$lambda * params$epsilon)/((params$epsilon + params$mu)*(params$gamma + params$mu)),
                     finally = 0)
    if (r0_2 > 1){
      ss_2 <- 1/r0_2
      ii_2 <- (params$mu*(1-ss_2))/(params$lambda * ss_2)
      ee_2 <- ((params$gamma + params$mu)/params$epsilon)*ii_2
      rr_2 <- 1-ss_2-ii_2-ee_2
    } else {
      ss_2 <- 1
      ii_2 <- 0
      ee_2 <- 0
      rr_2 <- 0
    }
    
    eq <- data.frame(S = c(ss_1, ss_2), E = c(ee_1, ee_2), I = c(ii_1, ii_2),
               R = c(rr_1, rr_2), R0 = c(r0_1, r0_2), 
               row.names = c('Scenario 1', 'Scenario 2'))
    kable(eq, 
          digits = 3,
          booktabs = TRUE) %>% 
      kable_styling(latex_options = c("HOLD_position", "striped"),
                    font_size=15)
    
  })
  
  

}

  

#main ----------
shinyApp(ui = ui, server = server)


