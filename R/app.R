# ==============
# Robert Smith & Paul Schneider
# University of Sheffield
# contact: rasmith3@sheffield.ac.uk
# ==============

## app.R ##
library(shiny)
library(shinydashboard)

source("wrapper.R")

ui <- fluidPage(
  
  titlePanel("Sick Sicker Model in Shiny"),
  
  # SIDEBAR
  sidebarLayout(
    
    sidebarPanel(
    
    numericInput(inputId = "SI_c_Trt",
                 label = "Treatment Cost",
                 value = 200,
                 min = 0,
                 max = 400),
    
    numericInput(inputId = "SI_n_sim",
                 label = "PSA runs",
                 value = 1000,
                 min = 0,
                 max = 400),
    
    sliderInput(inputId = "SI_n_age_init",
                label = "initial age",
                min = 10,
                max = 80,
                value = 25),
    

    actionButton("run_model","Run / update model")
                       
                ),  # close sidebarPanel
    
  mainPanel(
    
    h3("Results Table"),
    
    tableOutput("icer_table"),
    
    h3("Cost-effectiveness Plane"),
    
    plotOutput("CE_plane")
    
            ) # close mainpanel    
      
      ) # close sidebarlayout
                     
  ) # close UI fluidpage


server <- function(input, output){
  
 observeEvent(input$run_model,
              ignoreNULL = F, {
    
    
  df_model_res = f_wrapper(shiny_c_Trt = input$SI_c_Trt,
                           n_age_init = input$SI_n_age_init,
                           n_sim = input$SI_n_sim)
    
  
  #  }) # actiobbutton end
    # ICER TABLE
    output$icer_table <- renderTable({
      
      df_res_table <- data.frame(
        
        Option =  c("Treatment","No Treatment"),
        
        QALYs  =  c(mean(df_model_res$QALY_Trt),mean(df_model_res$QALY_NoTrt)),
        
        Costs  =  c(mean(df_model_res$Cost_Trt),mean(df_model_res$Cost_NoTrt)),
        
        Inc.QALYs = c(mean(df_model_res$QALY_Trt) - mean(df_model_res$QALY_NoTrt),NA),
        
        Inc.Costs = c(mean(df_model_res$Cost_Trt) - mean(df_model_res$Cost_NoTrt),NA),
        
        ICER = c(mean(df_model_res$ICER),NA)
      )
      
      df_res_table[,2:6] <- round(df_res_table[,2:6],digits = 2) 
      
      df_res_table

      }) # table plot end.
    
    
   #  CE PLANE
    output$CE_plane <- renderPlot({
      
      #model_res = model_res() # calling reactive function
      df_model_res$inc_C <- df_model_res$Cost_Trt - df_model_res$Cost_NoTrt
      df_model_res$inc_Q <- df_model_res$QALY_Trt - df_model_res$QALY_NoTrt
      
      # create cost effectiveness plane plot
      plot(x = df_model_res$inc_Q,
           y = df_model_res$inc_C,
           xlab = "Incremental QALYs",
           ylab = "Incremental Costs",
           xlim = c(min(df_model_res$inc_Q,df_model_res$inc_Q*-1),
                    max(df_model_res$inc_Q,df_model_res$inc_Q*-1)),
           ylim = c(min(df_model_res$inc_C,df_model_res$inc_C*-1),
                    max(df_model_res$inc_C,df_model_res$inc_C*-1)),
           abline(h = 0,v=0)
          ) # plot end 
    }) # renderplot end
    
 }) # Observe Event End
    
    
  } # Server end
  
  



## ----- run app------

shinyApp(ui, server)