# ==============
# Making Markov Models Shiny 
# Robert Smith & Paul Schneider
# University of Sheffield
# contact: rasmith3@sheffield.ac.uk
# ==============

## app.R ##

# install.packages("shiny") # necessary if you don't already have the function 'shiny' installed.

# we need the function shiny installed, this loads it from the library.
library(shiny)             

# source the wrapper function.
source("./wrapper.R")   

#================================================================
#                   Create User Interface
#================================================================

ui <- fluidPage(    # create user interface using fluidpage function
  
  titlePanel("Sick Sicker Model in Shiny"),   # title of app
  
  # SIDEBAR
  sidebarLayout(    # indicates layout is going to be a sidebar-layout
    
    sidebarPanel( # open sidebar panel
      
      numericInput(inputId = "SI_c_Trt",      # id of input, used in server
                   label = "Treatment Cost",  # label next to numeric input
                   value = 200,               # initial value
                   min = 0,                   # minimum value allowed
                   max = 400),                # maximum value allowed
      
      numericInput(inputId = "SI_n_sim",      # id of input, used in server
                   label = "PSA runs",        # label next to numeric input
                   value = 1000,              # initial value
                   min = 0,                   # minimum value allowed
                   max = 400),                # maximum value allowed
      
      sliderInput(inputId = "SI_n_age_init",  # id of input, used in server
                  label = "Initial Age",      # label next to numeric input
                  value = 25,                 # initial value
                  min = 10,                   # minimum value allowed
                  max = 80),                  # maximum value allowed
      
      
      actionButton(inputId = "run_model",     # id of action button, used in server
                   label   = "Run model")     # action button label (on button)
      
    ),  # close sidebarPanel
    
    mainPanel(                                # open main panel
      
      h3("Results Table"),                    # heading (results table)                
      
      tableOutput(outputId = "SO_icer_table"),   # tableOutput id = icer_table, from server
      
      h3("Cost-effectiveness Plane"),         # heading (Cost effectiveness plane)
      
      plotOutput(outputId = "SO_CE_plane")       # plotOutput id = CE_plane, from server
      
    ) # close mainpanel    
    
  ) # close sidebarlayout
  
) # close UI fluidpage


#================================================================
#                     Create Server Function
#================================================================

server <- function(input, output){   # server = function with two inputs
  
  observeEvent(input$run_model,       # when action button pressed ...
               ignoreNULL = F, {
                 
                 # Run model wrapper function with the Shiny inputs and store as data-frame 
                 df_model_res = f_wrapper(c_Trt = input$SI_c_Trt,
                                          n_age_init = input$SI_n_age_init,
                                          n_sim = input$SI_n_sim)
                 
                 
                 #--- CREATE COST EFFECTIVENESS PLANE ---#
                 output$SO_icer_table <- renderTable({ # this continuously updates table
                   
                   df_res_table <- data.frame( # create dataframe
                     
                     Option =  c("Treatment","No Treatment"), 
                     
                     QALYs  =  c(mean(df_model_res$QALY_Trt),mean(df_model_res$QALY_NoTrt)),
                     
                     Costs  =  c(mean(df_model_res$Cost_Trt),mean(df_model_res$Cost_NoTrt)),
                     
                     Inc.QALYs = c(mean(df_model_res$QALY_Trt) - mean(df_model_res$QALY_NoTrt),NA),
                     
                     Inc.Costs = c(mean(df_model_res$Cost_Trt) - mean(df_model_res$Cost_NoTrt),NA),
                     
                     ICER = c(mean(df_model_res$ICER),NA)
                   )
                   
                   # round the dataframe to two digits so looks tidier
                   df_res_table[,2:6] <- round(df_res_table[,2:6],digits = 2) 
                   
                   #print the dataframe
                   df_res_table
                   
                 }) # table plot end.
                 
                 
                 #---  CREATE COST EFFECTIVENESS PLANE ---#
                 output$SO_CE_plane <- renderPlot({ # render plot repeatedly updates.
                   
                   # calculate incremental costs and qalys from results dataframe
                   df_model_res$inc_C <- df_model_res$Cost_Trt - df_model_res$Cost_NoTrt
                   df_model_res$inc_Q <- df_model_res$QALY_Trt - df_model_res$QALY_NoTrt
                   
                   # create cost effectiveness plane plot
                   plot(x = df_model_res$inc_Q, # x axis incremental QALYS
                        y = df_model_res$inc_C, # y axis incremental Costs
                        #label axes
                        xlab = "Incremental QALYs", 
                        ylab = "Incremental Costs", 
                        
                        # set xlimits and ylimits for plot.
                        xlim = c(min(df_model_res$inc_Q,df_model_res$inc_Q*-1),
                                 max(df_model_res$inc_Q,df_model_res$inc_Q*-1)),
                        ylim = c(min(df_model_res$inc_C,df_model_res$inc_C*-1),
                                 max(df_model_res$inc_C,df_model_res$inc_C*-1)),
                        # include y and y axis lines.
                        abline(h = 0,v=0)
                   ) # plot end 
                 }) # renderplot end
                 
               }) # Observe Event End
  
  
} # Server end

  



## ----- run app------

shinyApp(ui, server)