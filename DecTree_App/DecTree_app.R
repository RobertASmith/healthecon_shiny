# ==============
# Making Decision Models Shiny 
# Robert Smith & Paul Schneider
# University of Sheffield
# contact: rasmith3@sheffield.ac.uk
# 18/03/2020
# Simple decision tree example
# ==============

## app.R ##

# install.packages("shiny") # necessary if you don't already have the function 'shiny' installed.

# we need the function shiny installed, this loads it from the library.
library(shiny)             

# source the wrapper function.
source("../DecTree_app/f_model.R")   

#================================================================
#                   Create User Interface
#================================================================

ui <- fluidPage(    # create user interface using fluidpage function
  
  titlePanel("Decision Tree Model in Shiny: Treatment A vs Treatment B"),   # title of app
  
  # SIDEBAR
  sidebarLayout(    # indicates layout is going to be a sidebar-layout
    
    sidebarPanel( # open sidebar panel
      
      selectInput(inputId = "SI_n_draws",                       # id of input, used in server
                   label = "Number of Draws from Distribution", # label above dropdown
                   choices = list(100, 500, 1000, 2000, 5000),  # choices in dropdown box
                   selected = 1000,                             # initial selection
                   multiple = F),                               # allow multiple selections?               
      
      numericInput(inputId = "SI_c_B_tr_mean",      # id of input, used in server
                   label = "Mean for log-normal distribution of costs of B",  # label next to numeric input
                   value = 1000,               # initial value
                   min = 0,                   # minimum value allowed
                   max = 2000),                # maximum value allowed
      
      numericInput(inputId = "SI_c_B_tr_sd",      # id of input, used in server
                   label = "Standard deviation for log-normal distribution of costs of B",  # label next to numeric input
                   value = 100,               # initial value
                   min = 0,                   # minimum value allowed
                   max = 200),                # maximum value allowed
      
      sliderInput(inputId = "SI_rr_mean",  # id of input, used in server
                  label = "Mean for distribution of 'relative risk' of B vs A",      # label next to numeric input
                  value = 0.95,                 # initial value
                  min = 0.5,                   # minimum value allowed
                  max = 1),                  # maximum value allowed
      
      sliderInput(inputId = "SI_rr_sd",  # id of input, used in server
                  label = "Standard deviation for distribution of rr",      # label next to numeric input
                  value = 0.01,                 # initial value
                  min = 0.0000001,              # minimum value allowed
                  max = 0.03),                  # maximum value allowed
      
      sliderInput(inputId = "SI_p_A_S1",  # id of input, used in server
                  label = "Probability of S1 with treatment A",      # label next to numeric input
                  value = 0.05,                 # initial value
                  min = 0.01,                   # minimum value allowed
                  max = 0.08),                  # maximum value allowed
      
      actionButton(inputId = "run_model",     # id of action button, used in server
                   label   = "Run model")     # action button label (on button)
      
    ),  # close sidebarPanel
    
    mainPanel(                                # open main panel
      
      h3("Results Table"),                    # heading (results table)                
      
      tableOutput(outputId = "SO_icer_table"),   # tableOutput id = icer_table, from server
      
      br(),br(),br(),
      h3("Cost-effectiveness Plane: Treatment B vs Treatment A"),         # heading (Cost effectiveness plane)
      
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
                 df_model_res = f_model(n_draws     = input$SI_n_draws,      # number of PSA iterations
                                        c_B_tr_mean = input$SI_c_B_tr_mean,  # mean costs of B
                                        c_B_tr_sd   = input$SI_c_B_tr_sd,     # max costs of B
                                        rr_mean     = input$SI_rr_mean,           # Mean Relative risk under B
                                        rr_sd       = input$SI_rr_sd,            # standard deviation of rr
                                        p_A_S1      = input$SI_p_A_S1        # prob. of S1 with treatment A
                                        )
                
                 
                 #--- CREATE COST EFFECTIVENESS PLANE ---#
                 output$SO_icer_table <- renderTable({ # this continuously updates table
                   
                   df_res_table <- data.frame( # create dataframe
                     
                     Option =  c("Treatment A", "Treatment B"), 
                     
                     QALYs  =  c(mean(df_model_res$u_A),mean(df_model_res$u_B)),
                     
                     Costs  =  c(mean(df_model_res$c_A),mean(df_model_res$c_B)),
                     
                     Inc.QALYs = c(NA, mean(df_model_res$u_incr)),
                     
                     Inc.Costs = c(NA, mean(df_model_res$c_incr)),
                     
                     ICER = c(NA, mean(df_model_res$c_incr)/mean(df_model_res$u_incr))
                   )
                   
                   # round the dataframe to two digits so looks tidier
                   df_res_table[,2:6] <- round(df_res_table[,2:6],digits = 2) 
                   
                   #print the dataframe
                   df_res_table
                   
                 }) # table plot end.
                 
                 
                 #---  CREATE COST EFFECTIVENESS PLANE ---#
                 output$SO_CE_plane <- renderPlot({ # render plot repeatedly updates.
                   
                  # baseR plot
                  plot(x=df_model_res$u_incr,
                        y=df_model_res$c_incr,
                        cex=0.8,
                        col="cadetblue4",
                        xlab = "Incremental Life Years",
                        ylab = "incremental Costs",
                        xlim = c(min(df_model_res$u_incr,df_model_res$u_incr*-1),
                                max(df_model_res$u_incr,df_model_res$u_incr*-1)),
                        ylim = c(min(df_model_res$c_incr,df_model_res$c_incr*-1),
                                max(df_model_res$c_incr,df_model_res$c_incr*-1)))
                   abline(h=0,v=0) # straight vertical and horizontal lines at zero.
                   abline(coef= c(0,25000),lty=2) # dotted line at lambda.
                    
                 }) # renderplot end
                 
               }) # Observe Event End
  
  
} # Server end

  



## ----- run app------

shinyApp(ui, server)