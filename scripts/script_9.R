# ==============
# 
# SCRIPT 9 - HTML formatting
#
# Making Markov Models Shiny 
# Robert Smith & Paul Schneider
# University of Sheffield
# contact: rasmith3@sheffield.ac.uk
# ==============

rm(list = ls())
# install.packages("shiny") # necessary if you don't already have the function 'shiny' installed.

# we need the function shiny installed, this loads it from the library.
library(shiny)  
library(ggplot2)

# source the wrapper function.
source("./App/wrapper.R") 
source("./scripts/fun_script_8.R")

#================================================================
#                   Create User Interface
#================================================================

ui <- fluidPage(    # create user interface using fluid-page function
  
  titlePanel("Sick Sicker Model in Shiny"),   # title of app
  
  # SIDEBAR
  sidebarLayout(    # indicates layout is going to be a sidebar-layout
    
    sidebarPanel( # open sidebar panel
      
      numericInput(inputId = "SI_c_Trt",      # id of input, used in server
                   label = "Treatment Cost",  # label next to numeric input
                   value = 200,               # initial value
                   min = 0,                   # minimum value allowed
                   max = 400),                # maximum value allowed
      
      HTML("<br>"), # add a line break
      
      numericInput(inputId = "SI_n_sim",      # id of input, used in server
                   label = "PSA runs",        # label next to numeric input
                   value = 1000,              # initial value
                   min = 0,                   # minimum value allowed
                   max = 400),                # maximum value allowed
      
      HTML("<br>"),  # add a line break
      
      sliderInput(inputId = "SI_n_age_init",  # id of input, used in server
                  label = "Initial Age",      # label next to numeric input
                  value = 25,                 # initial value
                  min = 10,                   # minimum value allowed
                  max = 80),                  # maximum value allowed
      
      HTML("<br>"),  # add a line break
      
      
      actionButton(inputId = "run_model",     # id of action button, used in server
                   label   = "Run model",     # action button label (on button)
                   style = "width: 150px;     
                            background-color: #DC6E64;  
                            border-radius: 40px;  
                            color: white;         
                            font-weight: bold;    
                            font-size:100%;       
                            text-align:center")  # customised button    
      
    ),  # close sidebar panel
    
    mainPanel(                                # open main panel
      
      h3("Results Table"),                    # heading (results table)
      
      HTML("<br>"),
      
      tableOutput(outputId = "SO_icer_table"),   # tableOutput id = icer_table, from server
      
      h3("Cost-effectiveness Plane"),         # heading (Cost effectiveness plane)
      
      plotOutput(outputId = "SO_CE_plane"),       # plotOutput id = CE_plane, from server
      
      downloadButton('cep')
      
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
                   
                   df_res_table <- data.frame( # create data-frame
                     
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
                   
                   # use function ce_plot frrom fun_script_8 file to create plot.
                   plot <- ce_plot(results = df_model_res)
                   
                   # save cost-effectiveness plane for download
                   ceP_download <<-  reactive({plot})
                   
                   # output plot from function.
                   plot
                   
                 }) # render plot end
                 
                 
                 # cost effectiveness plane fig. download ----
                 output$cep = downloadHandler(
                   filename = 'ce_plane.png',    # select file name
                   content = function(file) {
                     device <- function(..., width, height) {
                       grDevices::png(..., 
                                      width = width, 
                                      height = height,
                                      res = 300, 
                                      units = "in")
                     }
                     ggsave(file, 
                            plot = ceP_download(), # need to remember to have "()" after the ceP_download we created above!
                            device = device)
                   })
                 
               }) # Observe Event End
  
  
} # Server end





## ----- run app------

shinyApp(ui, server)

