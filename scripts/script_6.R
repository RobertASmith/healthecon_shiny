# ==============
# 
# SCRIPT 6
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

#================================================================
#source the function from an R script:
#================================================================

source("./scripts/fun_script_6.R")

#================================================================
#                   Create User Interface
#================================================================

ui <- fluidPage(    # create user interface using fluidpage function
  
  titlePanel("More complex model"),   # title of app
  
  # SIDEBAR
  sidebarLayout(    # indicates layout is going to be a sidebar-layout
    
    sidebarPanel( # open sidebar panel
      
      numericInput(inputId = "SI_X",      # id of input, used in server
                   label = "Number X",  # label next to numeric input
                   value = 200,               # initial value
                   min = 0,                   # minimum value allowed
                   max = 400),                # maximum value allowed
      
      numericInput(inputId = "SI_Y",      # id of input, used in server
                   label = "Number Y",        # label next to numeric input
                   value = 1000,              # initial value
                   min = 0,                   # minimum value allowed
                   max = 400),                # maximum value allowed
      
      sliderInput(inputId = "SI_Z",  # id of input, used in server
                  label = "Number Z",      # label next to numeric input
                  value = 25,                 # initial value
                  min = 10,                   # minimum value allowed
                  max = 1000),                  # maximum value allowed
      
      
      actionButton(inputId = "run_model",     # id of action button, used in server
                   label   = "Run model")     # action button label (on button)
      
    ),  # close sidebarPanel
    
    mainPanel(                                # open main panel
      
      h3("Basic Plot"),                           # heading (Cost effectiveness plane)
      
      plotOutput(outputId = "SO_plot")       # plot Output id = CE_plane, from server
      
      
      
    ) # close main panel    
    
  ) # close sidebar layout
  
) # close UI fluid page


#================================================================
#                     Create Server Function
#================================================================

server <- function(input, output){   # server = function with two inputs
  
  
  
  observeEvent(input$run_model,       # when action button pressed ...
               ignoreNULL = F, {
                 
                df = fun_shiny(x = input$SI_X,
                               y = input$SI_Y,
                               z = input$SI_Z)

                 output$SO_plot <- renderPlot({ # render plot repeatedly updates.
                 
                   plot(x = df$x,
                        y = df$y) # plot end
                   
                 }) # renderplot end  
                 
                 
               }) # Observe Event End  
  
  
} # Server end





## ----- run app------

shinyApp(ui, server)