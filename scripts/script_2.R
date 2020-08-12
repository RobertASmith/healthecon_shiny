# ==============
# 
# SCRIPT 2
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
#                   Create User Interface
#================================================================

ui <- fluidPage(    # create user interface using fluidpage function
  
  titlePanel("Multiply a number by 10"),   # title of app
  
  # SIDEBAR
  sidebarLayout(    # indicates layout is going to be a sidebar-layout
    
    sidebarPanel( # open sidebar panel
      
      numericInput(inputId = "num",      # id of input, used in server
                   label = "number",  # label next to numeric input
                   value = 200,               # initial value
                   min = 0,                   # minimum value allowed
                   max = 100000)                # maximum value allowed
    ),  # close sidebarPanel
    
    mainPanel(                                # open main panel
      
      textOutput(outputId = "printvalue")                    # heading (results table)                
      
    ) # close mainpanel    
    
  ) # close sidebarlayout
  
) # close UI fluidpage


#================================================================
#                     Create Server Function
#================================================================

server <- function(input, output){   # server = function with two inputs
  
                 
                 #--- CREATE NUMBER IN SERVER ---#
                 output$printvalue <- renderText({
                   
                   paste("number x 10 = ", input$num * 10)
                   
                 }) # render Text end.
                 
  
  
} # Server end





## ----- run app------

shinyApp(ui, server)