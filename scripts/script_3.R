# ==============
# 
# SCRIPT 3
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
  
  titlePanel("Multiply a number by 10 when button pressed"),   # title of app
  
  # SIDEBAR
  sidebarLayout(    # indicates layout is going to be a sidebar-layout
    
    sidebarPanel( # open sidebar panel
      
      numericInput(inputId = "num",           # id of input, used in server
                   label = "number",          # label next to numeric input
                   value = 200,               # initial value
                   min = 0,                   # minimum value allowed
                   max = 100000),              # maximum value allowed
                   
                   
      actionButton(inputId = "run_model",     # id of action button, used in server
                   label   = "Run model") # action button end
      
      ), # sidebar panel end
    
    mainPanel(                                # open main panel
      
      textOutput(outputId = "printvalue")     # text output                
      
    ) # close main panel    
    
  ) # close sidebar layout
  
) # close UI fluid page


#================================================================
#                     Create Server Function
#================================================================

server <- function(input, output){   # server = function with two inputs
  
  observeEvent(input$run_model,       # when action button pressed ...
               ignoreNULL = F, {
                 
     number_times_10 = input$num * 10
  
  #--- CREATE NUMBER IN SERVER ---#
  output$printvalue <- renderText({
    
    paste("number x 10 = ",number_times_10)
    
  }) # render Text end.
  
  
}) # Observe Event End  
  
  
} # Server end





## ----- run app------

shinyApp(ui, server)
