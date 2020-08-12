# ==============
# 
# SCRIPT 4
#
# Making Markov Models Shiny 
# Robert Smith & Paul Schneider
# University of Sheffield
# contact: rasmith3@sheffield.ac.uk
# ==============

# install.packages("shiny") # necessary if you don't already have the function 'shiny' installed.

# we need the function shiny installed, this loads it from the library.
library(shiny)             

#================================================================
#                          Simple Function                      
#================================================================

fun_shiny <- function(x,y,z){
  
  max_minus_mean = max(c(x,y,z)) - mean(c(x,y,z))
  
  max_minus_mean
  
}



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
                  max = 80),                  # maximum value allowed
      
      
      actionButton(inputId = "run_model",     # id of action button, used in server
                   label   = "Run model")     # action button label (on button)
      
    ),  # close sidebarPanel
    
    mainPanel(                                # open main panel
      
      h3("Results"),                    # heading (results)                
    
      textOutput(outputId = "printvalue")                    # heading (results table)                
        
      
    ) # close mainpanel    
    
  ) # close sidebarlayout
  
) # close UI fluidpage


#================================================================
#                     Create Server Function
#================================================================

server <- function(input, output){   # server = function with two inputs
  

  
  observeEvent(input$run_model,       # when action button pressed ...
               ignoreNULL = F, {
                 
                 num = fun_shiny(x = input$SI_X,
                                 y = input$SI_Y,
                                 z = input$SI_Z)
                 
                 
                 #--- CREATE NUMBER IN SERVER ---#
                 output$printvalue <- renderText({
                   
                   paste("The difference between the maximum and the mean of x, y, and z is:", round(num,1))
                   
                 }) # render Text end.
                 
                 
               }) # Observe Event End  
  
  
} # Server end





## ----- run app------

shinyApp(ui, server)
