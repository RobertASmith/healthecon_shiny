# ==============
# 
# SCRIPT 1
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
#                   Create User Interface
#================================================================

ui <- fluidPage(    # create user interface using fluid-page function
  
  "hello"
  
) # close UI fluidpage


#================================================================
#                     Create Server Function
#================================================================

server <- function(input, output){   
  
  # server function

  
} # Server end





## ----- run app------

shinyApp(ui, server)