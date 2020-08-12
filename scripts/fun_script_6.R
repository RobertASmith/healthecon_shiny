#================================================================
#                          Plot Function - for script 6                      
#================================================================

fun_shiny <- function(x,y,z){
  
  df = data.frame(x =  rnorm(n = z, mean = x, sd = x/10),
                  y =  rnorm(n = z, mean = y, sd = y/5))
                
  df
  
}


