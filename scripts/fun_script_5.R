#================================================================
#                          Simple Function                      
#================================================================

fun_shiny <- function(x,y,z){
  
  max_minus_mean = max(c(x,y,z)) - mean(c(x,y,z))
  
  max_minus_mean
  
}