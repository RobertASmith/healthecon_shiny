## Sick Sicker example app files 

This is an example application which aims to be as simple a functioning example as possible. The full tutorial based on this simple application is here: 
https://github.com/RobertASmith/healthecon_shiny/tree/master/Tutorial.
For those completely new to Shiny a good place to start is here: https://shiny.rstudio.com/articles/basics.html

This folder contains two files:
- app.R
- wrapper.R

If you are unfamiliar with GitHub and wish to download them to your desktop, ensure both files are within the same folder.

### app.R
This file contains a script with the user interface and server functions. In the top right corner there is a button 'Run App'. If you click this the application should then pop up.
In order to understand this code you will need to know: 
- what a reactive function does: https://shiny.rstudio.com/articles/reactivity-overview.html
- how an action button works: https://shiny.rstudio.com/reference/shiny/latest/observeEvent.html

### wrapper.R
This file contains a script with the model wrapper function in it. The functions *f_gen_psa* and *f_MM_sicksicker* are created within this function for simplicity. In more complex examples with longer code chunks functions would generally be created externally and sourced.
In order to understand this code you will need to know: 
- What a function is: https://www.tutorialgateway.org/functions-in-r-programming/
- What matrix multiplication (%*%) is https://stat.ethz.ch/R-manual/R-patched/library/base/html/matmult.html ; https://www.dummies.com/education/math/calculus/how-to-multiply-matrices-by-each-other/





