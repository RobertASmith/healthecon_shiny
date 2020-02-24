# Making Markov Models Shiny: A Tutorial

Robert Smith and Paul Schneider

ScHARR, University of Sheffield

21 February 2020

*Please cite this tutorial as: XXXXXXXXXXXXXXXXX.

This tutorial aims to provide a step by step guide in the use of R-Shiny for health economic modelling. It uses a previously published 4-state markov model, the Sick-Sicker model (Krijkamp et al., 2020; Alarid-Escudero et al., 2020) as a case-study, using the DARTH coding framework when creating the shiny app (Alarid-Escudero et al., 2019). We have adapted it such that it has one purpose for this tutorial, to create PSA outputs.

The diagram below shows the general model structure.

![Sick Sicker Diagram](https://github.com/RobertASmith/healthecon_shiny/blob/master/Tutorial/Sick%20Sicker%20Diagram.PNG)

In this case there are two functions witin the model, the first gen_psa creates a set of psa inputs, the second runs the markov model for a specific set of PSA inputs.

## Functions

### Creating PSA inputs.

```{r echo = T, eval = F}
f_gen_psa <- function(n_sim = 1000, SI_c_Trt){
   
  df_psa <- data.frame(
    
    # Transition probabilities (per cycle)
    p_HS1   = rbeta(n_sim, 30, 170),        # prob Healthy -> Sick
    p_S1H   = rbeta(n_sim, 60, 60) ,        # prob Sick    -> Healthy
    p_S1S2  = rbeta(n_sim, 84, 716),        # prob Sick    -> Sicker
    
    p_HD    = rbeta(n_sim, 10, 1990)      ,  # prob Healthy -> Dead
    hr_S1   = rlnorm(n_sim, log(3),  0.01),  # rate ratio death S1 vs healthy
    hr_S2   = rlnorm(n_sim, log(10), 0.02),  # rate ratio death S2 vs healthy 
    
    # Cost vectors with length n_sim
    c_H   = rgamma(n_sim, shape = 100, scale = 20)    , # cost p/cycle in state H
    c_S1  = rgamma(n_sim, shape = 177.8, scale = 22.5), # cost p/cycle in state S1
    c_S2  = rgamma(n_sim, shape = 225, scale = 66.7)  , # cost p/cycle in state S2
    c_D   = 0                                         , # cost p/cycle in state D
    c_Trt = SI_c_Trt,                                   # cost p/cycle of treatment
    
    # Utility vectors with length n_sim 
    u_H   = rtruncnorm(n_sim, mean =    1, sd = 0.01, b = 1), # utility when healthy
    u_S1  = rtruncnorm(n_sim, mean = 0.75, sd = 0.02, b = 1), # utility when sick
    u_S2  = rtruncnorm(n_sim, mean = 0.50, sd = 0.03, b = 1), # utility when sicker
    u_D   = 0                                               , # utility when dead
    u_Trt = rtruncnorm(n_sim, mean = 0.95, sd = 0.02, b = 1)  # utility when being treated
  )
  
  return(df_psa)
}

```

### Running the model for a specific set of PSA inputs

```{r echo = T, eval = F}

f_MM_sicksicker <- function(params) {
  with(as.list(params), {
    
    # compute internal paramters as a function of external parameter
    r_HD    = - log(1 - p_HD) # rate of death in healthy
    r_S1D   = hr_S1 * r_HD 	  # rate of death in sick
    r_S2D   = hr_S2 * r_HD  	# rate of death in sicker
    p_S1D   = 1 - exp(-r_S1D) # probability to die in sick
    p_S2D   = 1 - exp(-r_S2D) # probability to die in sicker
    
    # calculate discount weight for each cycle based on discount rate d_r
    v_dwe <- v_dwc <- 1 / (1 + d_r) ^ (0:n_t) 
    
    # create transition probability matrix for NO treatment
    m_P <- matrix(0,
                  nrow = n_states, ncol = n_states,
                  dimnames = list(v_n, v_n))
    # fill in the transition probability array
    ### From Healthy
    m_P["H", "H"]  <- 1 - (p_HS1 + p_HD)
    m_P["H", "S1"] <- p_HS1
    m_P["H", "D"]  <- p_HD
    ### From Sick
    m_P["S1", "H"]  <- p_S1H
    m_P["S1", "S1"] <- 1 - (p_S1H + p_S1S2 + p_S1D)
    m_P["S1", "S2"] <- p_S1S2
    m_P["S1", "D"]  <- p_S1D
    ### From Sicker
    m_P["S2", "S2"] <- 1 - p_S2D
    m_P["S2", "D"]  <- p_S2D
    ### From Dead
    m_P["D", "D"] <- 1
    
    # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    m_TR <- matrix(NA, nrow = n_t + 1 , ncol = n_states, 
                   dimnames = list(0:n_t, v_n))     
    
    m_TR[1, ] <- c(1, 0, 0, 0)          # initialize Markov trace
    
    ############# PROCESS ###########################################
      
    
    for (t in 1:n_t){ # throughout the number of cycles
      # estimate the Markov trace for cycle the next cycle (t + 1)
      m_TR[t + 1, ] <- m_TR[t, ] %*% m_P           
    }
    
    ############ OUTPUT  ###########################################
    
    # create vectors of utility and costs for each state
    v_u_trt    <- c(u_H, u_Trt, u_S2, u_D)
    v_u_no_trt <- c(u_H, u_S1, u_S2, u_D)
    
    v_c_trt    <- c(c_H, c_S1 + c_Trt, c_S2 + c_Trt, c_D)
    v_c_no_trt <- c(c_H, c_S1, c_S2, c_D)
    
    # estimate mean QALys and costs
    v_E_no_trt <- m_TR %*% v_u_no_trt
    v_E_trt    <- m_TR %*% v_u_trt
    
    v_C_no_trt <- m_TR %*% v_c_no_trt
    v_C_trt    <- m_TR %*% v_c_trt
    
    ### discount costs and QALYs
    te_no_trt <- t(v_E_no_trt) %*% v_dwe  # 1x31 %*% 31x1 -> 1x1
    te_trt    <- t(v_E_trt) %*% v_dwe
    
    tc_no_trt <- t(v_C_no_trt) %*% v_dwc
    tc_trt    <- t(v_C_trt)    %*% v_dwc
    
    results <- c("Cost_NoTrt" = tc_no_trt, 
                 "Cost_Trt"   = tc_trt, 
                 "QALY_NoTrt" = te_no_trt, 
                 "QALY_Trt"   = te_trt,
                 "ICER"       = (tc_trt - tc_no_trt)/(te_trt - te_no_trt))
    
    return(results)
  }
  )
}
```

## Creating a Wrapper
When using a web application it is likely that the user will want to be able to change parameter inputs and re-run the model. In order to make this simple, we recommend wrapping the entire PSA process into a function. We call this function *wrapper* and use *f_* to denote that this is a function.

The wrapper function has as its inputs all the things which we may wish to vary using shiny. We set the default values to those of the base model in any report/publication. The moedl then generates PSA inputs using the *f_gen_psa* function, creates a table of results, and finally loops through the PSA, running the model with each set of PSA inputs (a row from df_psa) in turn. The function then returns the results (the costs and qalys for treatment and no treatment for each PSA run).

```{r, echo = T, eval = F}
f_wrapper <- function(
  
  #================================================================
  #                         Shiny inputs
  #================================================================
  n_age_init = 25,   # age at baseline default is 25
  n_age_max  = 110,  # maximum age of follow up default is 110
  d_r     = 0.035,   # discount rate for costs & QALYS (NICE 3.5%)
  n_sim   = 1000,    # number of simulations default 1000
  c_Trt   = 50 # cost of treatment deault 50

  ){
  
  
  
  #================================================================
  #                         Non-shiny inputs
  #================================================================
  n_t <- n_age_max - n_age_init # time horizon, number of cycles
  v_n <- c("H", "S1", "S2", "D") # the 4 health states of the model:
  n_states <- length(v_n) # number of health states 

  
  #================================================================
  #                         Create PSA Inputs
  #================================================================
  
  df_psa <- f_gen_psa(n_sim = n_sim, c_Trt)


  #================================================================
  #                         RUN PSA
  #================================================================

  # Initialize matrix of results outcomes
  df_out <- matrix(NaN, 
                 nrow = n_sim, 
                 ncol = 5,
                 dimnames = list(1:n_sim,c("Cost_NoTrt", "Cost_Trt",
                                           "QALY_NoTrt", "QALY_Trt",
                                           "ICER")))
  # loop through psa inputs running the model for each.
  for(i in 1:n_sim){
    # store results in one row of results matrix
    df_out[i,] <- f_MM_sicksicker(df_psa[i, ])
    
    # display the progress of the simulation
    cat('\r', paste(round(i/n_sim * 100), "% done", sep = " "))       
  }
  
  df_out <- as.data.frame(df_out) # convert matrix to dataframe
  
  return(df_out) # output the dataframe from the function
  
  }

```

## Integrating into Shiny

### Set-up

The set-up is relatively simple, load the shiny package from the library so that we can use the shiny functions, the source the file to load the wrapper function (which includes all model functionality).

```{r, echo = T, eval = F}

# we need the function shiny installed, this loads it from the library.
library(shiny)             

# source the wrapper function.
source("../R/wrapper.R")   

```

### User Interface

The user interface is extremely flexible, we show the code for a very simple structure (fluidpage) with a sidebar containing inputs and a main panel containing outputs. We have done very little formatting in order to minimize the quantity of code while maintaining all functionality. In order to get an aesthetically pleasing application we would have much more sophisticated formatting, relying on css, HTML and javascript.

This example user interface below is made up of two components, a titlepanel and a sidebar layout. The sidebarLayout function has within it a sidebar and a main panel. These are all contained within the *fluidpage* function which creates the ui. 

The title panel contains the title "Sick Sicker Model in Shiny", the sidebar panel contains two numeric inputs and a slider input ("Treatment Cost","PSA runs","initial age") and an Action Button ("Run / update model"). 

The values of the inputs have ids which are used by the server function, we denote these with an SI to indicate they are Shiny Inputs ("SI_c_Trt","SI_n_sim","SI_n_age_init").

The action button also has an id, this is not an input into the model wrapper (f_wrapper) so we leave out the SI and call it "run_model".

The main panel contains two objects which have been output from the server: tableOutput("icer_table") is a table of results, and plotOutput("CE_plane") is a cost-effectiveness plane plot. It is important that the format (e.g. tableOutput) matches the format of the object from the server "icer_table". The two h3() functions are simply headings which appear as "Results Table" and "Cost-effectiveness Plane".

```{r, echo = T, eval = F}

ui <- fluidPage(
  
  titlePanel("Sick Sicker Model in Shiny"),   # title of app
  
  # SIDEBAR
  sidebarLayout(    # indicates layout is going to be a sidebar-layout
    
    sidebarPanel( # open sidebar panel
    
    numericInput(inputId = "SI_c_Trt",      # id of input, used in server
                 label = "Treatment Cost",  # label next to numeric input
                 value = 200,               # initial value
                 min = 0,                   # minimum value allowed
                 max = 400),                # maximum value allowed
    
    numericInput(inputId = "SI_n_sim",      # id of input, used in server
                 label = "PSA runs",        # label next to numeric input
                 value = 1000,              # initial value
                 min = 0,                   # minimum value allowed
                 max = 400),                # maximum value allowed
    
    sliderInput(inputId = "SI_n_age_init",  # id of input, used in server
                label = "initial age",      # label next to numeric input
                value = 25,                 # initial value
                min = 10,                   # minimum value allowed
                max = 80),                  # maximum value allowed
    

    actionButton(inputId = "run_model",     # id of action button, used in server
                 label   = "Run model")     # action button label (on button)
                       
                ),  # close sidebarPanel
    
  mainPanel(                                # open main panel
    
    h3("Results Table"),                    # heading (results table)                
    
    tableOutput(outputId = "icer_table"),   # tableOutput id = icer_table, from server
    
    h3("Cost-effectiveness Plane"),         # heading (Cost effectiveness plane)
    
    plotOutput(outputId = "CE_plane")       # plotOutput id = CE_plane, from server
    
            ) # close mainpanel    
      
      ) # close sidebarlayout
                     
  ) # close UI fluidpage

```

### Server

The server is marginally more complicated than the user interface. It is created by a function with inputs and outputs.

The observe event indicates that when the action button (run_model) is pressed the code within the curly brackets is run. 

The first thing that happens is that the model wrapper function *f_wrapper* is used, with the inputs from the numeric inputs and sliders in the user interface ("SI_c_Trt","SI_n_age_init","SI_n_sim"). The *input$* prefix indicates that the objects have come from the user interface. The results of the model are stored as the dataframe object *df_model_res*.

The ICER table is then created and output (note the prefix *output$*) in the object *icer_table*. See previous section on the user interface and note that the tableOutput was reliant on "icer_table". The function renderTable rerenders the table continuously so that table always reflects the values from the data-frame of results created above. In this simple example we have created a table of results using code within the *renderTable* function. In reality we would generally use a function which creates a publication quality table which is aesthetically pleasing. There are numerous packages which provide this functionality.

The Cost-effectiveness plane is created in a similar process, using the *renderPlot* function to continuously update a plot which is created using baseR plot function using icers calcualted from the results dataframe *df_model_res*. For aesthetic purposes we recommend this is replaced by a ggplot plot which has much improved functionality.

```{r, echo = T, eval = F}

server <- function(input, output){   # server = function with two inputs
  
 observeEvent(input$run_model,       # when action button pressed ...
              ignoreNULL = F, {
    
  # use the model wrapper function to run the model with the inputs from the UI  
  df_model_res = f_wrapper(c_Trt = input$SI_c_Trt,
                           n_age_init = input$SI_n_age_init,
                           n_sim = input$SI_n_sim)
    
  
    # CREATE ICER TABLE 
    output$icer_table <- renderTable({
      
      df_res_table <- data.frame(
        
        Option =  c("Treatment","No Treatment"),
        
        QALYs  =  c(mean(df_model_res$QALY_Trt),mean(df_model_res$QALY_NoTrt)),
        
        Costs  =  c(mean(df_model_res$Cost_Trt),mean(df_model_res$Cost_NoTrt)),
        
        Inc.QALYs = c(mean(df_model_res$QALY_Trt) - mean(df_model_res$QALY_NoTrt),NA),
        
        Inc.Costs = c(mean(df_model_res$Cost_Trt) - mean(df_model_res$Cost_NoTrt),NA),
        
        ICER = c(mean(df_model_res$ICER),NA)
      )
      
      df_res_table[,2:6] <- round(df_res_table[,2:6],digits = 2) 
      
      df_res_table

      }) # table plot end.
    
    
   #  CREATE COST EFFECTIVENESS PLANE
    output$CE_plane <- renderPlot({ # render plot repeatedly updates.
      
      # calculate incremental costs and qalys from results dataframe
      df_model_res$inc_C <- df_model_res$Cost_Trt - df_model_res$Cost_NoTrt
      df_model_res$inc_Q <- df_model_res$QALY_Trt - df_model_res$QALY_NoTrt
      
      # create cost effectiveness plane plot
      plot(x = df_model_res$inc_Q,
           y = df_model_res$inc_C,
           xlab = "Incremental QALYs",
           ylab = "Incremental Costs",
           xlim = c(min(df_model_res$inc_Q,df_model_res$inc_Q*-1),
                    max(df_model_res$inc_Q,df_model_res$inc_Q*-1)),
           ylim = c(min(df_model_res$inc_C,df_model_res$inc_C*-1),
                    max(df_model_res$inc_C,df_model_res$inc_C*-1)),
           abline(h = 0,v=0)
          ) # plot end 
    }) # renderplot end
    
 }) # Observe Event End
    
    
  } # Server end
  
```  

### Running the app internally
The app can be run within the R file using the function *shinyApp* which depends on the *ui* and *server* which have been created and described above. Running this creates a shiny application in the local environment (e.g. your desktop).

```

## ----- run app------

shinyApp(ui, server)
```

### Deploying the app on the web
There are numerous methods of deploying applications on the web. Since our applications are used by a relatively small number of people we have typically used the Shinyapps.io platform provided through R-Studio. The method of deploying apps on Shinyapps.io is explained in detail [here](https://shiny.rstudio.com/articles/shinyapps.html). For apps which are expected to be used extensively or which require very large datasets it is likely to be much cheaper to use another service.

### Conclusion

The aim of this tutorial was to provide a useful reference for those hoping to create a user interface for a health economic model created in R. It is our hope that more health economic models will be created open source, and open access so that other economists can critique, learn from and adapt these models. The creation of user interfaces for these apps should improve transparency further, allowing stakeholders and third parties to conduct their own sensitivity analysis. We are particularly interested in the use of R and R-Shiny in public health decision modelling, so feel free to get in touch with any queries.

# References

Alarid-Escudero, F., Krijkamp, E.M., Enns, E.A., Hunink, M.G., Pechlivanoglou, P. and Jalal, H., 2020. Cohort state-transition models in R: From conceptualization to implementation. arXiv preprint arXiv:2001.07824.

Decision Analysis in R for Technologies in Health (DARTH) workgroup.Decision Analysis in R for Technologies in Health (DARTH) workgroup. http://darthworkgroup.com/. Accessed 19/02/2020. 

Krijkamp, E.M., Alarid-Escudero, F., Enns, E.A., Jalal, H.J., Hunink, M.M. and Pechlivanoglou, P., 2018. Microsimulation modeling for health decision sciences using R: a tutorial. Medical Decision Making, 38(3), pp.400-422.
