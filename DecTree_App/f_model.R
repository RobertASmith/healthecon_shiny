# DECISION TREE MODEL WRAPPER FUNCTION

# Robert Smith, University of Sheffield, 18.03.2020
# Paul Schneider, University of Sheffield 11.03.2020

# Aim: Simplest model possible with uncertainty included.


f_model = function(n_draws     = 1000,      # number of iterations
                   c_B_tr_mean = 1000,      # mean for distribution of costs of B
                   c_B_tr_sd   = 100,       # standard deviation for distribution of costs of B
                   rr_mean     = 0.95,      # mean for distribution of 'relative risk' of B vs A
                   rr_sd       = 0.01,      # standard deviation for distribution of rr 
                   p_A_S1      = 0.05       # probability of S1 with treatment A
) {              
  
  # -- Transition Probabilities -- #   
  
  p_A_S2 = 1 - p_A_S1     # probability of S2 with treatment A     
  
  # vector of relative risk from random log-normal distribution.
  rr = rlnorm(n = n_draws,
              meanlog =  log(rr_mean^2 / sqrt(rr_sd^2 + rr_mean^2)),
              sdlog =    sqrt(log(1 + (rr_sd^2 / rr_mean^2)))
  )
  
  p_B_S1 = p_A_S1 * rr    # probability of S1 with treatment B
  
  p_B_S2 = 1-p_B_S1       # prob. of S2 with treatment B
  
  # -- Costs -- #
  
  c_S1 = 100             # cost of S1
  c_S2 = 300             # cost of S2
  
  # One off costs for treatment B are from a normal distribution.
  c_B_tr = rnorm(n = n_draws, 
                 mean = c_B_tr_mean,
                 sd = c_B_tr_sd) 
  
  # -- Life Years -- #
  
  ly_S1 = 3              # life expectancy in S1
  ly_S2 = 30             # life expectancy in S2
  
  # -- Calculate Outcomes -- #
  
  # Costs
  c_A = p_A_S1 * c_S1 + p_A_S2 * c_S2                     # total cost under A
  c_B = p_B_S1 * (c_S1+c_B_tr) + p_B_S2 * (c_S2+c_B_tr)   # total cost under B
  c_incr = c_B - c_A                                      # incremental costs
  
  u_A =  p_A_S1 * ly_S1 + p_A_S2 * ly_S2                      # total life years under A
  u_B =  p_B_S1 * ly_S1 + p_B_S2 * ly_S2                      # total life years under B
  u_incr = u_B - u_A                                          # incremental life years
  
  # -- Calc ICER -- #  
  icer = c_incr / u_incr   # resulting ICER
  
  df_res = data.frame(u_A, u_B, c_A, c_B, c_incr, u_incr, icer) # convert to dataframe
  
  return(df_res)   # return dataframe.
}
