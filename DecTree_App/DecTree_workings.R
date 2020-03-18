
# rm(list=ls())     


# ------------------------------------------------------------------------------
# 1. DETERMINISTIC SIMPLE CODE
# probabilities    
p_A_S1 = 0.1            # prob. of S1 with treatment A
p_A_S2 = 1-p_A_S1       # prob. of S2 with treatment A

rr = 0.95               # Relative risk of S1 with B, comapred to A
p_B_S1 = p_A_S1 * rr    # prob. of S1 with treatment B
p_B_S2 = 1-p_B_S1       # prob. of S2 with treatment B

# costs
c_S1 = 100             # cost of S1
c_S2 = 300             # cost of S2
c_B_tr = 1000          # one off cost for treatment B

# life years
ly_S1 = 3              # life expectancy in S1
ly_S2 = 10             # life expectancy in S2

# outcomes
c_A = p_A_S1 * c_S1 + p_A_S2 * c_S2                     # total cost under A
c_B = p_B_S1 * (c_S1+c_B_tr) + p_B_S2 * (c_S2+c_B_tr)   # total cost under B

u_A =  p_A_S1 * ly_S1 + p_A_S2 * ly_S2                      # total life years under A
u_B =  p_B_S1 * ly_S1 + p_B_S2 * ly_S2                      # total life years under B

icer = (c_B - u_A) / (u_B -u_A)   # resulting ICER

icer


# ------------------------------------------------------------------------------
# 2. PROBABILISTIC SIMPLE CODE WITH LOOP

n_draws = 1000   # number of PSA iterations

# probabilities   
p_A_S1 = rbeta(n_draws,100,900)            # prob. of S1 with treatment A
p_A_S2 = 1-p_A_S1       # prob. of S2 with treatment A

rr = 0.95               # Mean Relative risk under B
s = 0.03               # standard deviation of rr
# drawing from random normal
location <- log(rr^2 / sqrt(s^2 + rr^2))
shape <- sqrt(log(1 + (s^2 / rr^2)))
rr = rlnorm(n=n_draws, location, shape)

p_B_S1 = p_A_S1 * rr    # prob. of S1 with treatment B
p_B_S2 = 1-p_B_S1       # prob. of S2 with treatment B

# costs
c_S1 = rgamma(n_draws,shape= 100,scale = 100)             # cost of S1
c_S2 = rgamma(n_draws,shape= 300,scale = 100)             # cost of S2
c_B_tr = 1000          # one off cost for treatment B

# life years
ly_S1 = 3              # life expectancy in S1
ly_S2 = 10             # life expectancy in S2

# outcomes
c_A = p_A_S1 * c_S1 + p_A_S2 * c_S2                     # total cost under A
c_B = p_B_S1 * (c_S1+c_B_tr) + p_B_S2 * (c_S2+c_B_tr)   # total cost under B

u_A =  p_A_S1 * ly_S1 + p_A_S2 * ly_S2                      # total life years under A
u_B =  p_B_S1 * ly_S1 + p_B_S2 * ly_S2                      # total life years under B

icer = (c_B - u_A) / (u_B -u_A)   # resulting ICER

mean(icer)



# ------------------------------------------------------------------------------
# 3. WRAPPER

f_model = function(n_draws = 1000,      # number of PSA iterations
                   c_B_tr_mean = 1000,  # min costs of B
                   c_B_tr_sd = 100,     # max costs of B
                   rr = 0.95,           # Mean Relative risk under B
                   s = 0.03,            # standard deviation of rr
                   p_A_S1 = 0.05        # prob. of S1 with treatment A
) {              
  # probabilities   
  # p_A_S1 = 0.05           # prob. of S1 with treatment A
  p_A_S2 = 1-p_A_S1       # prob. of S2 with treatment A
  
  # rr = 0.95               # Mean Relative risk under B
  # s = 0.03               # standard deviation of rr
  # drawing from random normal
  location <- log(rr^2 / sqrt(s^2 + rr^2))
  shape <- sqrt(log(1 + (s^2 / rr^2)))
  rr = rlnorm(n=n_draws, location, shape)
  
  p_B_S1 = p_A_S1 * rr    # prob. of S1 with treatment B
  p_B_S2 = 1-p_B_S1       # prob. of S2 with treatment B
  
  # costs
  c_S1 = 100             # cost of S1
  c_S2 = 300             # cost of S2
  c_B_tr = rnorm(n = n_draws,mean = c_B_tr_mean,sd = c_B_tr_sd) # one off cost for treatment B
  
  # life years
  ly_S1 = 3              # life expectancy in S1
  ly_S2 = 30             # life expectancy in S2
  
  # outcomes
  c_A = p_A_S1 * c_S1 + p_A_S2 * c_S2                     # total cost under A
  c_B = p_B_S1 * (c_S1+c_B_tr) + p_B_S2 * (c_S2+c_B_tr)   # total cost under B
  c_incr = c_B-c_A
  
  u_A =  p_A_S1 * ly_S1 + p_A_S2 * ly_S2                      # total life years under A
  u_B =  p_B_S1 * ly_S1 + p_B_S2 * ly_S2                      # total life years under B
  u_incr = u_B-u_A
  
  icer = c_incr / u_incr   # resulting ICER
  # mean(icer)
  
  df_res = data.frame(u_A, u_B, c_A, c_B, c_incr, u_incr, icer)
  return(df_res)
}


df_model_res = f_model(n_draws = 1000,
                 rr = 0.95,s = 0.01,
                 c_B_tr_mean = 1000,c_B_tr_sd = 200,
                 p_A_S1 = 0.05)

plot(x=df_model_res$u_incr,
     y=df_model_res$c_incr,
     cex=0.8,
     col="cadetblue4",
     xlab = "Incremental Life Years",
     ylab = "incremental Costs",
     xlim = range(c(0,df_model_res$u_incr)),
     ylim = range(c(0,df_model_res$c_incr)))
abline(h=0,v=0)
abline(coef= c(0,20000),lty=2)






#   
#   
# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
#     
#     # Create DecTree data.frame for auto + plot
#      DecTree = data.frame(level_1 = c(".",".",".","."),
#                     level_2 = c("A","A","B","B"),
#                     level_3 = c("State 1","State 2","State 1","State 2"),
#                     p = c(0.1,0.9,0.2,0.8),
#                     cost = c(1000+500,2000+500,1000,2000),
#                     lyears = c(5,1,5,1)
#                     )
#      DecTree$pathString = paste(DecTree$level_1,DecTree$level_2,DecTree$level_3,sep="/")
#      DecTreeN <- FromDataFrameTable(DecTree)
#      # print(DecTreeN, "p", "cost","lyears")
#      # plot(DecTreeN)  
#      
#      SetGraphStyle(DecTreeN, rankdir = "LR")
#      SetEdgeStyle(DecTreeN, arrowhead = "vee")
#      Do(DecTreeN$leaves, function(node) SetEdgeStyle(node,label = paste0(node$p)))
#      node_label_maker = function(node){paste0(node$name,"\n",
#                                            "cost =", node$cost,"\n",
#                                            "life years =", node$lyears,"\n")}
#      
#      Do(DecTreeN$leaves, function(node) SetNodeStyle(node,
#                                                      label = node_label_maker(node),
#                                                      fontsize=6,
#                                                      shape="circle",fixedsize=T))
#      plot(DecTreeN)