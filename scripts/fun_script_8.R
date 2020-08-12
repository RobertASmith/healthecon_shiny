# plot function

ce_plot <- function(results = df_model_res){

  # calculate incremental costs and qalys from results data-frame in function input.
df_plot <- data.frame(inc_C = results$Cost_Trt - results$Cost_NoTrt,
                      inc_Q = results$QALY_Trt - results$QALY_NoTrt)

# now use this plotting data-frame to create a very simple ggplot.
plot <- ggplot(data = df_plot,
               aes(x = inc_Q,  # x axis incremental QALYS
                   y = inc_C)  # y axis incremental Costs
               ) +
               
               geom_point() +
                 
                 # titles
                 labs(
                   title = "Cost-effectiveness Plane",
                   subtitle = paste0("Results of Probabilistic Sensitivity Analysis:")
                 ) +
                 
                 # set xlimits and ylimits for plot.
                 xlim(c(
                   min(df_plot$inc_Q, df_plot$inc_Q * -1),
                   max(df_plot$inc_Q, df_plot$inc_Q * -1)
                 )) +
                 
                 ylim(c(
                   min(df_plot$inc_C, df_plot$inc_C * -1),
                   max(df_plot$inc_C, df_plot$inc_C * -1)
                 ))
               
plot # output the plot from the function.

}