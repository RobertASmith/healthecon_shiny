# ==============
# Robert Smith & Paul Schneider
# University of Sheffield
# contact: rasmith3@sheffield.ac.uk
# ==============

## app.R ##
library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)

ui <- dashboardPage(
  
  # HEADER
  dashboardHeader(title = "Sick-Sicker Model",  titleWidth = 300,
                  
                  dropdownMenu(type = "notifications",badgeStatus = NULL,icon=icon("info"),
                               notificationItem(text = "Info text display",icon = icon(""))
                  )
                  
  ),
  
  # SIDEBAR
  dashboardSidebar(width = 300,
                   sidebarMenu(
                     tags$head( 
                       tags$style(HTML(".treeview-menu  { 
                                       font-size: 12px; 
                                       padding-top:0px;
                                       padding-bottom:0px;
                                       }
                                       .main-sidebar { 
                                       font-size: 17px; 
                                       }
                                       #inline_input  label{ 
                                       display: table-cell; 
                                       text-align: center; 
                                       vertical-align: middle; 
                                       margin-right:100px
                                       padding-right:100px
                                       } 
                                       ")) 
                       
                       ),
                     
                     # menu itmes: cost inputs
                     menuItem("Costs", tabName="Costs",
                              HTML("<br>"),h4("Costs input"),
                              fluidRow(
                                column(offset = 1,width = 4,div(style = "margin-top:15px"),h4("Costs 1")),
                                column(width = 6, numericInput(inputId = "i.costs1",label = NULL,value = 200.42))
                              ),
                              fluidRow(
                                column(offset = 1,width = 4,div(style = "margin-top:15px"),h4("Costs 2")),
                                column(width = 6, numericInput(inputId = "i.costs2",label = NULL,value = 200.44))
                              ),
                              fluidRow(
                                column(offset = 1,width = 4,div(style = "margin-top:15px"),h4("Costs 3")),
                                column(width = 6, numericInput(inputId = "i.costs3",label = NULL,value = 200.44))
                              ),
                              fluidRow(
                                column(offset = 1,width = 4,div(style = "margin-top:15px"),h4("Costs 4")),
                                column(width = 6, numericInput(inputId = "i.costs4",label = NULL,value = 191.63))
                              )
                     ),
                     
                     # menu itmes: other inputs
                     menuItem("Other", tabName="Other",
                              HTML("<br>"),h4("Other model inputs"),
                              fluidRow(column(offset = 1,width = 11,
                                              sliderInput(inputId = "i.age",label = "Initial Age",min = 18,max=100,value = 70),
                                              sliderInput(inputId = "i.cyclelen",label = "Time horizon (years)",min = 1,max=100,value = 30),
                                              sliderInput(inputId = "i.samples",label = "PSA samples",min = 1,max=5000,value = 10)
                              ))),
                     
                     # menu itme: action button
                     hr(),
                     actionButton("run_model","Run / update model")
                       )
                       ),
  
  
  ### ---  DASHBOARD BODY ------------------
  dashboardBody(
    fluidRow(
      # table box
      box(width = 12, dataTableOutput("icer_table"))#,
      
      # # highlights box
      # tabBox(width = 4,
      #            id = "tabset1", # height = "250px",
      #            tabPanel("Tab1"),
      #            tabPanel("Tab2", "Tab content 2")
      #     )
    ),
    # ce-plane box
    fluidRow(
      tabBox(width = 9,title = "",
             id = "tabset1", # height = "250px",
             
             tabPanel("CE Plane", 
                      plotOutput("CE_plane",height = 450),
                      fluidPage(downloadButton('cep'))), # end of tab panel
             
             tabPanel("CEAC", 
                      plotOutput("CEAC", height = 450),
                      fluidPage(downloadButton('ceac')))
      ), # box1 end
      box(width = 3,br(),
          radioButtons("i.comparitor","Reference",
                       selected = "Coumarin (INR 2-3)",
                       choices = treatment_choices
          ),
          checkboxGroupInput("i.treatments","Treatments",
                             choices = treatment_choices),
          numericInput(inputId = "i.thresh",label = "WTP Threshold",value = 20000),
          checkboxInput(inputId = "i.ellipse",label = "Ellipse",value = F)
          
      ) # box2 end
      
    ) # end of row
    
  ) # dashboard body close
                     ) # UI close


server <- function(input, output,session) {
  
  
  observeEvent(input$run_model,ignoreNULL = F, {
    
        model_res = wrapper()
      }
      model_res
    }) # actionbutton end
    
    # ICER TABLE
    output$icer_table <- renderDataTable({
      display_results.table(model_output = model_res()$model.outputs)
    })
    
    
    # CE PLANCE
    output$CE_plane <- renderPlot({
      
      p1 = plot_ce.plane.plot(model_output = model_res()$model.outputs,
                              comparitor = input$i.comparitor,
                              treatment = input$i.treatments,
                              thresh = input$i.thresh,
                              show_ellipse = input$i.ellipse)
      
      p1
      
    })
    
    
    # CE-PLANE DOWNLOAD FUNCTIONALITY
    ceP_download <- plot_ce.plane.plot(model_output = model_res()$model.outputs)
    
    output$cep = downloadHandler(
      filename = 'ce_plane.png',
      content = function(file) {
        device <- function(..., width, height) {
          grDevices::png(..., width = width, height = height,
                         res = 300, units = "in")
        }
        ggsave(file, 
               plot = ceP_download, 
               device = device)
      })
    
    
    # ----------------------------
    # NB: DO WE WANT TO USE CHECKBOX INPUTS FOR CEAC TOO?
    
    # CEAC 
    output$CEAC <- renderPlot({
      plot_ceac.plot(model_output = model_res()$model.outputs)
      
    })
    
    
    # CEAC DOWNLOAD FUNCTIONALITY
    ceac_download <- plot_ceac.plot(model_output = model_res()$model.outputs)
    
    output$ceac = downloadHandler(
      filename = 'ceac.png',
      content = function(file) {
        device <- function(..., width, height) {
          grDevices::png(..., width = width, height = height,
                         res = 300, units = "in")
        }
        ggsave(file, 
               plot = ceac_download, 
               device = device)
      })
    
    
    
  }) # Server end
  
  
}


## ----- run app------

shinyApp(ui, server)