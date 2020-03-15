rm(list)

library(shiny)
library(rmarkdown)

shinyApp(
  ui = fluidPage(
    sliderInput(inputId ="slider", 
                label =  "Slider", 
                min = 1, 
                max = 100, 
                value = 50),
    numericInput(inputId = "number", 
                label = "Number", 
                value = 0.4,
                min = 0, 
                max = 1),
    downloadButton("report", "Generate report")
  ),
  server = function(input, output) {
    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "report.pdf",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- "newreport.Rmd"
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(n = input$slider, 
                       t = input$number)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, 
                          output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )
  }
)
