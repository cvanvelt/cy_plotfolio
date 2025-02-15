centroidApp <- function(cl.center.df, edges.df = NULL, show.size = TRUE, max.size = 5) {
  library(shiny)
  library(plotly)
  library(ggplot2)
  library(viridis)
  library(knitr)
  library(dplyr)
  
  dataset <- cl.center.df
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput('xvar', 'X', names(dataset), selected = 'x'),
        selectInput('yvar', 'Y', names(dataset), selected = 'y'),
        selectInput('level', 'GROUP', names(dataset), selected = 'subclass_label'),
        selectInput('color', 'COLOR', names(dataset), selected = 'subclass_color'),
        textInput('NewGroup', 'Enter new group designation'),
        actionButton('Change', 'Change Group Assignment'),
        actionButton('exit', 'Return to R and write data')
      ),
      mainPanel(
        plotlyOutput('plot', width = '100%', height = '800px'),
        verbatimTextOutput('brush')
      )
    )
  )
  
  server <- function(input, output, session) {
    values <- reactiveValues(vv = NULL, dataset = dataset)
    
    data.sel <- reactive({
      req(input$xvar, input$yvar, input$level, input$color)
      values$dataset[, c(input$xvar, input$yvar, input$level, input$color)]
    })
    
    edges.sel <- reactive({
      req(poly.Edges)
      poly.Edges
    })
    
    output$plot <- renderPlotly({
      g1 <- data.sel()
      
      p <- ggplot() + 
        # Add edges as polygons (grouped by "Group")
        geom_polygon(data = edges.df, aes(x = x, y = y, group = Group), 
                     fill = "gray70", alpha = 0.4) +
        
        # Add nodes (points)
        geom_point(data = g1, aes_string(x = input$xvar, y = input$yvar, 
                                         color = input$color, 
                                         size = if ("cluster_size" %in% names(dataset)) dataset$cluster_size else 1), 
                   alpha = 0.8, shape = 19,show_guide = FALSE ) +
        scale_size_area(trans = "sqrt", max_size = max.size) +
        labs(x = input$xvar, y = input$yvar) +
        coord_fixed() +
        theme_bw()
      
      # Apply color scale conditions
      if (is.numeric(g1[[input$color]])) {
        p <- p + scale_color_viridis(option = "magma")
      } else if (all(grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", g1[[input$color]]))) {
        p <- p + scale_color_identity()
      }
      
      # Convert to plotly
      plotly_obj <- ggplotly(p, source = "plot_source") %>% 
        layout(dragmode = 'lasso')
      
      event_register(plotly_obj, 'plotly_selected')  # Ensure selection works
      
      plotly_obj
    })

    output$brush <- renderPrint({
      d <- event_data("plotly_selected", source = "plot_source")
      if (is.null(d)) return("Click and drag events (i.e., select/lasso) appear here (double-click to clear)")
      
      selected_points <- data.frame(
        x = round(d$x, 3),
        y = round(d$y, 3)
      )
      
      selected_data <- values$dataset %>%
        mutate(x = round(.data[[input$xvar]], 3), y = round(.data[[input$yvar]], 3)) %>%
        inner_join(selected_points, by = c("x", "y"))
      
      if (nrow(selected_data) > 0) {
        values$vv <- selected_data
        print(kable(selected_data))
      } else {
        "No points selected."
      }
    })
    
    observeEvent(input$Change, {
      req(values$vv, input$NewGroup)
      row_indices <- rownames(values$dataset) %in% rownames(values$vv)
      values$dataset[row_indices, "new"] <<- input$NewGroup
    })
    
    observeEvent(input$exit, {
      stopApp(values$dataset)
    })
  }
  
  runApp(shinyApp(ui, server))
}
