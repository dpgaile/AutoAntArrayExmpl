library(shiny)
library(AutoAntigen)
library(datasets)
library(gplots)


shinyServer(function(input, output){

    # Slide1_IgG=read.table("data/Slide1_532.gpr",sep="\t",skip=32,header=T)
    # Slide1_IgM=read.table("data/Slide1_635.gpr",sep="\t",skip=32,header=T)
    # Slide2_IgG=read.table("data/Slide2_532.gpr",sep="\t",skip=32,header=T)
    # Slide2_IgM=read.table("data/Slide2_635.gpr",sep="\t",skip=32,header=T)

    data_sets <- c("Slide1_IgG","Slide1_IgM","Slide2_IgG","Slide2_IgM")

    # Drop-down selection box for which data set
    output$choose_dataset <- renderUI({
        selectInput("dataset", "Data set", as.list(data_sets))
    })

    # Choose column based on dataset
    output$choose_columns <- renderUI({
        # If missing input, return to avoid error later in function
        if(is.null(input$dataset))
            return()

        # Or choose your own file
        inFile <- input$uploadfile

        # Get the data set with the appropriate name
        if (is.null(inFile)) {
            dat <- get(input$dataset)
            #dat <- get("Slide1_IgG")
        } else {
            dat <- read.table(inFile$datapath,sep="\t",skip=32,header=T)
        }

        colnames <- names(dat)[9:37]

        # Create the menu and select the first one by default
        selectInput("columns", "Choose columns",
                    choices  = colnames,
                    selected = colnames[1])
    })

    # Choose coltype based on dataset
    output$choose_coltype <- renderUI({
        # If missing input, return to avoid error later in function
        if(is.null(input$dataset))
            return()

        # Or choose your own file
        inFile <- input$uploadfile

        # Get the data set
        if (is.null(inFile)) {
            dat <- get(input$dataset)
        } else {
            dat <- read.table(inFile$datapath,sep="\t",skip=32,header=T)
        }

        coltype <- ifelse(grepl("532", colnames(dat)[9]), "green", "red")

        # Create the menu and select the correct one by default
        selectInput("coltype", "Choose color",
                    choices  = c("green", "red"),
                    selected = coltype)
    })

    # datasetInput <- reactive({
    #     switch(input$dataset,
    #            "Slide1_IgG"=rawS1g,
    #            "Slide1_IgM"=rawS1r,
    #            "Slide2_IgG"=rawS2g,
    #            "Slide2_IgM"=rawS2r
    #            )
    # })


    # output$colSelector <- renderUI({
    #     col <- colnames(datasetInput())[9:37]
    #     selectInput("col", "Choose column:", col)
    # })


    output$plot1 <- renderPlot({
        # If missing input, return to avoid error later in function
        if(is.null(input$dataset))
            return()

        # Or choose your own file
        inFile <- input$uploadfile

        # Get the data set
        if (is.null(inFile)) {
            dat <- get(input$dataset)
        } else {
            dat <- read.table(inFile$datapath,sep="\t",skip=32,header=T)
        }

        # Make sure columns are correct for data set (when data set changes, the
        # columns will initially be for the previous data set)
        if (is.null(input$columns) || !(input$columns %in% names(dat)))
            return()

        QuickPlot(dataset=dat,
                   coltype=input$coltype,
                   col=input$columns
                  )
    })

    output$plot2 <- renderPlot({
        # If missing input, return to avoid error later in function
        if(is.null(input$dataset))
            return()

        # Or choose your own file
        inFile <- input$uploadfile

        # Get the data set
        if (is.null(inFile)) {
            dat <- get(input$dataset)
        } else {
            dat <- read.table(inFile$datapath,sep="\t",skip=32,header=T)
        }

        # Make sure columns are correct for data set (when data set changes, the
        # columns will initially be for the previous data set)
        if (is.null(input$columns) || !(input$columns %in% names(dat)))
            return()

        quickPlot(dataset=dat,
                  coltype=input$coltype,
                  col=input$columns,
                  block = input$block
                  )
    })


    # Heatmap

    output$heat1 <- renderPlot({
        # If missing input, return to avoid error later in function
        if(is.null(input$dataset))
            return()

        # Or choose your own file
        inFile <- input$uploadfile

        # Get the data set
        if (is.null(inFile)) {
            dat <- get(input$dataset)
        } else {
            dat <- read.table(inFile$datapath,sep="\t",skip=32,header=T)
        }

        # Make sure columns are correct for data set (when data set changes, the
        # columns will initially be for the previous data set)
        if (is.null(input$columns) || !(input$columns %in% names(dat)))
            return()

        QuickHeat(dataset=dat,
                  useRank=input$useRank,
                  col=input$columns
        )
    })


    output$heat2 <- renderPlot({
        # If missing input, return to avoid error later in function
        if(is.null(input$dataset))
            return()

        # Or choose your own file
        inFile <- input$uploadfile

        # Get the data set
        if (is.null(inFile)) {
            dat <- get(input$dataset)
        } else {
            dat <- read.table(inFile$datapath,sep="\t",skip=32,header=T)
        }

        # Make sure columns are correct for data set (when data set changes, the
        # columns will initially be for the previous data set)
        if (is.null(input$columns) || !(input$columns %in% names(dat)))
            return()

        quickHeat(dataset=dat,
                  useRank=input$useRank,
                  col=input$columns,
                  block = input$block
        )
    })



})
