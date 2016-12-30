library(shiny)

shinyUI(fluidPage(
    # headerPanel("Plot for GenePix"),
    titlePanel(h1("Plot for GenePix")),

    fluidRow(
        column(
            4,

            wellPanel(

                radioButtons("choose_source",
                            label = "Choose data source",
                            choices = list("Pre-loaded data"="p", "choose your own file"="o"),
                            selected = "p"),

                uiOutput("choose_dataset"),


                uiOutput("choose_columns"),

                numericInput("block",
                             label = h4("Block:"),
                             value = 1)


            ),
            # br(),
            tabsetPanel(

                tabPanel("Simple plot",
                    uiOutput("choose_coltype")
                ),

                tabPanel("Heatmap",
                    selectInput("useRank",
                                label = "Use rank matrix or not",
                                choices = list("No rank, use original" = "n",
                                               "Rank whole matrix" = "a",
                                               "Rank with sample" = "s"),
                                selected = "n")
                )

            )
        ),
        column(8,
               tabsetPanel(

                   tabPanel("Simple plot",
                            fluidRow(
                                column(12,
                                       h2("Main plot"), plotOutput("plot1", height = "640px")
                                ),
                                column(12,
                                       h2("Block plot"), plotOutput("plot2",height = "480px")
                                )
                            )
                   ),

                   tabPanel("Heatmap",
                            fluidRow(
                                column(12,
                                       h2("Main plot"), plotOutput("heat1", height = "640px")
                                ),
                                column(12,
                                       h2("Block plot"), plotOutput("heat2",height = "480px")
                                )
                            )

                   )

               )

        )
    )

))


