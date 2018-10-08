library(shiny)
ui <- pageWithSidebar(
  headerPanel("Teste dos Postos Sinalizados de Wilcoxon"),
  sidebarPanel(
    numericInput(inputId = "p1",
                 label = "P(X>0)",
                 value = 0.7,
                 min   = 0,
                 max   = 1,
                 step  = 0.05
    ),
    numericInput(inputId = "p2",
                 label = "P(X1+X2>0)",
                 value = 0.8,
                 min   = 0,
                 max   = 1,
                 step  = 0.05
    ),
    numericInput(inputId = "p3",
                 label = "P(X1+X2>0 & X1+X3>0)",
                 value = 0.75,
                 min   = 0,
                 max   = 1,
                 step  = 0.05
    ),
    numericInput(inputId = "sig.level",
                 label = "nível de significância",
                 value = 0.05,
                 min   = 0,
                 max   = 1,
                 step  = 0.01
    ),
    numericInput(inputId = "power",
                 label = "poder",
                 value = 0.8,
                 min   = 0,
                 max   = 1,
                 step  = 0.01
    ),
    selectInput(inputId = "alternative",
                label   = "tipo:",
                c("bilateral" = "two.sided", "unilateral inferior" = "less", "unilateral superior" = "greater")
    )
  ),
  mainPanel(
    h3("O Tamanho da Amostra"),
    textOutput("resposta")
  )
)
server <- function(input, output) {
  output$resposta <- renderText({
    if (input$alternative == "two.sided") {
      z.a <- qnorm(1 - input$sig.level / 2)
    }
    else {
      z.a <- qnorm(1 - input$sig.level)
    }
    z.b <- qnorm(1 - input$power)
    if (input$alternative == "two.sided") {
      if (input$p1 < 0.5) {
        f <- function(x) (x * (x + 1)) / 4 - z.a * sqrt((x * (x + 1) * (2 * x + 1)) / 24) - x * (input$p1 + (x - 1) / 2 * input$p2) +
          z.b * sqrt((x * input$p1 * (1 - input$p1) + (x * (x - 1)) / 2 * (2 * (input$p1 - input$p2) ^ 2 + 3 * input$p2 * (1 - input$p2)) + x * (x - 1) * (x - 2) * (input$p3 - input$p2 ^ 2)))
        n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
      }
      if (input$p1 > 0.5) {
        f <- function(x) (x * (x + 1)) / 4 + z.a * sqrt((x * (x + 1) * (2 * x + 1)) / 24) - x * (input$p1 + (x - 1) / 2 * input$p2) -
          z.b * sqrt((x * input$p1 * (1 - input$p1) + (x * (x - 1)) / 2 * (2 * (input$p1 - input$p2) ^ 2 + 3 * input$p2 * (1 - input$p2)) + x * (x - 1) * (x - 2) * (input$p3 - input$p2 ^ 2)))
        n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
      }
    }
    if (input$alternative == "less") {
      f <- function(x) (x * (x + 1)) / 4 - z.a * sqrt((x * (x + 1) * (2 * x + 1)) / 24) - x * (input$p1 + (x - 1) / 2 * input$p2) +
        z.b * sqrt((x * input$p1 * (1 - input$p1) + (x * (x - 1)) / 2 * (2 * (input$p1 - input$p2) ^ 2 + 3 * input$p2 * (1 - input$p2)) + x * (x - 1) * (x - 2) * (input$p3 - input$p2 ^ 2)))
      n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
    }
    if (input$alternative == "greater") {
      f <- function(x) (x * (x + 1)) / 4 + z.a * sqrt((x * (x + 1) * (2 * x + 1)) / 24) - x * (input$p1 + (x - 1) / 2 * input$p2) -
        z.b * sqrt((x * input$p1 * (1 - input$p1) + (x * (x - 1)) / 2 * (2 * (input$p1 - input$p2) ^ 2 + 3 * input$p2 * (1 - input$p2)) + x * (x - 1) * (x - 2) * (input$p3 - input$p2 ^ 2)))
      n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
    }
    paste0("Recomenda-se que a amostra tenha ", ceiling(n), " elementos.")
  })
}
shinyApp(ui = ui, server = server)