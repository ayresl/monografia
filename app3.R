library(shiny)
ui <- pageWithSidebar(
  headerPanel("Teste da Soma dos Postos de Wilcoxon"),
  sidebarPanel(
    numericInput(inputId = "p1",
                 label = "P( X < Y )",
                 value = 0.7,
                 min   = 0,
                 max   = 1,
                 step  = 0.05
    ),
    numericInput(inputId = "p2",
                 label = HTML("P( X < Y<sub>1</sub> ∩ X < Y<sub>2</sub> )"),
                 value = 0.8,
                 min   = 0,
                 max   = 1,
                 step  = 0.05
    ),
    numericInput(inputId = "p3",
                 label = HTML("P( X<sub>1</sub> < Y ∩ X<sub>2</sub> < Y )"),
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
    ),
    checkboxInput(inputId = "unequal",
                  label   = "tamanhos de amostra desiguais",
                  value   = FALSE
    ),
    conditionalPanel(
      condition = "input.unequal == true",
        numericInput(inputId = "a",
                     label = "razão m/n",
                     value = 1,
                     min = 0,
                     step = 0.1
        )
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
        f <- function(x) (x * (input$a * x + x + 1)) / 2 - z.a * sqrt((input$a * x ^ 2 * (input$a * x + x + 1)) / 12) - input$a * x ^ 2 * input$p1 - (x * (x + 1)) / 2 +
          z.b * sqrt(input$a * x ^ 2 * (input$p1 * (1 - input$p1) + (x - 1) * (input$p2 - input$p1 ^ 2)+(input$a * x - 1) * (input$p3 - input$p1 ^ 2)))
        n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
      }
      if (input$p1 > 0.5) {
        f <- function(x) (x * (input$a * x + x + 1)) / 2 + z.a * sqrt((input$a * x ^ 2 * (input$a * x + x + 1)) / 12) - input$a * x ^ 2 * input$p1 - (x * (x + 1)) / 2 -
          z.b * sqrt(input$a * x ^ 2 * (input$p1 * (1 - input$p1) + (x - 1)*(input$p2 - input$p1 ^ 2) + (input$a * x - 1) * (input$p3 - input$p1 ^ 2)))
        n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
      }
    }
    if (input$alternative == "less") {
      f <- function(x) (x * (input$a * x + x + 1)) / 2 - z.a * sqrt((input$a * x ^ 2 * (input$a * x + x + 1)) / 12) - input$a * x ^ 2 * input$p1 - (x * (x + 1)) / 2 +
        z.b * sqrt(input$a * x ^ 2 * (input$p1 * (1 - input$p1) + (x - 1) * (input$p2 - input$p1 ^ 2) + (input$a * x - 1) * (input$p3 - input$p1 ^ 2)))
      n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
    }
    if (input$alternative == "greater") {
      f <- function(x) (x * (input$a * x + x + 1)) / 2 + z.a * sqrt((input$a * x ^ 2 * (input$a * x + x + 1)) / 12) - input$a * x ^ 2 * input$p1 - (x * (x + 1)) / 2 -
        z.b * sqrt(input$a * x ^ 2 * (input$p1 * (1 - input$p1) + (x - 1) * (input$p2 - input$p1 ^ 2) + (input$a * x - 1) * (input$p3 - input$p1 ^ 2)))
      n <- uniroot(f, c(2, 1e+07), tol = .Machine$double.eps ^ 0.25)$root
    }
    if (input$unequal == TRUE) {
      paste0("Recomenda-se que as amostras tenham ", ceiling(input$a * n), " e ", ceiling(n), " elementos.")
    }
    else {
      paste0("Recomenda-se que as amostras tenham ", ceiling(n), " elementos cada.")
    }
  })
}
shinyApp(ui = ui, server = server)