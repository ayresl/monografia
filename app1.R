library(shiny)
ui <- pageWithSidebar(
headerPanel("Teste de Aderência do Qui-Quadrado"),
sidebarPanel(
numericInput(inputId = "w",
label = "tamanho do efeito",
value = 0.3,
min   = 0,
step  = 0.05
),
numericInput(inputId = "k",
label = "número de categorias",
value = 7,
min   = 1,
step  = 1
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
)
),
mainPanel(
h3("O Tamanho da Amostra"),
textOutput("resposta")
)
)
server <- function(input, output) {
output$resposta <- renderText({
f <- function(x) qchisq(input$sig.level, input$k - 1, lower.tail = FALSE) - qchisq(input$power, input$k - 1, ncp = x, lower.tail = FALSE)
lambda <- uniroot(f, c(1 + 1e-10, 1e+04))$root
res <- lambda / input$w ^ 2
paste0("Recomenda-se que a amostra tenha ", ceiling(res), " elementos.")
})
}
shinyApp(ui = ui, server = server)