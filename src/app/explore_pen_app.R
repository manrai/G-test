library(reticulate)
library(plotly)
library(purrr)
library(stringr)

# py_install(c("scipy",
#              "pandas",
#              "numpy",
#              "sqlite",
#              "networkx",
#              "python-graphviz",
#              "pydot",
#              "pygraphviz"))

source_python("scripts/compute_pens.py")
pens <- compute_pens("data/bdc-test.db", "data/variants.db")
save.csv(pens, file = "data/pens.csv")

# # If you've already computed penetrances, comment out the 
# # previous three lines and uncomment this line
# pens <- read.csv("data/pens.csv")

# Jiggle discrete variables for better plotting
pens$cohort1_noisy = as.integer(pens$cohort1) + runif(nrow(pens), min = -0.2, max = 0.2)
pens$cohort2_noisy = as.integer(pens$cohort2) + runif(nrow(pens), min = -0.2, max = 0.2)
pens$cohort3_noisy = as.integer(pens$cohort3) + runif(nrow(pens), min = -0.2, max = 0.2)
pens$cohort4_noisy = as.integer(pens$cohort4) + runif(nrow(pens), min = -0.2, max = 0.2)
pens$Quadrant_noisy = as.integer(pens$Quadrant) + runif(nrow(pens), min = -0.2, max = 0.2)

create_label <- function(r) {
    return(paste("Penetrance:", toString(r[2]), "<br>",
                 "Mode:", toString(r[1]), "<br>",
                 "Cohort 1:", toString(r[3]), "<br>",
                 "Cohort 2:", toString(r[4]), "<br>",
                 "Cohort 3:", toString(r[5]), "<br>",
                 "Cohort 4:", toString(r[6]), "<br>",
                 "LVWT threshold:", toString(r[7]), "<br>",
                 "SBP threshold:", toString(r[8]), "<br>",
                 "Quadrant:", toString(r[9])
    ))
}
pens$label <- map(split(pens, seq(nrow(pens))), create_label)

fig <- plot_ly(pens, 
               x = ~`cohort1_noisy`, 
               y = ~`cohort1_noisy`,
               xaxis = list(title = ""),
               yaxis = list(title = ""),
               hoverinfo = 'text',
               text = ~paste('Penetrance: ', penetrance, '<br>',
                             'Mode: ', mode, '<br>',
                             "Cohort 1:", cohort1, "<br>",
                             "Cohort 2:", cohort2, "<br>",
                             "Cohort 3:", cohort3, "<br>",
                             "Cohort 4:", cohort4, "<br>",
                             "LVWT threshold:", `LVWT threshold`, "<br>",
                             "SBP threshold:", `SBP threshold`, "<br>",
                             "Quadrant:", Quadrant
               ),
               type = 'scatter',
               mode = 'marker',
               marker = list(color = ~penetrance,
                             size = ~penetrance,
                             opacity = 0.5,
                             sizeref = max(filter(pens, mode == "map")$penetrance) / 50,
                             showscale = TRUE),
               transforms = list(list(type = 'filter',
                                      target = ~mode,
                                      operation = '=',
                                      value = "map"
               ))
) %>%
    layout(xaxis = list(title = ""), yaxis = list(title = "")) %>%
    layout(
        updatemenus = list(
            list(type = "dropdown",
                 x = 0.5,
                 y = 1,
                 buttons = list(
                     list(method = "restyle",
                          label = "Cohort 1",
                          args = list(
                              "x", list(pens$cohort1_noisy)
                          )
                     ),
                     list(method = "restyle",
                          label = "Cohort 2",
                          args = list(
                              "x", list(pens$cohort2_noisy)
                          )
                     ),
                     list(method = "restyle",
                          label = "Cohort 3",
                          args = list(
                              "x", list(pens$cohort3_noisy)
                          )
                     ),
                     list(method = "restyle",
                          label = "Cohort 4",
                          args = list(
                              "x", list(pens$cohort4_noisy)
                          )
                     ),
                     list(method = "restyle",
                          label = "LVWT threshold",
                          args = list(
                              "x", list(pens$`LVWT threshold`)
                          )
                     ),
                     list(method = "restyle",
                          label = "SBP threshold",
                          args = list(
                              "x", list(pens$`SBP threshold`)
                          )
                     ),
                     list(method = "restyle",
                          label = "Quadrant",
                          args = list(
                              "x", list(pens$Quadrant_noisy)
                          )
                     )
                 )
            ),
            list(type = "dropdown",
                 y = 0.5,
                 buttons = list(
                     list(method = "restyle",
                          label = "Cohort 1",
                          args = list(
                              "y", list(pens$cohort1_noisy)
                          )
                     ),
                     list(method = "restyle",
                          label = "Cohort 2",
                          args = list(
                              "y", list(pens$cohort2_noisy)
                          )
                     ),
                     list(method = "restyle",
                          label = "Cohort 3",
                          args = list(
                              "y", list(pens$cohort3_noisy)
                          )
                     ),
                     list(method = "restyle",
                          label = "Cohort 4",
                          args = list(
                              "y", list(pens$cohort4_noisy)
                          )
                     ),
                     list(method = "restyle",
                          label = "LVWT threshold",
                          args = list(
                              "y", list(pens$`LVWT threshold`)
                          )
                     ),
                     list(method = "restyle",
                          label = "SBP threshold",
                          args = list(
                              "y", list(pens$`SBP threshold`)
                          )
                     ),
                     list(method = "restyle",
                          label = "Quadrant",
                          args = list(
                              "y", list(pens$Quadrant_noisy)
                          )
                     )
                 )
            ),
            list(type = "dropdown",
                 buttons = list(
                     list(method = "restyle",
                          args = list(list(
                              "transforms[0].value" = list("map"),
                              "marker.sizeref" = list(max(filter(pens, mode == "map")$penetrance) / 50)
                          )),
                          label = "MAP"),
                     list(method = "restyle",
                          args = list(list(
                              "transforms[0].value" = list("95"),
                              "marker.sizeref" = list(max(filter(pens, mode == "95")$penetrance) / 50)
                          )),
                          label = "95"),
                     list(method = "restyle",
                          args = list(list(
                              "transforms[0].value" = list("5"),
                              "marker.sizeref" = list(max(filter(pens, mode == "5")$penetrance) / 50)
                          )),
                          label = "5")
                 )
            )
        )
    )

source_python("scripts/tree.py")
graph_info <- compute_graph("data/pens.csv")

nodes2label <- function(node) {
    return(str_replace_all(node, "\\], \\[", "]<br>["))
}
graph_info$labels = lapply(graph_info$nodes, nodes2label)

fig2 <- plot_ly(
    x=graph_info$edge_x, 
    y=graph_info$edge_y,
    line=list(width=0, color="#888"),
    hoverinfo="none",
    mode="lines"
) %>%
    add_trace(
        x=graph_info$edge_x,
        y=graph_info$edge_y,
        line=list(width=1, color="#768"),
        hoverinfo="none",
        mode="lines"
    ) %>%
    add_trace(
        x=graph_info$node_x, 
        y=graph_info$node_y,
        mode='markers',
        hoverinfo="text",
        text=paste(graph_info$penetrances, "<br>",
                   graph_info$labels),
        marker=list(
            color=graph_info$penetrances,
            showscale = TRUE,
            size=30,
            line=list(width=2, color='Black')
        )) %>%
    add_trace(
        x=-80,
        y=unique(graph_info$node_y) - 70,
        mode="text",
        text=c("Cohort 1", "Cohort 2", "Cohort 3", "Cohort 4", "")
    ) %>%
    add_trace(
        x=c(600, 1050),
        y=unique(graph_info$node_y)[1] - 70,
        mode="text",
        text=c("include", "exclude")
    ) %>%
    layout(xaxis=list(zeroline = FALSE,
                      showline = FALSE,
                      showticklabels = FALSE,
                      showgrid = FALSE),
           yaxis=list(zeroline = FALSE,
                      showline = FALSE,
                      showticklabels = FALSE,
                      showgrid = FALSE),
           showlegend = FALSE
    )

#### UNDER CONSTRUCTION - try to add other hyperparams in a dropdown?
# lvwt_list = list()
# for (lt in unique(pens$`LVWT threshold`)) {
#     button = list(
#         method = "restyle",
#         args = list(list(
#             "marker.color" = lapply(graph_info$nodes, 
#                                     function(node) return(compute_penetrance_tree(node, )))
#         ))
#     )
# }

ui <- fluidPage(
    headerPanel('Penetrance Estimates'),
    mainPanel(
        plotlyOutput('plot'),
        plotlyOutput('tree')
    )
)

server <- function(input, output) {
    output$plot <- renderPlotly(
        fig
    )
    output$tree <- renderPlotly(
        fig2
    )
}

shinyApp(ui,server)