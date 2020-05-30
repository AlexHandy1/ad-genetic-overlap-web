#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Step 3 AD PRS to Protein Results"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("apoe", strong("APOE"), 
                        choices = c("With-APOE","No-APOE"), selected = "With-APOE")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("prs")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    #Load data
    with_apoe_res <- read.csv("Step3_PRS_Results_With_APOE.csv", header=T)
    no_apoe_res <- read.csv("Step3_PRS_Results_No_APOE.csv", header=T)
    
    #Prepare data
    with_apoe_res <- with_apoe_res %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
    names(with_apoe_res)[1] <- "Protein"
    no_apoe_res <- no_apoe_res %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
    names(no_apoe_res)[1] <- "Protein"
    
    with_apoe_res <- with_apoe_res %>% mutate(P_MinusLog10 = -log10(P), Significant = if_else(P < 0.05, "Y","N"))
    no_apoe_res <- no_apoe_res %>% mutate(P_MinusLog10 = -log10(P), Significant = if_else(P < 0.05, "Y","N"))
    
    output$prs <- renderPlot({
        if(input$apoe == "With-APOE") {
            ggplot(with_apoe_res, aes(Protein, P_MinusLog10, label=Threshold, fill=Significant)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "red") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="blue", "N"="gray" ), guide = FALSE )
        } else {
            ggplot(no_apoe_res, aes(Protein, P_MinusLog10, label=Threshold, fill=Significant)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "red") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="blue", "N"="gray" ), guide = FALSE )
        }
    })
    
    # output$table <- renderTable({
    #     tableInput <- if(input$apoe == "With-APOE") {
    #         with_apoe_res   
    #     } else {
    #         no_apoe_res
    #     }
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)
