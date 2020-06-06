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
    titlePanel(""),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("step", strong("Analysis Step"), 
                        choices = c("2","3"), selected = "2"),
            
            conditionalPanel(
                condition = "input.step == '2'",
                selectInput("chart", strong("Chart"), 
                            choices = c("heritability","individual-sample-prs","meta-analysis-prs"), selected = "heritability"),
            ), 
            
            conditionalPanel(
                condition = "input.step == '2' || input.step == '3'",
                selectInput("apoe", strong("APOE"), 
                            choices = c("With-APOE","No-APOE"), selected = "With-APOE"),
            ), 
            
            conditionalPanel(
                condition = "input.step == '2' || input.step == '3'",
                selectInput("group", strong("Group"), 
                            choices = c("All","Males", "Females", "65.And.Over", "70.And.Over", "75.And.Over", "80.And.Over"), selected = "All"),
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
           conditionalPanel(
               condition = "input.step == '2' & input.chart == 'heritability'",
               plotOutput("h2")
           ),
           
           conditionalPanel(
               condition = "input.step == '2' & input.chart == 'individual-sample-prs'",
               plotOutput("prs_indiv")
           ),
           
           conditionalPanel(
               condition = "input.step == '2' & input.chart == 'meta-analysis-prs'",
               plotOutput("prs_meta")
           )
        )
    )
)

server <- function(input, output) {
    
    #DATA PREPARATION#
    
    #Load step 2 h2 data
    h2_res <- read.csv("Step2_Results/protein_h2_results.csv", header=T)
    
    #Prepare step 2 h2 data
    h2_res <- h2_res %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything())
    h2_res <- h2_res %>% select("Protein Short Code", "h2_no_se", "h2_se")
    names(h2_res) <- c("Protein", "h2", "h2_se")
    #h2_res <- h2_res %>% arrange(desc(h2))    

    #Load step 2 prs data
    
    #individual-sample-prs
    step2_prs_all_with_apoe_res <- read.csv("Step2_Results/protein_prs_all_with_apoe_indiv_sample_res.csv", header=T)
    step2_prs_all_no_apoe_res <- read.csv("Step2_Results/protein_prs_all_no_apoe_indiv_sample_res.csv", header=T)    
    
    #meta-analysis-prs
    step2_prs_all_with_apoe_meta <- read.csv("Step2_Results/protein_prs_all_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_all_no_apoe_meta <- read.csv("Step2_Results/protein_prs_all_no_apoe_meta_analysis_res.csv", header=T)    
    
    #Prepare step 2 prs data
    
    #individual-sample-prs
    step2_prs_all_with_apoe_res <- step2_prs_all_with_apoe_res %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
    names(step2_prs_all_with_apoe_res)[1] <- "Protein"
    step2_prs_all_with_apoe_res <- step2_prs_all_with_apoe_res %>% mutate("Protein_Sample" = paste(Protein, Sample))
    
    step2_prs_all_no_apoe_res <- step2_prs_all_no_apoe_res %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
    names(step2_prs_all_no_apoe_res)[1] <- "Protein"
    step2_prs_all_no_apoe_res <- step2_prs_all_no_apoe_res %>% mutate("Protein_Sample" = paste(Protein, Sample))
    
    #meta-analysis-prs
    
    step2_prs_all_with_apoe_meta  <- step2_prs_all_with_apoe_meta  %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
    names(step2_prs_all_with_apoe_meta)[1] <- "Protein"
    step2_prs_all_with_apoe_meta_min_p <- step2_prs_all_with_apoe_meta %>% group_by(Protein) %>% filter(p == min(p))
    
    step2_prs_all_no_apoe_meta  <- step2_prs_all_no_apoe_meta  %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
    names(step2_prs_all_no_apoe_meta)[1] <- "Protein"
    step2_prs_all_no_apoe_meta_min_p <- step2_prs_all_no_apoe_meta %>% group_by(Protein) %>% filter(p == min(p))
    
    #OUTPUT#
    
    #consider adding legend (few attempts inc. geom_text and adding aes layer to hline. Conclusion is may require chart redesign)
    output$h2 <- renderPlot({
        ggplot(data = h2_res) + geom_pointrange(mapping = aes(x=reorder(Protein,h2), y=h2, ymin=h2-h2_se, ymax=h2+h2_se)) + geom_hline(yintercept=1, linetype="dashed", color = "red") + geom_hline(yintercept=mean(h2_res$h2), linetype="dashed", color = "blue") + geom_hline(yintercept=0, linetype="dashed", color = "red") + xlab("Protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0))
    })
    
    output$prs_indiv <- renderPlot({
        
        if(input$apoe == "With-APOE" & input$group == "All"){
            ggplot(step2_prs_all_with_apoe_res, aes(Protein, P_MinusLog10, fill=Sample)) + geom_bar(stat = "identity", position=position_dodge()) + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "red") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("-log10 p-value") + xlab("Protein") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red"))))
        } else if (input$apoe == "No-APOE" & input$group == "All") {
            ggplot(step2_prs_all_no_apoe_res, aes(Protein, P_MinusLog10, fill=Sample)) + geom_bar(stat = "identity", position=position_dodge()) + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "red") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("-log10 p-value") + xlab("Protein") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red"))))
        } else {
        }
    })
    
    output$prs_meta <- renderPlot({
        if(input$apoe == "With-APOE" & input$group == "All") {
            ggplot(step2_prs_all_with_apoe_meta_min_p, aes(Protein, P_MinusLog10, label=Threshold, fill=Significant)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "red") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="blue", "N"="gray" ), guide = FALSE )
        } else if (input$apoe == "No-APOE" & input$group == "All") {
            ggplot(step2_prs_all_no_apoe_meta_min_p, aes(Protein, P_MinusLog10, label=Threshold, fill=Significant)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "red") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("red")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="blue", "N"="gray" ), guide = FALSE )
        } else {
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
