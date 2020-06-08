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

select_data <- function(apoe, group, results){
    if (apoe == "With-APOE" & group == "All"){ results[[1]] }
    else if (apoe == "No-APOE" & group == "All") { results[[2]] }
    else if (apoe == "With-APOE" & group == "Males") { results[[3]] }
    else if (apoe == "No-APOE" & group == "Males") { results[[4]] }
    else if (apoe == "With-APOE" & group == "Females") { results[[5]] }
    else if (apoe == "No-APOE" & group == "Females") { results[[6]] }
    else if (apoe == "With-APOE" & group == "65.And.Over") { results[[7]] }
    else if (apoe == "No-APOE" & group == "65.And.Over") { results[[8]] }
    else if (apoe == "With-APOE" & group == "70.And.Over") { results[[9]] }
    else if (apoe == "No-APOE" & group == "70.And.Over") { results[[10]] }
    else if (apoe == "With-APOE" & group == "75.And.Over") { results[[11]] }
    else if (apoe == "No-APOE" & group == "75.And.Over") { results[[12]] }
    else if (apoe == "With-APOE" & group == "80.And.Over") { results[[13]] }
    else if (apoe == "No-APOE" & group == "80.And.Over") { results[[14]] }
    else { }
}

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
                condition = "input.step == '2' | input.step == '3'",
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
               condition = "input.step == '2' & input.chart == 'individual-sample-prs'",
               tableOutput("prs_indiv_table")
           ),
           
           
           conditionalPanel(
               condition = "input.step == '2' & input.chart == 'meta-analysis-prs'",
               plotOutput("prs_meta")
           ),
           
           conditionalPanel(
               condition = "input.step == '2' & input.chart == 'meta-analysis-prs'",
               tableOutput("prs_meta_table")
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

    #Load step 2 prs data
    
    #individual-sample-prs
    step2_prs_all_with_apoe_res <- read.csv("Step2_Results/protein_prs_all_with_apoe_indiv_sample_res.csv", header=T)
    step2_prs_all_no_apoe_res <- read.csv("Step2_Results/protein_prs_all_no_apoe_indiv_sample_res.csv", header=T)    
   
    step2_prs_males_with_apoe_res <- read.csv("Step2_Results/protein_prs_Males_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_males_no_apoe_res <- read.csv("Step2_Results/protein_prs_Males_no_apoe_indiv_sample_res.csv", header=T)    
    
    step2_prs_females_with_apoe_res <- read.csv("Step2_Results/protein_prs_Females_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_females_no_apoe_res <- read.csv("Step2_Results/protein_prs_Females_no_apoe_indiv_sample_res.csv", header=T)    
    
    step2_prs_65_with_apoe_res <- read.csv("Step2_Results/protein_prs_65.And.Over_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_65_no_apoe_res <- read.csv("Step2_Results/protein_prs_65.And.Over_no_apoe_indiv_sample_res.csv", header=T)    
    
    step2_prs_70_with_apoe_res <- read.csv("Step2_Results/protein_prs_70.And.Over_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_70_no_apoe_res <- read.csv("Step2_Results/protein_prs_70.And.Over_no_apoe_indiv_sample_res.csv", header=T)    
    
    step2_prs_75_with_apoe_res <- read.csv("Step2_Results/protein_prs_75.And.Over_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_75_no_apoe_res <- read.csv("Step2_Results/protein_prs_75.And.Over_no_apoe_indiv_sample_res.csv", header=T)  
    
    step2_prs_80_with_apoe_res <- read.csv("Step2_Results/protein_prs_80.And.Over_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_80_no_apoe_res <- read.csv("Step2_Results/protein_prs_80.And.Over_no_apoe_indiv_sample_res.csv", header=T)
    
    step2_prs_res_indiv <- list(step2_prs_all_with_apoe_res, 
                      step2_prs_all_no_apoe_res, 
                      step2_prs_males_with_apoe_res, 
                      step2_prs_males_no_apoe_res, 
                      step2_prs_females_with_apoe_res, 
                      step2_prs_females_no_apoe_res,
                      step2_prs_65_with_apoe_res,
                      step2_prs_65_no_apoe_res,
                      step2_prs_70_with_apoe_res,
                      step2_prs_70_no_apoe_res,
                      step2_prs_75_with_apoe_res,
                      step2_prs_75_no_apoe_res,
                      step2_prs_80_with_apoe_res,
                      step2_prs_80_no_apoe_res)

    
    #meta-analysis-prs
    step2_prs_all_with_apoe_meta <- read.csv("Step2_Results/protein_prs_all_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_all_no_apoe_meta <- read.csv("Step2_Results/protein_prs_all_no_apoe_meta_analysis_res.csv", header=T)    
    
    step2_prs_males_with_apoe_meta <- read.csv("Step2_Results/protein_prs_Males_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_males_no_apoe_meta <- read.csv("Step2_Results/protein_prs_Males_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_females_with_apoe_meta <- read.csv("Step2_Results/protein_prs_Females_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_females_no_apoe_meta <- read.csv("Step2_Results/protein_prs_Females_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_65_with_apoe_meta <- read.csv("Step2_Results/protein_prs_65.And.Over_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_65_no_apoe_meta <- read.csv("Step2_Results/protein_prs_65.And.Over_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_70_with_apoe_meta <- read.csv("Step2_Results/protein_prs_70.And.Over_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_70_no_apoe_meta <- read.csv("Step2_Results/protein_prs_70.And.Over_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_75_with_apoe_meta <- read.csv("Step2_Results/protein_prs_75.And.Over_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_75_no_apoe_meta <- read.csv("Step2_Results/protein_prs_75.And.Over_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_80_with_apoe_meta <- read.csv("Step2_Results/protein_prs_80.And.Over_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_80_no_apoe_meta <- read.csv("Step2_Results/protein_prs_80.And.Over_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_res_meta <- list(step2_prs_all_with_apoe_meta,
                               step2_prs_all_no_apoe_meta, 
                      step2_prs_males_with_apoe_meta, 
                      step2_prs_males_no_apoe_meta, 
                      step2_prs_females_with_apoe_meta, 
                      step2_prs_females_no_apoe_meta,
                      step2_prs_65_with_apoe_meta,
                      step2_prs_65_no_apoe_meta,
                      step2_prs_70_with_apoe_meta,
                      step2_prs_70_no_apoe_meta,
                      step2_prs_75_with_apoe_meta,
                      step2_prs_75_no_apoe_meta,
                      step2_prs_80_with_apoe_meta,
                      step2_prs_80_no_apoe_meta)
    
    #OUTPUT#
    
    #consider adding legend (few attempts inc. geom_text and adding aes layer to hline. Conclusion is may require chart redesign)
    output$h2 <- renderPlot({
        ggplot(data = h2_res) + geom_pointrange(mapping = aes(x=reorder(Protein,h2), y=h2, ymin=h2-h2_se, ymax=h2+h2_se)) + geom_hline(yintercept=1, linetype="dashed", color = "red") + geom_hline(yintercept=mean(h2_res$h2), linetype="dashed", color = "blue") + geom_hline(yintercept=0, linetype="dashed", color = "red") + xlab("Protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0))
    })
    
    output$prs_indiv <- renderPlot({
        
        apoe <- input$apoe
        group <- input$group
        
        chart_data <- select_data(apoe, group, step2_prs_res_indiv)
        
        #issue with sample data labels for groups (review results preparation script)
        chart_data <- chart_data %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
        names(chart_data)[1] <- "Protein"
        chart_data <- chart_data %>% mutate("Protein_Sample" = paste(Protein, Sample))
        
        ggplot(chart_data, aes(Protein, P_MinusLog10, fill=Sample)) + geom_bar(stat = "identity", position=position_dodge()) + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "blue") + geom_hline(aes(yintercept=2.792392, linetype="Bonferroni corrected"), color = "green") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("-log10 p-value") + xlab("Protein") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("green", "blue"))))
        
    })
    
    output$prs_indiv_table <- renderTable({ 
        
        apoe <- input$apoe
        group <- input$group
        
        table_data <- select_data(apoe, group, step2_prs_res_indiv)
 
        table <- table_data %>% filter(Significant == "Y") %>% select(Protein, Sample, Threshold, P) %>% mutate_if(is.numeric, format, scientific=T, digits=1)
    
    })
    
    output$prs_meta <- renderPlot({
        
        apoe <- input$apoe
        group <- input$group
        
        chart_data <- select_data(apoe, group, step2_prs_res_meta)
        
        chart_data  <- chart_data  %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
        names(chart_data)[1] <- "Protein"
        chart_data <- chart_data %>% group_by(Protein) %>% filter(p == min(p))
        

        ggplot(chart_data, aes(Protein, P_MinusLog10, label=Threshold, fill=Significant)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "blue") + geom_hline(aes(yintercept=2.792392, linetype="Bonferroni corrected"), color = "green") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("green", "blue")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="blue", "N"="gray" ), guide = FALSE )
        
    })
    
    output$prs_meta_table <- renderTable({ 
        
        apoe <- input$apoe
        group <- input$group
        
        table_data <- select_data(apoe, group, step2_prs_res_meta)
        
        table_data  <- table_data  %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
        names(table_data)[1] <- "Protein"
        table_data <- table_data %>% group_by(Protein) %>% filter(p == min(p))
        
        table <- table_data %>% filter(Significant == "Y") %>% select(Protein, p) %>% mutate_if(is.numeric, format, scientific=T, digits=1)
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
