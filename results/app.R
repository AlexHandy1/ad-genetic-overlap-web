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
library(htmltools)

select_data <- function(apoe, group, results){
    if (apoe == "Yes" & group == "All"){ results[[1]] }
    else if (apoe == "No" & group == "All") { results[[2]] }
    else if (apoe == "Yes" & group == "Males") { results[[3]] }
    else if (apoe == "No" & group == "Males") { results[[4]] }
    else if (apoe == "Yes" & group == "Females") { results[[5]] }
    else if (apoe == "No" & group == "Females") { results[[6]] }
    else if (apoe == "Yes" & group == "65 and over") { results[[7]] }
    else if (apoe == "No" & group == "65 and over") { results[[8]] }
    else if (apoe == "Yes" & group == "70 and over") { results[[9]] }
    else if (apoe == "No" & group == "70 and over") { results[[10]] }
    else if (apoe == "Yes" & group == "75 and over") { results[[11]] }
    else if (apoe == "No" & group == "75 and over") { results[[12]] }
    else if (apoe == "Yes" & group == "80 and over") { results[[13]] }
    else if (apoe == "No" & group == "80 and over") { results[[14]] }
    else { }
}

ui <- fluidPage(
    includeCSS("custom.css"),
    
    titlePanel(""),

    sidebarLayout(
        sidebarPanel(
            selectInput("step", strong("Show results for: "), 
                        choices = c("Step 1: Protein shortlist","Step 2: Protein PRS to AD","Step 3: AD PRS to protein", "Step 4: Bi-directional MR"), selected = "Step 2: Protein PRS to AD"),
            
            conditionalPanel(
                condition = "input.step == 'Step 2: Protein PRS to AD'",
                selectInput("chart", strong("Analysis"), 
                            choices = c("Heritability","Individual sample PRS","Meta-analysis PRS"), selected = "Heritability"),
            ), 
            
            conditionalPanel(
                condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart != 'Heritability'",
                selectInput("apoe", strong("Include APOE SNPs in PRS"), 
                            choices = c("Yes","No"), selected = "Yes"),
            ), 
            
            conditionalPanel(
                condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart != 'Heritability'",
                selectInput("group", strong("Sample group"), 
                            choices = c("All","Males", "Females", "65 and over", "70 and over", "75 and over", "80 and over"), selected = "All"),
            ), 
            
            conditionalPanel(
                condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart == 'Meta-analysis PRS'",
                selectInput("adjp", strong("Adjust p-values with BH"), 
                            choices = c("Yes","No"), selected = "No"),
            )
        ),

        mainPanel(
           conditionalPanel(
               condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart == 'Heritability'",
               plotOutput("h2")
           ),
           
           conditionalPanel(
               condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart == 'Individual sample PRS'",
               plotOutput("prs_indiv")
           ),
           
           conditionalPanel(
               condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart == 'Meta-analysis PRS'",
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

    #Load step 2 prs data
    
    #Individual sample PRS
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

    
    #Meta-analysis PRS
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
        
        ggplot(chart_data, aes(Protein, P_MinusLog10, fill=Sample)) + geom_bar(stat = "identity", position=position_dodge()) + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "blue") + geom_hline(aes(yintercept=3.769551, linetype="Bonferroni corrected"), color = "green") + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("-log10 p-value") + xlab("Protein") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("green", "blue"))))
        
    })
    
    output$prs_meta <- renderPlot({
        
        apoe <- input$apoe
        group <- input$group
        adjusted_p <- input$adjp
        
        chart_data <- select_data(apoe, group, step2_prs_res_meta)
        
        
        chart_data  <- chart_data  %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
        names(chart_data)[1] <- "Protein"
        
        #Present with and without adjusted p values
        if (adjusted_p == "Yes") {
            chart_data <- chart_data %>% mutate(adjp = p.adjust(p, method="BH")) %>% mutate(adjp_MinusLog10 = -log10(adjp)) %>% mutate(SignificantAdjP = if_else(adjp < 0.1, "Y","N")) 
            chart_data <- chart_data %>% group_by(Protein) %>% filter(adjp == min(adjp)) %>% filter(p == min(p))
            ggplot(chart_data, aes(Protein, adjp_MinusLog10, label=Threshold, fill=SignificantAdjP)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1, linetype="FDR < 0.1"), color = "green") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("green4")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="green4", "N"="gray" ), guide = FALSE )
        } else {
           chart_data <- chart_data %>% mutate(SignificantBonf = if_else(p < 0.00017, "Y","N")) 
           chart_data <- chart_data %>% group_by(Protein) %>% filter(p == min(p))  
           ggplot(chart_data, aes(Protein, P_MinusLog10, label=Threshold, fill=SignificantBonf)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "blue") + geom_hline(aes(yintercept=3.769551, linetype="Bonferroni corrected"), color = "green") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("green", "blue")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="green4", "N"="gray" ), guide = FALSE) 
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
