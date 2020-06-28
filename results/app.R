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
library(DT)

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
            selectInput("intro", strong("Show me: "), 
                        choices = c("Project background","Results"), selected = "Project background"),
            
            conditionalPanel(
                condition = "input.intro == 'Results'",
                conditionalPanel(
                    condition = "input.intro == 'Results'",
                    selectInput("step", strong("Show results for: "), 
                            choices = c("Step 1: Protein shortlist","Step 2: Protein PRS to AD","Step 3: AD PRS to protein", "Step 4: Bi-directional MR"), selected = "Step 1: Protein shortlist"),
                ),
                
                conditionalPanel(
                    condition = "input.step == 'Step 2: Protein PRS to AD'",
                    selectInput("chart", strong("Analysis"), 
                                choices = c("Heritability","Individual sample PRS","Meta-analysis PRS"), selected = "Heritability"),
                ), 
                
                conditionalPanel(
                    condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart != 'Heritability' | input.step == 'Step 3: AD PRS to protein'",
                    selectInput("apoe", strong("Include APOE SNPs in PRS"), 
                                choices = c("Yes","No"), selected = "Yes"),
                ), 
                
                conditionalPanel(
                    condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart != 'Heritability' | input.step == 'Step 3: AD PRS to protein'",
                    selectInput("group", strong("Sample group"), 
                                choices = c("All","Males", "Females", "65 and over", "70 and over", "75 and over", "80 and over"), selected = "All"),
                ), 
                
                conditionalPanel(
                    condition = "input.step == 'Step 2: Protein PRS to AD' & input.chart == 'Meta-analysis PRS' | input.step == 'Step 3: AD PRS to protein'",
                    selectInput("adjp", strong("Adjust p-values with BH"), 
                                choices = c("Yes","No"), selected = "No"),
                ), 
                
                conditionalPanel(
                    condition = "input.step == 'Step 4: Bi-directional MR'",
                    selectInput("exposure", strong("Exposure"), 
                                choices = c("AD","Proteins"), selected = "Proteins"),
                ),
                
                conditionalPanel(
                    condition = "input.step == 'Step 4: Bi-directional MR' & input.exposure == 'Proteins'",
                    selectInput("threshold", strong("P-value threshold for protein SNP instruments"), 
                                choices = c("5e-08","5e-06"), selected = "5e-08"),
                ),
                
                conditionalPanel(
                    condition = "input.step == 'Step 4: Bi-directional MR'",
                    selectInput("outcome", strong("Outcome"), 
                                choices = c("AD","Proteins"), selected = "AD"),
                ),
                
                conditionalPanel(
                    condition = "input.step == 'Step 4: Bi-directional MR'",
                    selectInput("harmonisation", strong("Harmonisation assumption"), 
                                choices = c("All alleles forward strand","Palindromic SNPs inferred or excluded"), selected = "All alleles forward strand"),
                )
            )    
            
        ),

        mainPanel(
            conditionalPanel(
                condition = "input.intro == 'Project background'",
                htmlOutput("background")
            ),
            
            conditionalPanel(
                condition = "input.intro == 'Results'",
                conditionalPanel(
                    condition = "input.step == 'Step 1: Protein shortlist'",
                    DT::dataTableOutput('proteins')
                ),
                
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
               ),
               
               conditionalPanel(
                   condition = "input.step == 'Step 3: AD PRS to protein'",
                   plotOutput("step3_prs")
               ), 
               
               conditionalPanel(
                   condition = "input.step == 'Step 4: Bi-directional MR'",
                   DT::dataTableOutput("step4_mr_table")
               ), 
               
               conditionalPanel(
                   condition = "input.step == 'Step 4: Bi-directional MR' & input.exposure == 'AD' & input.outcome == 'AD' | input.exposure == 'Proteins' & input.outcome == 'Proteins'",
                   textOutput("step4_mr_error")
               )
            )
        )
    )
)

server <- function(input, output) {
    
    #DATA PREPARATION#
    
    #Load step 1 protein shortlist
    
    prots <- read.csv("step1_results/protein_shortlist.csv", header=T)
    
    #Load step 2 h2 data
    h2_res <- read.csv("step2_results/protein_h2_results.csv", header=T)
    
    #Prepare step 2 h2 data
    h2_res <- h2_res %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything())
    h2_res <- h2_res %>% select("Protein Short Code", "h2_no_se", "h2_se")
    names(h2_res) <- c("Protein", "h2", "h2_se")

    #Load step 2 prs data
    
    #Step 2 prs individual sample load
    
    step2_prs_all_with_apoe_res <- read.csv("step2_results/protein_prs_all_with_apoe_indiv_sample_res.csv", header=T)
    step2_prs_all_no_apoe_res <- read.csv("step2_results/protein_prs_all_no_apoe_indiv_sample_res.csv", header=T)    
   
    step2_prs_males_with_apoe_res <- read.csv("step2_results/protein_prs_Males_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_males_no_apoe_res <- read.csv("step2_results/protein_prs_Males_no_apoe_indiv_sample_res.csv", header=T)    
    
    step2_prs_females_with_apoe_res <- read.csv("step2_results/protein_prs_Females_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_females_no_apoe_res <- read.csv("step2_results/protein_prs_Females_no_apoe_indiv_sample_res.csv", header=T)    
    
    step2_prs_65_with_apoe_res <- read.csv("step2_results/protein_prs_65.And.Over_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_65_no_apoe_res <- read.csv("step2_results/protein_prs_65.And.Over_no_apoe_indiv_sample_res.csv", header=T)    
    
    step2_prs_70_with_apoe_res <- read.csv("step2_results/protein_prs_70.And.Over_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_70_no_apoe_res <- read.csv("step2_results/protein_prs_70.And.Over_no_apoe_indiv_sample_res.csv", header=T)    
    
    step2_prs_75_with_apoe_res <- read.csv("step2_results/protein_prs_75.And.Over_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_75_no_apoe_res <- read.csv("step2_results/protein_prs_75.And.Over_no_apoe_indiv_sample_res.csv", header=T)  
    
    step2_prs_80_with_apoe_res <- read.csv("step2_results/protein_prs_80.And.Over_with_apoe_indiv_sample_res.csv", header=T)    
    step2_prs_80_no_apoe_res <- read.csv("step2_results/protein_prs_80.And.Over_no_apoe_indiv_sample_res.csv", header=T)
    
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

    
    #Step 2 prs meta-analysis load
    
    step2_prs_all_with_apoe_meta <- read.csv("step2_results/protein_prs_all_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_all_no_apoe_meta <- read.csv("step2_results/protein_prs_all_no_apoe_meta_analysis_res.csv", header=T)    
    
    step2_prs_males_with_apoe_meta <- read.csv("step2_results/protein_prs_Males_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_males_no_apoe_meta <- read.csv("step2_results/protein_prs_Males_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_females_with_apoe_meta <- read.csv("step2_results/protein_prs_Females_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_females_no_apoe_meta <- read.csv("step2_results/protein_prs_Females_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_65_with_apoe_meta <- read.csv("step2_results/protein_prs_65.And.Over_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_65_no_apoe_meta <- read.csv("step2_results/protein_prs_65.And.Over_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_70_with_apoe_meta <- read.csv("step2_results/protein_prs_70.And.Over_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_70_no_apoe_meta <- read.csv("step2_results/protein_prs_70.And.Over_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_75_with_apoe_meta <- read.csv("step2_results/protein_prs_75.And.Over_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_75_no_apoe_meta <- read.csv("step2_results/protein_prs_75.And.Over_no_apoe_meta_analysis_res.csv", header=T)
    
    step2_prs_80_with_apoe_meta <- read.csv("step2_results/protein_prs_80.And.Over_with_apoe_meta_analysis_res.csv", header=T)
    step2_prs_80_no_apoe_meta <- read.csv("step2_results/protein_prs_80.And.Over_no_apoe_meta_analysis_res.csv", header=T)
    
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
    
    
    #Step 3 prs results load
    
    step3_prs_all_with_apoe <- read.csv("step3_results/ad_prs_all_with_apoe_res.csv", header=T)
    step3_prs_all_no_apoe<- read.csv("step3_results/ad_prs_all_no_apoe_res.csv", header=T)    
    
    step3_prs_males_with_apoe <- read.csv("step3_results/ad_prs_males_only_with_apoe_res.csv", header=T)
    step3_prs_males_no_apoe<- read.csv("step3_results/ad_prs_males_only_no_apoe_res.csv", header=T)
    
    step3_prs_females_with_apoe <- read.csv("step3_results/ad_prs_females_only_with_apoe_res.csv", header=T)
    step3_prs_females_no_apoe <- read.csv("step3_results/ad_prs_females_only_no_apoe_res.csv", header=T)
    
    step3_prs_65_with_apoe <- read.csv("step3_results/ad_prs_65_and_over_with_apoe_res.csv", header=T)
    step3_prs_65_no_apoe <- read.csv("step3_results/ad_prs_65_and_over_no_apoe_res.csv", header=T)
    
    step3_prs_70_with_apoe <- read.csv("step3_results/ad_prs_70_and_over_with_apoe_res.csv", header=T)
    step3_prs_70_no_apoe <- read.csv("step3_results/ad_prs_70_and_over_no_apoe_res.csv", header=T)
    
    step3_prs_75_with_apoe <- read.csv("step3_results/ad_prs_75_and_over_with_apoe_res.csv", header=T)
    step3_prs_75_no_apoe <- read.csv("step3_results/ad_prs_75_and_over_no_apoe_res.csv", header=T)
    
    step3_prs_80_with_apoe <- read.csv("step3_results/ad_prs_80_and_over_with_apoe_res.csv", header=T)
    step3_prs_80_no_apoe <- read.csv("step3_results/ad_prs_80_and_over_no_apoe_res.csv", header=T)
    
    step3_prs_res <- list(step3_prs_all_with_apoe,
                             step3_prs_all_no_apoe, 
                             step3_prs_males_with_apoe, 
                             step3_prs_males_no_apoe, 
                             step3_prs_females_with_apoe, 
                             step3_prs_females_no_apoe,
                             step3_prs_65_with_apoe,
                             step3_prs_65_no_apoe,
                             step3_prs_70_with_apoe,
                             step3_prs_70_no_apoe,
                             step3_prs_75_with_apoe,
                             step3_prs_75_no_apoe,
                             step3_prs_80_with_apoe,
                             step3_prs_80_no_apoe)
    
    
    #Step 4 mr results load
    
    step4_mr_protein_to_ad_5e08_1 <- read.csv("step4_results/protein_to_ad_mr_5e-08_1.csv", header=T)
    step4_mr_protein_to_ad_5e08_2 <- read.csv("step4_results/protein_to_ad_mr_5e-08_2.csv", header=T)
    
    step4_mr_protein_to_ad_5e06_1 <- read.csv("step4_results/protein_to_ad_mr_5e-06_1.csv", header=T)
    step4_mr_protein_to_ad_5e06_2 <- read.csv("step4_results/protein_to_ad_mr_5e-06_2.csv", header=T)
    
    step4_mr_ad_to_protein_1 <- read.csv("step4_results/ad_to_protein_mr_1.csv", header=T)
    step4_mr_ad_to_protein_2 <- read.csv("step4_results/ad_to_protein_mr_2.csv", header=T)
    
    #OUTPUT#
    
    getBackground<-function() {
        return(includeHTML("background.html"))
    }
    output$background <- renderUI({getBackground()})
    
    output$proteins <- DT::renderDataTable(
        DT::datatable(prots,caption = 'From 175 candidate proteins identified in the literature review, 31 were selected for PRS analysis based on replicability and genetic data availability.', options = list(pageLength=10,columnDefs = list(list(className = 'dt-center', targets = c(1:2)))), rownames = F,colnames = c('Protein', 'Short Code', '# of Studies Replicated In'))
    )
    
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
        
        #bonf threshold -> -log10(0.00017) -> 3.769551
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
           #bonf threshold -> -log10(0.00017) -> 3.769551
           chart_data <- chart_data %>% mutate(SignificantBonf = if_else(p < 0.00017, "Y","N")) 
           chart_data <- chart_data %>% group_by(Protein) %>% filter(p == min(p))  
           ggplot(chart_data, aes(Protein, P_MinusLog10, label=Threshold, fill=SignificantBonf)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "blue") + geom_hline(aes(yintercept=3.769551, linetype="Bonferroni corrected"), color = "green") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("green", "blue")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="green4", "N"="gray" ), guide = FALSE) 
        }
    })
    
    output$step3_prs <- renderPlot({
        
        somamer_ids <- c("ANGPT2.13660.76.3","APOB.2797.56.2","APOE.2937.10.2","NPPB.3723.1.2","SERPING1.4479.14.2", "C3.2755.8.2","C4A.C4B.4481.34.2","CLU.4542.24.2","CRP.4337.49.2","FGA.FGB.FGG.4907.56.1","CFH.4159.130.1","CSF3.8952.65.3","HP.3054.3.2","IGFBP2.2570.72.5","IL10.2773.50.2","IL13.3072.4.2","IL3.4717.55.2","CXCL8.3447.64.2","MMP9.2579.17.5","PLG.3710.49.2","RETN.3046.31.1","APCS.2474.54.5","TNC.4155.3.2","TNF.5936.53.3","TF.4162.54.2","VCAM1.2967.8.1")
        
        apoe <- input$apoe
        group <- input$group
        adjusted_p <- input$adjp
        
        chart_data <- select_data(apoe, group, step3_prs_res)
        chart_data <- chart_data %>% group_by(Pheno) %>% filter(R2 == max(R2))
        chart_data$Protein <- somamer_ids
        chart_data  <- chart_data  %>% mutate("Protein Short Code" = gsub("\\..*","",Protein)) %>% select(1, `Protein Short Code`, everything()) %>% select(-Protein)
        #update index to 2, so keep original phenotype but show protein short code
        names(chart_data)[2] <- "Protein"
        
        #Present with and without adjusted p values
        if (adjusted_p == "Yes") {
            ggplot(chart_data, aes(Protein, adjp_MinusLog10, label=Threshold, fill=SignificantAdjP)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1, linetype="FDR < 0.1"), color = "green") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("green4")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="green4", "N"="gray" ), guide = FALSE )
        } else {
            #bonf threshold -> -log10(0.00019) -> 3.721246
            ggplot(chart_data, aes(Protein, P_MinusLog10, label=Threshold, fill=SignificantBonf)) + geom_bar(stat = "identity") + geom_hline(aes(yintercept=1.3, linetype="nominal p < 0.05"), color = "blue") + geom_hline(aes(yintercept=3.721246, linetype="Bonferroni corrected"), color = "green") + expand_limits(y = c(0, 2.5)) + geom_text(vjust=-0.5, size=2) + ylab("-log10 p-value") + scale_linetype_manual(name = "Significance", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("green", "blue")))) + labs(caption = "P-value threshold for most significant PRS model displayed above bar for each protein") + theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.caption = element_text(hjust = 0, face= "italic")) + scale_fill_manual( values = c( "Y"="green4", "N"="gray" ), guide = FALSE) 
        }
    })
    
    output$step4_mr_table <- 
        DT::renderDataTable({
            exposure <- input$exposure
            outcome <- input$outcome
            threshold <- input$threshold
            harmonisation <- input$harmonisation
            caption_text <- "6 proteins (Apolipoprotein E (isoform ε3), Apolipoprotein B-100, C-reactive protein, Vitamin D-binding protein, Insulin-like growth factor-binding protein 2 and Angiopoietin-2) were selected for MR analysis."
            headers <- c("Protein", "N SNPs", "MR Method", "OR", "95% CI Lower", "95% CI Upper", "P-value", "Egger Intercept p-value", "Cochran's Q p-value")
            headers_decimals <- c("OR", "95% CI Lower", "95% CI Upper", "P-value", "Egger Intercept p-value", "Cochran's Q p-value")
            
            #consider adding removal of egger p-value intercept values for non egger p-value methods
            if (exposure == "Proteins" & outcome == "AD") {
                if (threshold == "5e-08" & harmonisation == "All alleles forward strand"){
                    colnames(step4_mr_protein_to_ad_5e08_1) <- headers
                    DT::datatable(step4_mr_protein_to_ad_5e08_1, caption = caption_text, rownames = F, options = list(pageLength=10,columnDefs = list(list(className = 'dt-center', targets = c(1,4:8))))) %>% DT::formatRound(headers_decimals, 2)
                } else if (threshold == "5e-08" & harmonisation == "Palindromic SNPs inferred or excluded"){
                    colnames(step4_mr_protein_to_ad_5e08_2) <- headers
                    DT::datatable(step4_mr_protein_to_ad_5e08_2, caption = caption_text, rownames = F, options = list(pageLength=10,columnDefs = list(list(className = 'dt-center', targets = c(1,4:8))))) %>% DT::formatRound(headers_decimals, 2)
                } else if (threshold == "5e-06" & harmonisation == "All alleles forward strand") {
                    colnames(step4_mr_protein_to_ad_5e06_1) <- headers
                    DT::datatable(step4_mr_protein_to_ad_5e06_1, caption = caption_text, rownames = F, options = list(pageLength=10,columnDefs = list(list(className = 'dt-center', targets = c(1,4:8))))) %>% DT::formatRound(headers_decimals, 2)
                } else {
                    colnames(step4_mr_protein_to_ad_5e06_2) <- headers
                    DT::datatable(step4_mr_protein_to_ad_5e06_2, caption = caption_text, rownames = F, options = list(pageLength=10,columnDefs = list(list(className = 'dt-center', targets = c(1,4:8))))) %>% DT::formatRound(headers_decimals, 2)
                }
            } else if (exposure == "AD" & outcome == "Proteins") {
                if (harmonisation == "All alleles forward strand"){
                    colnames(step4_mr_ad_to_protein_1) <- headers
                    DT::datatable(step4_mr_ad_to_protein_1, caption = caption_text, rownames = F, options = list(pageLength=10,columnDefs = list(list(className = 'dt-center', targets = c(1,4:8))))) %>% DT::formatRound(headers_decimals, 2)
                } else {
                    colnames(step4_mr_ad_to_protein_2) <- headers
                    DT::datatable(step4_mr_ad_to_protein_2, caption = caption_text, rownames = F, options = list(pageLength=10,columnDefs = list(list(className = 'dt-center', targets = c(1,4:8))))) %>% DT::formatRound(headers_decimals, 2)
                }
            } else {
            }
    })
    
    output$step4_mr_error <- renderText({
        "Can't compare an exposure or outcome against itself, please select again"
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
