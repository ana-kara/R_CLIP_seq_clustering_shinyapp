
library(shiny)
options(shiny.maxRequestSize=10000*1024^2)

ui <- fluidPage(

    
    titlePanel("CLIP shiny app"),

    
    sidebarLayout(
        sidebarPanel(
            helpText("Input of atleast 1 BED file is required, while input of correspondng BAM file is optional"),
            
            fileInput("bed_file", p("BED input"),multiple = TRUE),
            
            fileInput("bam_file", p("BAM input (required for jaccard & spearman dissimilarity matrices)"),multiple = TRUE, accept = ".bam"),
            
            radioButtons("checkGroup", 
                               h3("Checkbox group"), 
                               choices = list("Jaccard distance" = "jaccard", 
                                              "Bray distance" = "bray", 
                                              "Spearman distance" = "spearman"),
                               selected = 1),
            actionButton(inputId = "click",label = "Run")),
        mainPanel(textOutput("erica"),textOutput("text_warning"),plotOutput("pam_graph"),plotOutput("hist_graph")) ))




 server <- function(input, output) {
     
     
     
     library(rlist)
     library(factoextra)
     library(cluster)
     library(vegan)
     library(tibble)
     library(ggpubr)
     library(varhandle)
     library(tidyverse)
     library(valr)
     library(lobstr)
     
     
     cat_list <- list.load("chr_start_end.rds")
     
     bb <- read.delim("all_sorted_merged_peaks", header = F)
     bb$V1 <- unfactor(bb$V1)
     colnames(bb) <- c("chrom", "start", "end","strand")
     
     
     
     
     gr1 <-   as.data.frame(bb) %>% 
         select(chrom, start, end, strand) 
     
     
     
     
     
     les <- eventReactive(input$click,{
         
         
         
         #cat_list <- list.load("chr_start_end.rds")
         
         #bb <- read.delim("all_sorted_merged_peaks", header = F)
         #bb$V1 <- unfactor(bb$V1)
         #colnames(bb) <- c("chrom", "start", "end","strand")
         
         
         
         
         #gr1 <-   as.data.frame(bb) %>% 
             #select(chrom, start, end, strand) 
         #rm(bb)
         
         temp_files <- list()
         for(x in 1:length(input$bed_file[,"datapath"])){
             
             aa <- read.delim(input$bed_file[x,"datapath"], header = F)
             aa$V1 <- unfactor(aa$V1)
             names(aa)[names(aa) == 'V1'] <- "chrom"
             names(aa)[names(aa) == 'V2'] <- "start"
             names(aa)[names(aa) == 'V3'] <- "end"
             names(aa)[names(aa) == 'V6'] <- "strand"
             gr2 <- as.data.frame(aa) %>% 
                 select(chrom, start, end, strand) 
             
             abc <- valr::bed_intersect(gr2, gr1, suffix = c("_gr2", "_gr1"))
             abc_2 <- abc[which(abc[,4] == abc[,7]),]
             
             cat_list_indeces <- sapply(1:nrow(abc_2), function(x) which(cat_list %in% paste(abc_2[x,1],abc_2[x,5],abc_2[x,6])))
             temp_files[[input$bed_file[x,"name"]]] <- cat_list_indeces
         }
         rm(cat_list,gr1,aa,gr2,cat_list_indeces)
         return(temp_files)}) 
     
     data1 <- eventReactive(input$click,{
         
         
         
         
         ordered_protein_names <- list.load("ordered_protein_names.rds")
         bam_present <- length(input$bam_file)
         
         
         all_indexes_2 <- list.load("all_indexes_HepG2_K562.rds")
         all_indexes_3 <- c(all_indexes_2,les())
         rm(all_indexes_2)
         
         
         
         final_K562_HepG2_table <- matrix(rep(0, 250713), nrow =  250713, ncol = length(all_indexes_3))
         colnames(final_K562_HepG2_table) <- c(ordered_protein_names,input$bed_file[,"name"])
         
         
         if(bam_present > 0){
             library(GenomicAlignments)
             library(GenomicFeatures)
             library(BiocParallel)
             library(Rsamtools)
             
             
             
             bam_dir_HepG2 <- input$bam_file[,"datapath"]
             merged_file_granges_HepG2_K562 <- list.load("HepG2_K562_granges_list.rds")
             bamfiles_K562_HepG2_both <- BamFileList(bam_dir_HepG2, yieldSize=2000000)
             
             
             
             register(MulticoreParam(2))
             
             se_K562_HepG2_both_bams_2 <- summarizeOverlaps(features=merged_file_granges_HepG2_K562,
                                                            reads=bamfiles_K562_HepG2_both,
                                                            mode="Union",
                                                            singleEnd=FALSE,
                                                            ignore.strand=FALSE,
                                                            fragments=TRUE )
             
             
             rm(bam_dir_HepG2,merged_file_granges_HepG2_K562,bamfiles_K562_HepG2_both)
             se_K562_HepG2_both_bams <- list.load("se_K562_HepG2_idr_merged_peaks_both_bams.rds")
             unordered_protein_names <- list.load("unordered_protein_names.rds")
             colnames(se_K562_HepG2_both_bams) <- unordered_protein_names
             rm(unordered_protein_names)
             
             se_K562_HepG2_both_bams <- se_K562_HepG2_both_bams[,ordered_protein_names]
             colnames(se_K562_HepG2_both_bams_2) <- input$bed_file[,"name"]
             se_K562_HepG2_both_bams_3 <- cbind(se_K562_HepG2_both_bams,se_K562_HepG2_both_bams_2)
             
             
             library(DESeq2)
             dds_K562_HepG2_both_bams <- DESeqDataSet(se_K562_HepG2_both_bams_3, design = ~ 1) # _3
             dds_K562_HepG2_both_bams <- DESeq2::estimateSizeFactors(dds_K562_HepG2_both_bams)
             rm(se_K562_HepG2_both_bams_3)
             normalized_counts_K562_HepG2_both_bams <- DESeq2::counts(dds_K562_HepG2_both_bams, normalized = TRUE)
             rm(dds_K562_HepG2_both_bams)
             
             for (index_item in 1:NROW(all_indexes_3)) {
                 last_counts <- round(normalized_counts_K562_HepG2_both_bams[all_indexes_3[[index_item]],index_item])
                 last_counts[last_counts == 0] <- 1
                 final_K562_HepG2_table[all_indexes_3[[index_item]],index_item] <-  last_counts
             }
             
         } else{
             for (index_item in 1:NROW(all_indexes_3)) {
                 final_K562_HepG2_table[all_indexes_3[[index_item]],index_item] <-  1
             }
         }
         rm(ordered_protein_names,all_indexes_3)
         final_K562_HepG2_table
     })
     
     
     
     output$pam_graph <- renderPlot({
         req(input$checkGroup)
         if(input$checkGroup == "jaccard"){
             
             dist_matrix <- as.matrix(vegdist(t(data1()), method= input$checkGroup,binary = T))
             distance_name <- "jaccard distance"
             
         } else if(input$checkGroup == "bray" & length(input$bam_file) > 0){
             
             dist_matrix <- as.matrix(vegdist(t(data1()), method= input$checkGroup))
             distance_name <- "bray distance"
             
         } else if(input$checkGroup == "spearman" & length(input$bam_file) > 0){
             
             dist_matrix <- 1 - cor(data1(), method= input$checkGroup)
             distance_name <- "spearman distance"
             
         } else{
             
             dist_matrix <- as.matrix(vegdist(t(data1()), method= "jaccard",binary = T)) 
             output$text_warning <- renderText({ "No BAM file(s) were provided. The default Jaccard distance was used instead" })
             distance_name <- "jaccard distance"
         }
         
         
         
         pam_result <- pam(dist_matrix, 4, diss = TRUE)
         
         
         
         
         
         mds.cor <- dist_matrix %>%
             cmdscale() %>%
             as_tibble()
         colnames(mds.cor) <- c("Dim.1", "Dim.2")
         mds.cor <- mds.cor %>%
             mutate(groups = as.factor(pam_result$cluster))
         
         
         
         ggscatter(mds.cor, x = "Dim.1", y = "Dim.2",
                   label = names(pam_result$cluster),
                   color = "groups",
                   title = paste0("Pam clustering of proteins of interest with ",distance_name),
                   palette = "jco",
                   size = 1,
                   ellipse = TRUE,
                   ellipse.type = "convex",
                   repel = TRUE)
         
     })
     output$hist_graph <- renderPlot({
         req(input$checkGroup)
         if(input$checkGroup == "jaccard"){
             
             hc1 <- hclust(vegdist(t(data1()), method= input$checkGroup,binary = T), method = "ward.D" )
             distance_name <- "jaccard distance"
             
         } else if(input$checkGroup == "bray" & length(input$bam_file) > 0){
             
             hc1 <- hclust(vegdist(t(data1()), method= input$checkGroup), method = "ward.D")
             distance_name <- "bray distance"
             
         } else if(input$checkGroup == "spearman" & length(input$bam_file) > 0){
             
             hc1 <- hclust(dist(1 - cor(data1(), method= input$checkGroup)), method = "ward.D")
             distance_name <- "spearman distance"
             
         } else{
             
             hc1 <- hclust(vegdist(t(data1()), method= "jaccard",binary = T), method = "ward.D" )
             distance_name <- "jaccard distance"
         }
         
         
         plot(hc1, cex = 0.6, hang = -1,main = paste0("Hclust clustering of proteins of interest with ",distance_name))
        
     })
     
     output$erica <- renderText({
         #req(data1())
         mem_used()})
 }
 
 
 
 
 shinyApp(ui = ui, server = server)
 
