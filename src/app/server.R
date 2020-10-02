shinyServer <- function(input, output, session) {
    
    output$formula1 <- renderUI({
        withMathJax(helpText("Probability density function: $$P(x)=\\frac{B(x+\\alpha,n-x+\\beta)\\binom{n}{x}}{B(\\alpha, \\beta)}$$"))
    })
    
    output$formula2 <- renderUI({
      withMathJax(helpText("Mean: $$\\mu=\\frac{n\\alpha}{\\alpha+\\beta}$$"))
    })
    
    output$formula3 <- renderUI({
      withMathJax(helpText("Variance: $$\\sigma^{2}=\\frac{n\\alpha\\beta(n+\\alpha+\\beta)}{(\\alpha+\\beta)^2(1+\\alpha+\\beta)}$$"))
    })
    
    
    output$formula4 <- renderUI({
      withMathJax(helpText("$$P(D|G)=\\frac{P(G|D)P(D)}{P(G)}$$"))
    })
    
    output$distPlot <- renderPlot({
        n <- 10
        x <- 6
        a <- 1
        b <- 1
        theta <- seq(0, 1, 0.001)
        
        prior <- dbeta(theta, a, b)  # dbeta(x, shape1, shape2, ncp = 0, log = FALSE) where x is vector of quantiles
        posterior <- dbeta(theta, a + x, b + n - x)  # X IS THE NEW INFORMATION FROM YOUR DATA FOR P(G|D)
        med <- qbeta(0.5, a + x, b + n - x)
        
        plot(theta, posterior, lty = 1, lwd = 1)
        
    })
    
    output$plot <- renderPlot({
      table <- fit_model()
     #   freq_table <- ggv(chr = fit_model()[1, 5], pos = fit_model()[1,6], db = "1000genomes", output = "table")
         freq_table <- ggv(chr=11, pos=47355475, db='1000genomes', output='table')
        envmap()
        addpie(freq_table)
    })
        
    output$value1 <- renderValueBox({
        valueBox(formatC(nrow(fit_model()), format = "d", big.mark = ","), 
            paste("Number of assertions"), icon = icon("stats", lib = "glyphicon"), 
            color = "purple")
    })
    output$value2 <- renderValueBox({
        valueBox("0.45-0.75", "Penetrance Range", icon = icon("globe", 
            lib = "glyphicon"), color = "green")
    })
    output$value3 <- renderValueBox({
        valueBox(formatC(fit_model()[1, 1], format = "d", big.mark = ","), 
            paste("Top variant rsID"), icon = icon("pushpin", lib = "glyphicon"), 
            color = "yellow")
    })
    
    output$mytable1 <- DT::renderDataTable({
        fit_model()
      #  table$AlleleID <- paste0("<a href='http://www.ncbi.nlm.nih.gov/clinvar/?term=", table$AlleleID, "[alleleid]'>", table$AlleleID, "</a>")
      #  datatable(table, escape = FALSE, options = list(pageLength = 10))
    }, server = FALSE)
    
        
    fit_model <- reactive({
        
        # Download data from database
        con <- dbConnect(RSQLite::SQLite(), "HCM.db")
        clinvar <- dbGetQuery(con, "Select * from clinvar_comprehensive")
        gnomad_cases <- dbGetQuery(con, "Select * from gnomad_comprehensive")
        colnames(clinvar) <- clinvar[1, ]
        colnames(gnomad_cases) <- gnomad_cases[1, ]
        clinvar$dbSNP.ID <- as.character(clinvar$dbSNP.ID)
        
        # Merge gnomAD cases and clinvar
        merged <- merge(clinvar, gnomad_cases, by.x = "dbSNP.ID", by.y = "rsID")
        merged_subset <- merged[!is.na(merged$dbSNP.ID), ]
        
        # Identify allele frequency across ancestry groups
        merged_subset$Allele.Frequency.African <- as.numeric(merged_subset$Allele.Count.African)/as.numeric(merged_subset$Allele.Number.African)
        merged_subset$Allele.Frequency.Latino <- as.numeric(merged_subset$Allele.Count.Latino)/as.numeric(merged_subset$Allele.Number.Latino)
        merged_subset$Allele.Frequency.Ashkenazi.Jewish <- as.numeric(merged_subset$Allele.Count.Ashkenazi.Jewish)/as.numeric(merged_subset$Allele.Number.Ashkenazi.Jewish)
        merged_subset$Allele.Frequency.East.Asian <- as.numeric(merged_subset$Allele.Count.East.Asian)/as.numeric(merged_subset$Allele.Number.East.Asian)
        merged_subset$Allele.Frequency.European.Finnish <- as.numeric(merged_subset$Allele.Count.European..Finnish.)/as.numeric(merged_subset$Allele.Number.European..Finnish.)
        merged_subset$Allele.Frequency.European.Non.Finnish <- as.numeric(merged_subset$Allele.Count.European..non.Finnish.)/as.numeric(merged_subset$Allele.Number.European..non.Finnish.)
        merged_subset$Allele.Frequency.South.Asian <- as.numeric(merged_subset$Allele.Count.South.Asian)/as.numeric(merged_subset$Allele.Number.South.Asian)
        
        # Filter by LP/P/VUS
        merged_subset$Clinical.significance[merged_subset$Clinical.significance == 
            "Pathogenic/Likely pathogenic"] <- "Likely pathogenic"
        cs <- c("Likely pathogenic", "Pathogenic", "Uncertain significance")
        file <- merged_subset[merged_subset$Clinical.significance %in% 
            cs, ]
        
        myvars1 <- c("dbSNP.ID", "Clinical.significance", "Gene.s.", "AlleleID.s.", 
            "Chromosome", "Position", "Allele.Frequency", "Allele.Frequency.African", 
            "Allele.Frequency.Latino", "Allele.Frequency.Ashkenazi.Jewish", 
            "Allele.Frequency.East.Asian", "Allele.Frequency.European.Finnish", 
            "Allele.Frequency.European.Non.Finnish", "Allele.Frequency.South.Asian")
        file_complete <- file[myvars1]
        subset <- subset(file_complete, Gene.s. %in% input$checkGroup)
        subset <- subset(subset, Clinical.significance %in% input$checkbox)
        table <- arrange(subset, desc(Allele.Frequency))
        table <- table %>% mutate_if(is.numeric, round, digits = 5)
        table$Ancestry_Group_Variant_Most_Freq <- colnames(table[, 8:14])[apply(table[, 
            8:14], 1, which.max)]
        table$Ancestry_Group_Variant_Most_Freq <- substring(table$Ancestry_Group_Variant_Most_Freq, 
            18)
       	colnames(table) <- c("rsID", "Class", "Gene", "AlleleID", "Chr", "Position", "AF", "Afr", 
            "Lat", "AJ", "EA", "EF", "ENF", "SA", "MostFrequent")
        table
    })   
    
}
