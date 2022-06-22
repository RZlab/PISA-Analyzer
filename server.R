library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(ggrepel)
library(broom)
library(shinycssloaders)

options(shiny.usecairo = T)
options(dplyr.summarise.inform = F)
options(ggrepel.max.overlaps = Inf)

ref.condition <- NULL
ref.sample <- NULL
ref.drug <- NULL
ref.experiment <- NULL
ref.parameter <- NULL

theme_beusch <- function(){
  theme_bw() %+replace% 
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          legend.position = 'bottom',
          complete = TRUE)
}

read.data <- function(df = NULL, reference.sample = 'DMSO'){
  
  # normalize by expression if applicable
  if(length(grep('Exp', df$Parameter)) != 0){
    s <- df %>%
      filter(Parameter == 'Exp') %>%
      group_by(ID, Experiment, Drug) %>%
      summarise(value.Exp = mean(value))
    
    df.pisa <- df %>%
      filter(Parameter == 'PISA') %>%
      full_join(s, by = c('ID', 'Experiment', 'Drug')) %>%
      mutate(value = value - value.Exp) %>%
      dplyr::select(-value.Exp)
    
    df.exp <- df %>%
      filter(Parameter == 'Exp')
    
    df <- rbind(df.pisa, df.exp)  
  }
  
  # normalize to corresponding DMSO control
  s <- df %>%
    filter(Drug == reference.sample) %>%
    group_by(ID, Experiment, Sample, Parameter) %>%
    summarise(value.DMSO = mean(value))
  
  df <- df %>%
    full_join(s, by = c('ID', 'Experiment', 'Sample', 'Parameter')) %>%
    mutate('normalized.value' = value-value.DMSO)
  
  return(df)
}

statistical.analysis.1DPISA <- function(data = NULL,
                                        experiment.of.interest = NULL,
                                        drug.of.interest = NULL,
                                        parameter.of.interest = NULL,
                                        pvalue.cutoff = 0.05,
                                        FC.cutoff = 0.3,
                                        reference.sample = 'DMSO'){
  
  experiment <- unique(data$Experiment)
  if(!is.null(experiment.of.interest)){experiment <- experiment.of.interest}
  if(length(experiment) == 0){message('Please specify a valid drug or leave this input empty.')}
  
  drugs <- unique(data$Drug)
  drugs <- drugs[drugs != reference.sample]
  if(!is.null(drug.of.interest)){drugs <- drug.of.interest}
  if(length(drugs) == 0){message('Please specify a valid drug or leave this input empty.')}
  
  sample <- unique(data$Sample)
  if(!is.null(parameter.of.interest)){sample <- parameter.of.interest}
  if(length(sample) == 0){message('Please specify a valid Paramter of interest or leave this input empty.')}
  
  res <- list()
  n = 1
  for(i in 1:length(experiment)){
    t <- data %>% filter(Experiment == experiment[i])
    
    for(j in 1:length(sample)){
      tt <- t %>% filter(Sample == sample[j])
      
      for(k in 1:length(drugs)){
        ttt <- tt %>%
          filter(Drug == drugs[k] | Drug == reference.sample) %>%
          mutate('comparison' = ifelse(Drug == reference.sample, 0, 1))
        
        test <- ttt %>%
          group_by(ID) %>%
          do({
            s <- .
            
            p <- try(t.test(normalized.value ~ comparison, data = s, var.equal = T), silent = T)
            
            if(class(p) == 'htest'){
              r <- data.frame('Gene' = s$Gene[1],
                              'Peptides' = s$Peptides[1],
                              'Experiment' = experiment[i],
                              'Drug' = drugs[k],
                              'Sample' = sample[j])
              
              r <- cbind(r, tidy(p))
            }
            else{
              r <- data.frame()
            }
            
            r
          })
        
        # curate and write statistical table
        test <- test %>%
          mutate(FC = (estimate2 - estimate1)) %>%
          mutate('legend' = ifelse(p.value < pvalue.cutoff & abs(FC) > FC.cutoff, 'significant', 'non-significant')) %>%
          dplyr::select(ID, Gene, Peptides, Experiment, Drug, Sample, p.value, FC, legend)
        
        res[[n]] <- test
        n = n+1
      }
    }
  }
  res <- plyr::rbind.fill(res) %>% mutate('Parameter' = 'PISA')
  
  return(res)
}

statistical.analysis.1DPISA.cells <- function(data = NULL,
                                              experiment.of.interest = NULL,
                                              drug.of.interest = NULL,
                                              parameter.of.interest = NULL,
                                              pvalue.cutoff = 0.05,
                                              FC.cutoff = 0.3,
                                              reference.sample = 'DMSO'){
  
  experiment <- unique(data$Experiment)
  if(!is.null(experiment.of.interest)){experiment <- experiment.of.interest}
  if(length(experiment) == 0){message('Please specify a valid drug or leave this input empty.')}
  
  drugs <- unique(data$Drug)
  drugs <- drugs[drugs != reference.sample]
  if(!is.null(drug.of.interest)){drugs <- drug.of.interest}
  if(length(drugs) == 0){message('Please specify a valid drug or leave this input empty.')}
  
  sample <- unique(data$Sample)
  if(!is.null(parameter.of.interest)){sample <- parameter.of.interest}
  if(length(sample) == 0){message('Please specify a valid Paramter of interest or leave this input empty.')}
  
  # pisa
  res <- list()
  n = 1
  for(i in 1:length(experiment)){
    t <- data %>%
      filter(Parameter == 'PISA') %>%
      filter(Experiment == experiment[i])
    
    for(j in 1:length(sample)){
      tt <- t %>% filter(Sample == sample[j])
      
      for(k in 1:length(drugs)){
        ttt <- tt %>%
          filter(Drug == drugs[k] | Drug == reference.sample) %>%
          mutate('comparison' = ifelse(Drug == reference.sample, 0, 1))
        
        test <- ttt %>%
          group_by(ID) %>%
          do({
            s <- .
            
            p <- try(t.test(normalized.value ~ comparison, data = s, var.equal = T), silent = T)
            
            if(class(p) == 'htest'){
              r <- data.frame('Gene' = s$Gene[1],
                              'Peptides' = s$Peptides[1],
                              'Experiment' = experiment[i],
                              'Drug' = drugs[k],
                              'Sample' = sample[j])
              
              r <- cbind(r, tidy(p))
            }
            else{
              r <- data.frame()
            }
            
            r
          })
        
        # curate and write statistical table
        test <- test %>%
          mutate(FC = (estimate2 - estimate1)) %>%
          mutate('legend' = ifelse(p.value < pvalue.cutoff & abs(FC) > FC.cutoff, 'significant', 'non-significant')) %>%
          dplyr::select(ID, Gene, Peptides, Experiment, Drug, Sample, p.value, FC, legend)
        
        res[[n]] <- test
        n = n+1
      }
    }
  }
  res.pisa <- plyr::rbind.fill(res) %>% mutate('Parameter' = 'PISA')
  
  # expression
  res <- list()
  n = 1
  for(i in 1:length(experiment)){
    t <- data %>%
      filter(Parameter == 'Exp') %>%
      filter(Experiment == experiment[i])
    
    for(j in 1:length(sample)){
      tt <- t %>% filter(Sample == sample[j])
      
      for(k in 1:length(drugs)){
        ttt <- tt %>%
          filter(Drug == drugs[k] | Drug == reference.sample) %>%
          mutate('comparison' = ifelse(Drug == reference.sample, 0, 1))
        
        test <- ttt %>%
          group_by(ID) %>%
          do({
            s <- .
            
            p <- try(t.test(normalized.value ~ comparison, data = s, var.equal = T), silent = T)
            
            if(class(p) == 'htest'){
              r <- data.frame('Gene' = s$Gene[1],
                              'Peptides' = s$Peptides[1],
                              'Experiment' = experiment[i],
                              'Drug' = drugs[k],
                              'Sample' = sample[j])
              
              r <- cbind(r, tidy(p))
            }
            else{
              r <- data.frame()
            }
            
            r
          })
        
        # curate and write statistical table
        test <- test %>%
          mutate(FC = (estimate2 - estimate1)) %>%
          mutate('legend' = ifelse(p.value < pvalue.cutoff & abs(FC) > FC.cutoff, 'significant', 'non-significant')) %>%
          dplyr::select(ID, Gene, Peptides, Experiment, Drug, Sample, p.value, FC, legend)
        
        res[[n]] <- test
        n = n+1
      }
    }
  }
  res.expression <- plyr::rbind.fill(res) %>% mutate('Parameter' = 'Exp')
  
  res <- rbind(res.pisa, res.expression)
  
  return(res)
}

statistical.analysis.2DPISA <- function(data = NULL,
                                        experiment.of.interest = NULL,
                                        drug.of.interest = NULL,
                                        parameter.of.interest = NULL,
                                        pvalue.cutoff = 0.05,
                                        FC.cutoff = 0.3,
                                        reference.sample = 'DMSO'){
  
  experiment <- unique(data$Experiment)
  if(!is.null(experiment.of.interest)){experiment <- experiment.of.interest}
  if(length(experiment) == 0){message('Please specify a valid drug or leave this input empty.')}
  
  drugs <- unique(data$Drug)
  drugs <- drugs[drugs != reference.sample]
  if(!is.null(drug.of.interest)){drugs <- drug.of.interest}
  if(length(drugs) == 0){message('Please specify a valid drug or leave this input empty.')}
  
  sample <- unique(data$Sample)
  if(!is.null(parameter.of.interest)){sample <- parameter.of.interest}
  if(length(sample) == 0){message('Please specify a valid Paramter of interest or leave this input empty.')}
  
  res <- list()
  n = 1
  for(i in 1:length(experiment)){
    t <- data %>% filter(Experiment == experiment[i])
    
    for(j in 1:length(sample)){
      tt <- t %>% filter(Sample == sample[j])
      
      for(k in 1:length(drugs)){
        ttt <- tt %>%
          filter(Drug == drugs[k] | Drug == reference.sample) %>%
          mutate('comparison' = ifelse(Drug == reference.sample, 0, 1))
        
        test <- ttt %>%
          group_by(ID) %>%
          do({
            s <- .
            
            p <- try(t.test(normalized.value ~ comparison, data = s, var.equal = T), silent = T)
            
            if(class(p) == 'htest'){
              r <- data.frame('Gene' = s$Gene[1],
                              'Peptides' = s$Peptides[1],
                              'Experiment' = experiment[i],
                              'Drug' = drugs[k],
                              'Sample' = sample[j])
              
              r <- cbind(r, tidy(p))
            }
            else{
              r <- data.frame()
            }
            
            r
          })
        
        # curate and write statistical table
        test <- test %>%
          mutate(FC = (estimate2 - estimate1)) %>%
          mutate('legend' = ifelse(p.value < pvalue.cutoff & abs(FC) > FC.cutoff, 'significant', 'non-significant')) %>%
          dplyr::select(ID, Gene, Peptides, Experiment, Drug, Sample, p.value, FC, legend)
        
        res[[n]] <- test
        n = n+1
      }
    }
  }
  res <- plyr::rbind.fill(res) %>% mutate('Parameter' = 'PISA')
  
  return(res)
}

analyze.cResT.data <- function(data = NULL,
                               dat.statistics = dat.statistics,
                               experiment.of.interest = NULL,
                               drug.of.interest = NULL,
                               parameter.of.interest = c('conc-PISA')){
  dat <- data %>% unite('protein.id', ID, Experiment, remove = FALSE)
  
  poi <- dat.statistics %>%
    filter(legend == 'significant') %>%
    unite('protein.id', ID, Experiment, remove = FALSE) %>%
    filter(grepl('Max', Drug))
  
  if(!is.null(experiment.of.interest)){poi <- poi %>% filter(Experiment %in% experiment.of.interest)}
  if(!is.null(drug.of.interest)){poi <- poi %>% filter(Drug %in% drug.of.interest)}
  if(!is.null(parameter.of.interest)){poi <- poi %>% filter(Sample %in% parameter.of.interest)}
  
  poi <- poi$protein.id
  
  # scale Pool values relative to max
  dat0 <- dat %>%
    group_by(protein.id, Sample, Drug) %>%
    summarise('mean.max.alteration' = mean(normalized.value),
              'sd.max.alteration' = sd(normalized.value)) %>%
    filter(grepl('Max', Drug)) %>%
    mutate('scaling.factor' = 1/mean.max.alteration) %>%
    dplyr::select(-Drug) 
  
  dat <- dat %>%
    full_join(dat0, by = c('protein.id', 'Sample')) %>%
    mutate('relative.normalized.value' = normalized.value*scaling.factor)
  
  # filter for protein with too much variation
  d <- dat %>% filter(sd.max.alteration < 3*mean.max.alteration)
  s <- d %>%
    group_by(protein.id) %>%
    summarise('n' = n())
  
  s <- s %>% filter(n == max(s$n))
  
  d <- d %>% filter(protein.id %in% s$protein.id)
  
  # remove too high pool values
  s <- d %>%
    filter(grepl('Pool', Drug)) %>%
    group_by(protein.id, Sample) %>%
    summarise(mean.relative.normalized.value = mean(relative.normalized.value)) %>%
    filter(mean.relative.normalized.value <= 1 & mean.relative.normalized.value >= 0) %>%
    group_by(protein.id) %>%
    summarise('n' = n())
  
  #s <- s %>% filter(n == 2)
  
  d <- d %>% filter(protein.id %in% s$protein.id)
  
  # filter based on statistical values
  d <- d %>% filter(protein.id %in% poi)
  
  # calculate 2D score
  dat.Concentration <- d %>%
    filter(grepl('Pool', Drug)) %>%
    filter(Sample == 'conc-PISA') %>%
    group_by(protein.id, ID, Gene, Peptides, Experiment, Drug) %>%
    summarise('Concentration.ratio' = mean(relative.normalized.value))
  
  dat.Residence <- d %>%
    filter(grepl('Pool', Drug)) %>%
    filter(Sample == 'cResT-PISA') %>%
    group_by(protein.id, ID, Gene, Peptides, Experiment, Drug) %>%
    summarise('Residence.ratio' = mean(relative.normalized.value))
  
  pb.score <- dat.Concentration %>%
    full_join(dat.Residence, by = c('protein.id', 'ID', 'Gene', 'Peptides', 'Experiment', 'Drug')) %>%
    mutate('Concentration.ratio' = ifelse(is.na(Concentration.ratio), 0, Concentration.ratio)) %>%
    mutate('Residence.ratio' = ifelse(is.na(Residence.ratio), 0, Residence.ratio))
  
  pb.score <- pb.score %>%
    filter(Concentration.ratio >= 0 & Concentration.ratio <= 1) %>%
    filter(Residence.ratio >= 0 & Residence.ratio <= 1)
  
  s <- pb.score %>%
    group_by(protein.id) %>%
    summarise('n' = n())
  ss <- s %>% filter(n == max(s$n))
  pb.score <- pb.score %>% filter(protein.id %in% ss$protein.id)
  
  # calculate 3D score
  dat <- dat %>% filter(protein.id %in% d$protein.id)
  
  if(sum(pb.score$Concentration.ratio) > 0){
    dat.alteration <- d %>%
      filter(grepl('Max', Drug)) %>%
      filter(Sample == 'conc-PISA') %>%
      group_by(protein.id, ID, Gene, Peptides, Experiment) %>%
      summarise('mean.max.alteration' = mean(mean.max.alteration)) %>%
      mutate('scaled.mean.max.alteration' = mean.max.alteration/max(abs(dat$mean.max.alteration)))
  }else{
    dat.alteration <- d %>%
      filter(grepl('Max', Drug)) %>%
      filter(Sample == 'cResT-PISA') %>%
      group_by(protein.id, ID, Gene, Peptides, Experiment) %>%
      summarise('mean.max.alteration' = mean(mean.max.alteration)) %>%
      mutate('scaled.mean.max.alteration' = mean.max.alteration/max(abs(dat$mean.max.alteration)))
  }
  
  pb.score <- pb.score %>%
    inner_join(dat.alteration, by = c('protein.id', 'ID', 'Gene', 'Peptides', 'Experiment')) %>%
    mutate('pb.score' = Concentration.ratio+Residence.ratio+scaled.mean.max.alteration) %>%
    dplyr::select(protein.id, ID, Gene, Peptides, Experiment, Concentration.ratio, Residence.ratio, scaled.mean.max.alteration, pb.score) %>%
    gather(experiment, pb.value, -c(protein.id:Experiment, pb.score)) %>%
    arrange(desc(pb.score))
  
  return(pb.score)
}

server <- function(input, output, session){
  
  #extending max file size
  options(shiny.maxRequestSize = 50*1024^2)
  
  #upload tsv file
  data.set <- reactive({
    inFile <- input$upload.file$datapath
    
    return(inFile)
  })
  
  data <- reactive({
    ds <- data.set()
    if(class(ds) == 'character'){
      df <- suppressMessages(read_tsv(ds))
    }else{
      df <- NULL
    }
    
    return(df)
  })
  
  dat.user <- reactive({
    df <- data()
    if(!is.null(df)){
      df <- read.data(df = df, reference.sample = input$reference.condition)
    }else{
      df <- NULL
    }
    
    return(df)
  })
  
  reference.condition <- reactive({
    o <- NULL
    if(!is.null(dat.user())){
      d <- dat.user()
      o <- append('', unique(d$Drug))
    }
    else{
      o <- NULL
    }
    
    return(o)
  })
  
  #reference.condition
  observe({
    ref.condition <- reference.condition()
    x <- input$reference.condition
    
    updateSelectInput(session, "reference.condition",
                      choices = ref.condition,
                      selected = x
    )
  })
  
  reference.sample <- reactive({
    o <- NULL
    if(!is.null(dat.user())){
      d <- dat.user()
      o <- unique(d$Sample)
    }
    else{
      o <- NULL
    }
    
    return(o)
  })
  
  observe({
    ref.sample <- reference.sample()
    x <- input$reference.sample
    
    updateSelectInput(session, "reference.sample",
                      choices = ref.sample,
                      selected = x
    )
  })
  
  reference.drug <- reactive({
    o <- NULL
    if(!is.null(dat.user())){
      d <- dat.user()
      o <- unique(d$Drug)
      o <- o[-which(o == input$reference.condition)] 
    }
    else{
      o <- NULL
    }
    
    return(o)
  })
  
  observe({
    ref.drug <- reference.drug()
    x <- input$reference.drug
    
    updateSelectInput(session, "reference.drug",
                      choices = ref.drug,
                      selected = x
    )
  })
  
  reference.experiment <- reactive({
    o <- NULL
    if(!is.null(dat.user())){
      d <- dat.user()
      o <- unique(d$Experiment)
    }
    else{
      o <- NULL
    }
    
    return(o)
  })
  
  observe({
    ref.experiment <- reference.experiment()
    x <- input$reference.experiment
    
    updateSelectInput(session, "reference.experiment",
                      choices = ref.experiment,
                      selected = x
    )
  })
  
  reference.parameter <- reactive({
    o <- NULL
    if(!is.null(dat.user())){
      d <- dat.user()
      o <- unique(d$Parameter)
    }
    else{
      o <- NULL
    }
    
    return(o)
  })
  
  observe({
    ref.parameter <- reference.parameter()
    x <- input$reference.parameter
    
    updateSelectInput(session, "reference.parameter",
                      choices = ref.parameter,
                      selected = x
    )
  })
  
  #process tsv file
  dat.statistic <- eventReactive(c(input$start.analysis),{
    dat <- dat.user()
    
    if(input$start.analysis[1] > 0){
      par <- unique(dat$Parameter)
      if(length(grep('Exp', par)) > 0){
        dat <- statistical.analysis.1DPISA.cells(data = dat,
                                                 experiment.of.interest = NULL,
                                                 drug.of.interest = NULL,
                                                 parameter.of.interest = NULL,
                                                 pvalue.cutoff = 0.05,
                                                 FC.cutoff = 0.3,
                                                 reference.sample = input$reference.condition)
      }
      drug <- unique(dat$Drug)
      if(length(grep('Pool', drug)) > 0){
        dat <- statistical.analysis.2DPISA(data = dat,
                                           experiment.of.interest = NULL,
                                           drug.of.interest = NULL,
                                           parameter.of.interest = NULL,
                                           pvalue.cutoff = 0.05,
                                           FC.cutoff = 0.3,
                                           reference.sample = input$reference.condition)
      }
      if(length(grep('Pool', drug)) == 0 & length(grep('Exp', par)) == 0){
        dat <- statistical.analysis.1DPISA(data = dat,
                                           experiment.of.interest = NULL,
                                           drug.of.interest = NULL,
                                           parameter.of.interest = NULL,
                                           pvalue.cutoff = 0.05,
                                           FC.cutoff = 0.3,
                                           reference.sample = input$reference.condition)
      }
    }else{
      dat <- NULL
    }
    
    return(dat)  
  })
  
  observe(dat.statistic(), suspended = F)
  
  dat.statistic.filtered <- reactive({
    dat <- dat.statistic()
    
    if(!is.null(dat)){
      dat <- dat %>% mutate('legend' = ifelse(p.value < input$p.value.cutoff & abs(FC) > input$delta.Sm.cutoff,
                                              'significant',
                                              'non-significant'))
    }
    else{
      dat <- NULL
    }
    return(dat)
  })
  
  observe(dat.statistic.filtered(), suspended = TRUE)
  
  #visualize results
  output$volcano.plot <- renderPlotly({
    req(dat.statistic.filtered)
    dat <- dat.statistic.filtered()

    if(!is.null(input$reference.sample)){
      input.reference.sample <- input$reference.sample
    }else{
      input.reference.sample <- unique(dat.user()$Sample)
    }
    
    if(!is.null(input$reference.experiment)){
      input.reference.experiment <- input$reference.experiment
    }else{
      input.reference.experiment <- unique(dat.user()$Experiment)
    }
    
    if(!is.null(input$reference.parameter)){
      input.reference.parameter <- input$reference.parameter
    }else{
      input.reference.parameter <- unique(dat.user()$Parameter)
    }
    
    if(!is.null(dat)){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        dat <- dat %>% filter(Sample == input.reference.sample)
        dat <- dat %>% filter(Drug == input$reference.drug)
        dat <- dat %>% filter(Experiment == input.reference.experiment)
        dat <- dat %>% filter(Parameter == input.reference.parameter)
        
        gg <- ggplotly(ggplot(dat, aes(x = FC, y = -log10(p.value), colour = legend,
                                       label = Gene, label2 = p.value)) +
                         geom_point() +
                         geom_hline(yintercept = -log10(input$p.value.cutoff), linetype = 'dashed') +
                         theme_beusch() +
                         scale_colour_manual(values = c('#3d3d3d', '#EB491E')) +
                         scale_y_continuous(limits = c(0, max(-log10(dat$p.value)*1.1)), expand = c(0, 0)) +
                         xlim(-max(abs(dat$FC)*1.1), max(abs(dat$FC)*1.1)) +
                         xlab('log2-scaled fold change') +
                         theme(legend.position = 'bottom'),  tooltip = c('Gene' ,'FC', 'p.value')) %>%
          config(displayModeBar = F) %>%
          layout(
            legend = list(
              orientation = "h",
              x = 0.2, y = -0.2
            )
          )
        
      }else{gg <- NULL}
      
    }else{gg <- NULL}
    
    return(gg)
  })
  
  volcano.plot <- reactive({
    req(dat.statistic.filtered)
    dat <- dat.statistic.filtered()
    
    if(!is.null(input$reference.sample)){
      input.reference.sample <- input$reference.sample
    }else{
      input.reference.sample <- unique(dat.user()$Sample)
    }
    
    if(!is.null(input$reference.experiment)){
      input.reference.experiment <- input$reference.experiment
    }else{
      input.reference.experiment <- unique(dat.user()$Experiment)
    }
    
    if(!is.null(input$reference.parameter)){
      input.reference.parameter <- input$reference.parameter
    }else{
      input.reference.parameter <- unique(dat.user()$Parameter)
    }
    
    if(!is.null(dat)){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        dat <- dat %>% filter(Sample == input.reference.sample)
        dat <- dat %>% filter(Drug == input$reference.drug)
        dat <- dat %>% filter(Experiment == input.reference.experiment)
        dat <- dat %>% filter(Parameter == input.reference.parameter)
        
        gg <- ggplot(dat, aes(x = FC, y = -log10(p.value), colour = legend)) +
          geom_point() +
          geom_text_repel(data = dat %>% filter(legend == 'significant'),
                          aes(label = `Gene`), size = 3, box.padding = unit(0.35, 'lines'),
                          point.padding = unit(0.3, 'lines')) +
          geom_hline(yintercept = -log10(input$p.value.cutoff), linetype = 'dashed') +
          theme_beusch() +
          scale_colour_manual(values = c('#3d3d3d', '#EB491E')) +
          scale_y_continuous(limits = c(0, max(-log10(dat$p.value)*1.1)), expand = c(0, 0)) +
          xlim(-max(abs(dat$FC)*1.1), max(abs(dat$FC)*1.1)) +
          xlab('log2-scaled fold change') +
          theme(legend.position = 'bottom')
        
      }else{gg <- NULL}
      
    }else{gg <- NULL}
    
    return(gg)
  })
  
  comparison.data <- reactive({
    req(dat.statistic.filtered)
    dat <- dat.statistic.filtered()
    
    if(!is.null(input$reference.sample)){
      input.reference.sample <- input$reference.sample
    }else{
      input.reference.sample <- unique(dat.user()$Sample)
    }
    
    if(!is.null(input$reference.experiment)){
      input.reference.experiment <- input$reference.experiment
    }else{
      input.reference.experiment <- unique(dat.user()$Experiment)
    }
    
    if(!is.null(input$reference.parameter)){
      input.reference.parameter <- input$reference.parameter
    }else{
      input.reference.parameter <- unique(dat.user()$Parameter)
    }
    
    if(length(grep('Exp', unique(dat$Parameter))) > 0){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        dat <- dat %>% filter(Sample == input.reference.sample)
        dat <- dat %>% filter(Drug == input$reference.drug)
        dat <- dat %>% filter(Experiment == input.reference.experiment)
        
        dat.pisa <- dat %>% filter(grepl('PISA', Parameter)) %>% filter(legend == 'significant')
        dat.Exp <- dat %>% filter(grepl('Exp', Parameter)) %>% filter(legend == 'significant')
        
        test <- dat %>%
          dplyr::select(ID, Gene, Peptides, Parameter, FC) %>%
          tidyr::spread(Parameter, FC) %>%
          mutate('legend' = 'none', 
                 'legend' = ifelse(ID %in% dat.pisa$ID, 'PISA', legend),
                 'legend' = ifelse(ID %in% dat.Exp$ID, 'Exp', legend),
                 'legend' = ifelse(ID %in% dat.pisa$ID & ID %in% dat.Exp$ID, 'both', legend))
      }
    }
    if(length(grep('Pool', unique(dat$Drug))) > 0){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        dat <- dat %>% filter(Sample == input.reference.sample)
        dat <- dat %>% filter(Experiment == input.reference.experiment)
        
        dat.max <- dat %>% filter(!grepl('Pool', Drug)) %>% filter(legend == 'significant')
        dat.pool <- dat %>% filter(grepl('Pool', Drug)) %>% filter(legend == 'significant')
        
        test <- dat %>%
          dplyr::select(ID, Gene, Peptides, Drug, FC) %>%
          tidyr::spread(Drug, FC) %>%
          mutate('legend' = 'none', 
                 'legend' = ifelse(ID %in% dat.max$ID, 'max', legend),
                 'legend' = ifelse(ID %in% dat.pool$ID, 'pool', legend),
                 'legend' = ifelse(ID %in% dat.max$ID & ID %in% dat.pool$ID, 'both', legend))
      }
    }
    #else{test <- NULL}
    if(length(grep('Exp', unique(dat$Parameter))) == 0 & length(grep('Pool', unique(dat$Drug))) == 0){test <- NULL}
    
    return(test)
  })
  
  #print data for scatter plot
  output$comparison.plot <- renderPlotly({
    req(comparison.data)
    dat <- dat.statistic.filtered()
    
    if(!is.null(input$reference.sample)){
      input.reference.sample <- input$reference.sample
    }else{
      input.reference.sample <- unique(dat.user()$Sample)
    }
    
    if(!is.null(input$reference.experiment)){
      input.reference.experiment <- input$reference.experiment
    }else{
      input.reference.experiment <- unique(dat.user()$Experiment)
    }
    
    if(!is.null(input$reference.parameter)){
      input.reference.parameter <- input$reference.parameter
    }else{
      input.reference.parameter <- unique(dat.user()$Parameter)
    }
    
    if(length(grep('Exp', unique(dat$Parameter))) > 0){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        test <- comparison.data()
        lim <- max(c(abs(test[,4]),abs(test[,5])))
        
        gg <- ggplotly(ggplot(test, aes_string(x = colnames(test)[4], y = colnames(test)[5])) +
                         geom_point(aes(colour = legend, label = Gene)) +
                         geom_hline(yintercept = 0, linetype = 'dashed') +
                         geom_vline(xintercept = 0, linetype = 'dashed') +
                         scale_colour_manual(values = c('both' = '#EB491E',
                                                        'PISA' = '#03b1fc',
                                                        'Exp' = '#3d3d3d',
                                                        'pool' = '#8003fc')) +
                         xlim(-lim, lim) +
                         ylim(-lim, lim) +
                         ylab('log2-scaled fold change PISA sample') +
                         xlab('log2-scaled fold change expression sample') +
                         theme_beusch() +
                         theme(legend.position = 'bottom'),
                       tooltip = c('Gene' , 'x', 'y')) %>%
          config(displayModeBar = F) %>%
          layout(
            legend = list(
              orientation = "h",
              x = 0.2, y = -0.2
            )
          )
        
      }else{
        gg <- NULL
      }
      
    }else{gg <- NULL}
    if(length(grep('Pool', unique(dat$Drug))) > 0){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        test <- comparison.data()
        lim <- max(c(abs(test[,4]),abs(test[,5])))

        gg <- ggplotly(ggplot(test, aes_string(x = colnames(test)[4], y = colnames(test)[5])) +
                         geom_point(aes(colour = legend, label = Gene)) +
                         geom_hline(yintercept = 0, linetype = 'dashed') +
                         geom_vline(xintercept = 0, linetype = 'dashed') +
                         scale_colour_manual(values = c('both' = '#EB491E',
                                                        'max' = '#03b1fc',
                                                        'none' = '#3d3d3d',
                                                        'pool' = '#8003fc')) +
                         xlim(-lim, lim) +
                         ylim(-lim, lim) +
                         ylab('log2-scaled fold change pooled sample') +
                         xlab('log2-scaled fold change max concentration sample') +
                         theme_beusch() +
                         theme(legend.position = 'bottom'),
                       tooltip = c('Gene' , 'x', 'y')) %>%
          config(displayModeBar = F) %>%
          layout(
            legend = list(
              orientation = "h",
              x = 0.2, y = -0.2
            )
          )
      }else{gg <- NULL}
    }
    
    return(gg)
  })
  
  comparison.plot <- reactive({
    req(comparison.data)
    dat <- dat.statistic.filtered()
    
    if(!is.null(input$reference.sample)){
      input.reference.sample <- input$reference.sample
    }else{
      input.reference.sample <- unique(dat.user()$Sample)
    }
    
    if(!is.null(input$reference.experiment)){
      input.reference.experiment <- input$reference.experiment
    }else{
      input.reference.experiment <- unique(dat.user()$Experiment)
    }
    
    if(!is.null(input$reference.parameter)){
      input.reference.parameter <- input$reference.parameter
    }else{
      input.reference.parameter <- unique(dat.user()$Parameter)
    }
    
    if(length(grep('Exp', unique(dat$Parameter))) > 0){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        test <- comparison.data()
        lim <- max(c(abs(test[,4]),abs(test[,5])))
        
        gg <- ggplot(test, aes_string(x = colnames(test)[4], y = colnames(test)[5])) +
          geom_point(aes(colour = legend, label = Gene)) +
          geom_text_repel(data = subset(test, legend != 'none'),
                          aes(label = `Gene`), size = 3, box.padding = unit(0.35, 'lines'),
                          point.padding = unit(0.3, 'lines')) +
          geom_hline(yintercept = 0, linetype = 'dashed') +
          geom_vline(xintercept = 0, linetype = 'dashed') +
          scale_colour_manual(values = c('both' = '#EB491E',
                                         'PISA' = '#03b1fc',
                                         'Exp' = '#3d3d3d',
                                         'pool' = '#8003fc')) +
          xlim(-lim, lim) +
          ylim(-lim, lim) +
          ylab('log2-scaled fold change PISA sample') +
          xlab('log2-scaled fold change expression sample') +
          theme_beusch() +
          theme(legend.position = 'bottom')
        
      }else{
        gg <- NULL
      }
      
    }else{gg <- NULL}
    if(length(grep('Pool', unique(dat$Drug))) > 0){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        test <- comparison.data()
        lim <- max(c(abs(test[,4]),abs(test[,5])))
        
        gg <- ggplot(test, aes_string(x = colnames(test)[4], y = colnames(test)[5])) +
          geom_point(aes(colour = legend, label = Gene)) +
          geom_text_repel(data = subset(test, legend != 'none'),
                          aes(label = `Gene`), size = 3, box.padding = unit(0.35, 'lines'),
                          point.padding = unit(0.3, 'lines')) +
          geom_hline(yintercept = 0, linetype = 'dashed') +
          geom_vline(xintercept = 0, linetype = 'dashed') +
          scale_colour_manual(values = c('both' = '#EB491E',
                                         'max' = '#03b1fc',
                                         'none' = '#3d3d3d',
                                         'pool' = '#8003fc')) +
          xlim(-lim, lim) +
          ylim(-lim, lim) +
          ylab('log2-scaled fold change pooled sample') +
          xlab('log2-scaled fold change max concentration sample') +
          theme_beusch() +
          theme(legend.position = 'bottom')
      }else{gg <- NULL}
    }
    
    return(gg)
  })
  
  # score data  
  score.data <- reactive({
    req(dat.statistic.filtered)
    dat <- dat.statistic.filtered()
    
    if(!is.null(input$reference.sample)){
      input.reference.sample <- input$reference.sample
    }else{
      input.reference.sample <- unique(dat.user()$Sample)
    }
    
    if(!is.null(input$reference.experiment)){
      input.reference.experiment <- input$reference.experiment
    }else{
      input.reference.experiment <- unique(dat.user()$Experiment)
    }
    
    if(!is.null(input$reference.parameter)){
      input.reference.parameter <- input$reference.parameter
    }else{
      input.reference.parameter <- unique(dat.user()$Parameter)
    }
    
    if(length(grep('Pool', unique(dat$Drug))) > 0){
      if(input.reference.sample != '' & input$reference.drug != '' &
         input.reference.experiment != '' & input.reference.parameter != ''){
        
        ref <- ifelse(length(grep('conc-PISA', unique(dat.user()$Sample))) > 0, 'conc-PISA', 'cResT-PISA')
        
        pb.score <- analyze.cResT.data(data = dat.user(),
                                       dat.statistics = dat,
                                       experiment.of.interest = NULL,
                                       drug.of.interest = NULL,
                                       parameter.of.interest = ref)
        
        pb.score <- pb.score %>%
          mutate(experiment = ifelse(experiment == 'Concentration.ratio', 'Conc-PISA', experiment)) %>%
          mutate(experiment = ifelse(experiment == 'Residence.ratio', 'ResT-PISA', experiment)) %>%
          mutate(experiment = ifelse(experiment == 'scaled.mean.max.alteration', 'Scaled deltaSm', experiment))
          
      }else{pb.score <- NULL}
    }else{pb.score <- NULL}

    return(pb.score)
  })
  
  score.plot <- reactive({
    req(score.data)
    dat <- score.data()
    
    if(!is.null(dat)){
      max.score = dat %>%
        filter(pb.value > 0) %>%
        distinct(experiment)
      
      gg <- ggplot(dat, aes(x = reorder(Gene, -pb.score), y = pb.value, fill = experiment)) +
        geom_bar(stat = 'identity') +
        facet_wrap(~Experiment) +
        scale_fill_manual(values = c('#438CCA', '#F89938', '#91288D')) +
        ylim(0, length(unique(max.score$experiment))) +
        theme_beusch()
    }else{
      gg <- NULL
    }
    
    return(gg)
  })
  
  observe(score.data(), suspended = F)
  
  output$score.plot <- renderPlot({
    req(score.plot())
    score.plot()
  })
  
  
  #add start button
  output$action <- renderUI({
    if(input$reference.condition != ''){
      actionButton('start.analysis', 'Start Analysis')
    }
  })
  
  output$user.input.1 <- renderUI({
    if(input$reference.condition != '' & length(unique(dat.user()$Sample)) > 1){
      selectInput('reference.sample',
                  'Please specify which experiment you want to display',
                  width = 400,
                  ref.sample,
                  multiple = FALSE)
    }
  })
  
  output$user.input.2 <- renderUI({
    if(input$reference.condition != ''){
      selectInput('reference.drug',
                  'Please specify which sample you want to compare in the analysis',
                  width = 400,
                  ref.drug,
                  multiple = FALSE)
    }
  })
  
  output$user.input.3 <- renderUI({
    if(input$reference.condition != '' & length(unique(dat.user()$Experiment)) > 1){
      selectInput('reference.experiment',
                  'Please specify which PISA style experimental type you want to display',
                  width = 400,
                  ref.experiment,
                  multiple = FALSE)
    }
  })
  
  output$user.input.4 <- renderUI({
    if(input$reference.condition != ''  & length(unique(dat.user()$Parameter)) > 1){
      selectInput('reference.parameter',
                  'Please specify the experimental paramteter for the volcano plot',
                  width = 400,
                  ref.parameter,
                  multiple = FALSE)
    }
  })
  
  #popup message
  observeEvent(input$start.analysis, priority = 1, {
    showModal(modalDialog(
      title = '',
      'Your data is getting processed. Please have patience! \n
        Depending on the size of your input data this could take a few seecond',
      easyClose = TRUE
    ))
  })
  
  
  ## DOWNLOADS
  #scatter plot
  output$download.volcano.source.data <- renderUI({
    req(dat.statistic.filtered())
    downloadButton('out.download.volcano.source.data', 'Download source data')
  })
  output$out.download.volcano.source.data <- downloadHandler(
    filename = function() {'volcano_plot_source.tsv'},
    content = function(file){write_tsv(dat.statistic.filtered(), file)
    })
  
  output$download.volcano.plot.svg <- renderUI({
    req(volcano.plot())
    downloadButton('out.download.volcano.plot.svg', 'Download plot as svg')
  })
  output$out.download.volcano.plot.svg <- downloadHandler(
    filename = function(){'volcano_plot.svg'},
    content = function(file){ggsave(file, volcano.plot(), device = 'svg')
    }
  )
  
  output$download.volcano.plot.pdf <- renderUI({
    req(volcano.plot())
    downloadButton('out.download.volcano.plot.pdf', 'Download plot as pdf')
  })
  output$out.download.volcano.plot.pdf <- downloadHandler(
    filename = function(){'volcano_plot.pdf'},
    content = function(file){ggsave(file, volcano.plot(), device = 'pdf')
    }
  )
  
  #
  output$download.scatter.source.data <- renderUI({
    req(comparison.data())
    downloadButton('out.download.scatter.source.data', 'Download source data')
  })
  output$out.download.scatter.source.data <- downloadHandler(
    filename = function() {'scatter_plot_source.tsv'},
    content = function(file){write_tsv(comparison.data(), file)
    })
  
  output$download.scatter.plot.svg <- renderUI({
    req(comparison.plot())
    downloadButton('out.download.scatter.plot.svg', 'Download plot as svg')
  })
  output$out.download.scatter.plot.svg <- downloadHandler(
    filename = function(){'scatter_plot.svg'},
    content = function(file){ggsave(file, comparison.plot(), device = 'svg')
    }
  )
  
  output$download.scatter.plot.pdf <- renderUI({
    req(comparison.plot())
    downloadButton('out.download.scatter.plot.pdf', 'Download plot as pdf')
  })
  output$out.download.scatter.plot.pdf <- downloadHandler(
    filename = function(){'scatter_plot.pdf'},
    content = function(file){ggsave(file, comparison.plot(), device = 'pdf')
    }
  )
  
  #
  output$download.score.source.data <- renderUI({
    req(score.data())
    downloadButton('out.download.score.source.data', 'Download source data')
  })
  output$out.download.score.source.data <- downloadHandler(
    filename = function() {'score_plot_source.tsv'},
    content = function(file){write_tsv(score.data(), file)
    })
  
  output$download.score.plot.svg <- renderUI({
    req(score.plot())
    downloadButton('out.download.score.plot.svg', 'Download plot as svg')
  })
  output$out.download.score.plot.svg <- downloadHandler(
    filename = function(){'score_plot.svg'},
    content = function(file){ggsave(file, score.plot(), device = 'svg')
    }
  )
  
  output$download.score.plot.pdf <- renderUI({
    req(score.plot())
    downloadButton('out.download.score.plot.pdf', 'Download plot as pdf')
  })
  output$out.download.score.plot.pdf <- downloadHandler(
    filename = function(){'score_plot.pdf'},
    content = function(file){ggsave(file, score.plot(), device = 'pdf')
    }
  )
  
  #link to paper
  url <- a('empty',
           href = 'https://github.com/beuschc')
  output$citation <- renderUI({
    tagList('Please cite:', url)
  })
}