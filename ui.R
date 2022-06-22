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

ui <- dashboardPage(
  title = 'PISA Analysis',
  dashboardHeader(
    title = img(src = 'PISAanalyzer.jpg', height = 100, align = 'left'),
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 350,
    
    fileInput('upload.file',
              'Please upload a tsv file with your data',
              accept = c('text/tsv', '.tsv')),
    
    # reference sample selection
    selectInput('reference.condition',
                'Please specify the reference condition for statistical analysise',
                ref.condition,
                multiple = FALSE),
    
    sliderInput('p.value.cutoff',
                'Specify a p-value which should be used as a cutoff for significance',
                min = 0.001, max = 0.2, value = 0.05, step = 0.001),
    
    sliderInput('delta.Sm.cutoff',
                paste0('Specify a ', paste0('\u0394'), 'Sm value which should be used as a cutoff of solubility alteration'),
                min = 0.05, max = 3, value = 0.3, step = 0.05),
   
    uiOutput('user.input.1'),
    uiOutput('user.input.2'),
    uiOutput('user.input.3'),
    uiOutput('user.input.4'),
    
    uiOutput('action'),
    
    p(strong(uiOutput('citation')))
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                                .skin-blue .main-header .navbar {
                                background-color: #870052;
                                }
                                /* logo */
                                .skin-blue .main-header .logo {
                                background-color: #ffffff ;
                                }
                                /* logo when hovered */
                                .skin-blue .main-header .logo:hover {
                                background-color: #ffffff ;
                                }
                                /* main sidebar */
                                .skin-blue .main-sidebar {
                                background-color: #808080;
                                }                         
                                /* active selected tab in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                                background-color: #808080;
                                }                          
                                /* other links in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                background-color: #808080;
                                }
                                /* other links in the sidebarmenu when hovered */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                                background-color: #808080;
                                }
                                /* toggle button when hovered  */                    
                                .skin-blue .main-header .navbar .sidebar-toggle:hover{
                                background-color: #870052;
                                }
                                .main-header { max-height: 100px; 
                                font-size:24px; 
                                font-weight:bold; 
                                line-height:24px;
                                }
                                .main-header .logo {
                                height: 100px;
                                font-size:24px; 
                                font-weight:bold; 
                                line-height:24px;align:
                                }
                                .skin-blue .sidebar a {
                                 color: #444;
                                 }
                                .main-sidebar {
                                float:top; margin-top:40px; padding-left:15px; padding-right:15px
                                }
                                '))),
    
    #result section      
    tabsetPanel(type = 'tabs',
                tabPanel('PISA analysis',
                         h1(''),
                         
                         plotlyOutput('volcano.plot') %>% withSpinner(),
                         h1(''),
                         splitLayout(uiOutput('download.volcano.source.data'),
                                     uiOutput('download.volcano.plot.pdf'),
                                     uiOutput('download.volcano.plot.svg')),
                         
                         
                         h1(''),
                         h1(''),
                         plotlyOutput('comparison.plot') %>% withSpinner(),
                         h1(''),
                         splitLayout(uiOutput('download.scatter.source.data'),
                                     uiOutput('download.scatter.plot.pdf'),
                                     uiOutput('download.scatter.plot.svg')),
                         
                         h1(''),
                         h1(''),
                         plotOutput('score.plot', width = '80%'),
                         h1(''),
                         splitLayout(uiOutput('download.score.source.data'),
                                     uiOutput('download.score.plot.pdf'),
                                     uiOutput('download.score.plot.svg')),
                         ),
                
                tabPanel('Instructions',
                         img(src = 'PISAanalyzer_instructions.jpg',
                             height = 1500,
                             align = 'middle')))
  )
)