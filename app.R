#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyhelper)
library(tibble)
library(markdown)
library(dplyr)


asthmaICER <- function (p.GA=0.25,
                        p.exa.Notr.notGA=0,
                        p.exa.withtr.notGA=0,
                        p.exa.Notr.GA.mean =0.55,
                        p.exa.withtr.GA.mean = 0.05,
                        Cost_treatemnt=2.49*2*30,
                        cost_exacmean=0.4*575+0.6*126,
                        Cost_genetic_id=199*(3796/7212)*(125.9/120.5),
                        wtp=50000
) {

  p.notGA=1-p.GA
  ####################### Input Model Parameters   ############################################
  res=matrix(NA,nrow=1000,ncol=7)
  colnames(res)<-c("run","notreatment_cost","only GA_cost","All_cost","Q_notreatment",
                   "Q_onlyGAs","Q_all")
  Strategies <- c("No treatment", "OnlyGAs", "all")

  BetaPar <- function(mu, s) {
    alpha <- ((1 - mu) / (s^2) - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }

  gammaPar <- function(m, s ) {
    a <- (m^2 ) / ( s ^2 )
    b <- m/ ( s ^2 )
    list ( a=a , b=b )
  }


  p.exa.Notr.GA.SD<-p.exa.Notr.GA.mean*(0.2/1.96)

  p.exa.withtr.GA.SD<- p.exa.withtr.GA.mean*(0.2/1.96)

  #############

  cost_exacSD=(0.2/1.96)* cost_exacmean



  set.seed(135)
  for (k in 1:1000){
    #########

    p.exa.Notr.GA.Par<-BetaPar(p.exa.Notr.GA.mean,p.exa.Notr.GA.SD)
    p.exa.Notr.GA<-rbeta(1,p.exa.Notr.GA.Par$alpha,p.exa.Notr.GA.Par$beta)

    ###########


    p.exa.withtr.GA.Par<-BetaPar(p.exa.withtr.GA.mean,p.exa.withtr.GA.SD)
    p.exa.withtr.GA<- rbeta(1,p.exa.withtr.GA.Par$alpha,p.exa.withtr.GA.Par$beta)
    ##############


    u_ex_med<-rbeta(1,0.51, 0.38)
    u_ex_ER<-rbeta(1,0.36, 0.25)
    qloss_exca  <- (1-(u_ex_med*0.6+u_ex_ER*0.4))*15/365
    qloss_notexca  <- 0*15/365

    ############

    cost_exacpar<-gammaPar(cost_exacmean, cost_exacSD)
    cost_exac<-rgamma(1,cost_exacpar$a,cost_exacpar$b)

    c.1 <- cost_exac
    c.2 <- 0
    c.3 <- cost_exac
    c.4 <-  0
    c.5 <- cost_exac+ Cost_treatemnt+Cost_genetic_id
    c.6 <-Cost_treatemnt+Cost_genetic_id
    c.7 <- cost_exac+ Cost_genetic_id
    c.8 <-Cost_genetic_id
    c.9 <-cost_exac+ Cost_treatemnt
    c.10 <-Cost_treatemnt
    c.11 <-cost_exac+ Cost_treatemnt
    c.12 <- Cost_treatemnt


    #################### ## ########################################

    C.notreatment<-c.1*p.exa.Notr.GA*0.25+c.2 *(1-p.exa.Notr.GA)*0.25+c.3*p.exa.Notr.notGA*0.75+c.4*(1-p.exa.Notr.notGA)*0.75
    C.OnlyGAs<-c.5*p.exa.withtr.GA*0.25 +c.6*(1-p.exa.withtr.GA)*0.25+c.7*p.exa.Notr.notGA*0.75 +c.8*(1-p.exa.Notr.notGA)*0.75
    C.all<-c.9*p.exa.withtr.GA*0.25+c.10*(1-p.exa.withtr.GA)*0.25+c.11*p.exa.withtr.notGA*0.75+c.12*(1-p.exa.withtr.notGA)*0.75

    q.notreatment<-qloss_exca*p.exa.Notr.GA*0.25+qloss_notexca *(1-p.exa.Notr.GA)*0.25+qloss_exca*p.exa.Notr.notGA*0.75+qloss_notexca*(1-p.exa.Notr.notGA)*0.75
    q.OnlyGAs<-qloss_exca*p.exa.withtr.GA*0.25 +qloss_notexca*(1-p.exa.withtr.GA)*0.25+qloss_exca*p.exa.Notr.notGA*0.75 +qloss_notexca *(1-p.exa.Notr.notGA)*0.75
    q.all<-qloss_exca*p.exa.withtr.GA*0.25+qloss_notexca*(1-p.exa.withtr.GA)*0.25+qloss_exca*p.exa.withtr.notGA*0.75+qloss_notexca*(1-p.exa.withtr.notGA)*0.75

    res[k,]=c(k,C.notreatment,C.OnlyGAs,C.all,q.notreatment,q.OnlyGAs,q.all)
  }


    return(res)


}



sequentialICER <- function(res) {
  res2 <- as_tibble(res) %>% summarise_all(mean, na.rm=TRUE)
  q.notreatment <- res2$Q_notreatment
  q.OnlyGAs     <- res2$Q_onlyGAs
  q.all         <- res2$Q_all

  C.notreatment <- res2$notreatment_cost
  C.OnlyGAs     <- res2$`only GA_cost`
  C.all         <- res2$All_cost


  LE <- round(c(q.notreatment,q.OnlyGAs, q.all), 5)
  seqLE <-c("Ref.",
            round(q.OnlyGAs-q.notreatment, 5),
            round(q.all-q.OnlyGAs, 5))


  C  <- round(c(C.notreatment, C.OnlyGAs, C.all),2)
  names(LE) <- names(C) <- c("notreatment", "OnlyGAs", "q.all")

  seqC <- c("Ref.",
            round(C.OnlyGAs-C.notreatment, 2),
            round(C.all-C.OnlyGAs, 2))

  seqICER <- c("Ref.",
               round((C.OnlyGAs-C.notreatment)/(-q.OnlyGAs+q.notreatment),0),
               round((C.all-C.OnlyGAs)/(-q.all+q.OnlyGAs),0))

  icer <- tibble (Strategy=c("No Intervention", "Prevention only for GA", "Prevention for all"),
                  Costs = C,
                  `QALY Loss` = LE,
                  `Sequential incremental costs` = seqC,
                  `Sequential incremental QALY loss` = seqLE,
                  `Sequential ICER` = seqICER
  )

  return(icer)



}

wtpProb <- function(res, wtp) {
  # plotting
  res_ICER_prob <- as_tibble(res)
  res_ICER_prob <- res_ICER_prob %>%
    mutate(ICER_OnlyGAs = (`only GA_cost`-notreatment_cost)/(-Q_onlyGAs + Q_notreatment)) %>%
    mutate(wtp_met = ICER_OnlyGAs<=wtp) %>%
    summarise_all(mean, na.rm=TRUE)

  wtp_p <- res_ICER_prob$wtp_met

  return(wtp_p)
}


# wtpPlot <- function(res) {
#   df <- tibble(wtp = seq(10000, 200000, by=5000), wtp_met=NA)
#   for (i in 1:dim(df)[1]) {
#
#   }
#
#
#        )
# }

# Define UI for application that draws a histogram
ui <- fluidPage(
#  theme = shinytheme("united"),
    # Application title
    titlePanel("The economics of precision health: preventing air pollution-induced exacerbation in asthma"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          # numericInput("p.GA",
          #              "Probability of Genetic Abnormality",
          #              min = 0,
          #              max = 1,
          #              value = 0.25,
          #              step = 0.01),

          numericInput("p.exa.Notr.notGA",
                       "Probability of Exacerbations without Tx -  No Genetic Abnormality",
                       min = 0,
                       max = 1,
                       value = 0,
                       step = 0.01) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Assumption"),

          numericInput("p.exa.withtr.notGA",
                       "Probability of Exacerbations with Tx - No Genetic Abnormality",
                       min = 0,
                       max = 1,
                       value = 0,
                       step = 0.01) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Assumption"),
          numericInput ("p.exa.Notr.GA.mean",
                        "Probability of Exacerbation without Tx - Genetic Abnormality",
                        min = 0,
                        max = 1,
                        value = 0.55,
                        step = 0.01) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = c("Orellano P, Quaranta N, Reynoso J, et al. Effect of outdoor air pollution on asthma exacerbations in children and
adults: systematic review and multilevel meta-analysis. PLoS One 2017; 12: e0174050.
23", "Zafari Z, Sadatsafavi M, Marra CA, et al. Cost-effectiveness of bronchial thermoplasty, omalizumab, and standard
therapy for moderate-to-severe allergic asthma. PLoS One 2016; 11: e0146003.")),
          numericInput ("p.exa.withtr.GA.mean",
                        "Probability of Exacerbation with Tx - Genetic Abnormality",
                        min = 0,
                        max = 1,
                        value = 0.05,
                        step = 0.01) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = c("Orellano P, Quaranta N, Reynoso J, et al. Effect of outdoor air pollution on asthma exacerbations in children and
adults: systematic review and multilevel meta-analysis. PLoS One 2017; 12: e0174050.
23", "Zafari Z, Sadatsafavi M, Marra CA, et al. Cost-effectiveness of bronchial thermoplasty, omalizumab, and standard
therapy for moderate-to-severe allergic asthma. PLoS One 2016; 11: e0146003.")),
          numericInput ("Cost_treatemnt",
                        "Cost of Treatment",
                        min = 0,
                        max = 1000,
                        value = 149.4,
                        step = 10
                          ) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Minelli C, Granell R, Newson R, et al. Glutathione-S-transferase genes and asthma phenotypes: a Human Genome
Epidemiology (HuGE) systematic review and meta-analysis including unpublished data. Int J Epidemiol 2010; 39:
539–562."),
          numericInput ("cost_exacmean",
                        "Mean Exacerbation Cost",
                        min = 0,
                        max = 1000,
                        value = 305.6,
                        step = 10
                        ) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Bielinski SJ, St. Sauver JL, Olson JE, et al. Are patients willing to incur out of pocket costs for pharmacogenomic
testing? Pharmacogenomics J 2017; 17: 1–3."),

          numericInput ("Cost_genetic_id",
                        "Cost for Genetic Test",
                        min = 0,
                        max = 1000,
                        value = 109.44,
                        step = 10) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Campbell JD, Spackman DE, Sullivan SD. The costs and consequences of omalizumab in uncontrolled asthma
from a USA payer perspective. Allergy 2010; 65: 1141–1148"),

          numericInput("wtp",
                       "Willingness to Pay",
                       min = 0,
                       max = 500000,
                       value = 50000,
                       step = 10000)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(type="tabs",
                      tabPanel("ICER",
                               tableOutput("ICER"),
                               textOutput("wtpProb")
                      ),
                      tabPanel("Terms",  includeMarkdown("./disclaimer.Rmd")),
                      tabPanel("About",  includeMarkdown("./about.Rmd"))#,
                               #imageOutput("logos")

          )


        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  observe_helpers(help_dir = "helpfiles")

    output$ICER <- renderTable({
      sequentialICER(
       asthmaICER(p.GA                 = input$p.GA,
                  p.exa.Notr.notGA     = input$p.exa.Notr.notGA,
                  p.exa.withtr.notGA   = input$p.exa.withtr.notGA,
                  wtp                  = input$wtp,
                  p.exa.Notr.GA.mean   = input$p.exa.Notr.GA.mean   ,
                  p.exa.withtr.GA.mean = input$p.exa.withtr.GA.mean ,
                  Cost_treatemnt       = input$Cost_treatemnt       ,
                  cost_exacmean        = input$cost_exacmean        ,
                  Cost_genetic_id      = input$Cost_genetic_id       ))
    }, digits = 5)

    output$wtpProb <- renderText({

      paste0("At a willingess to pay of $",
      input$wtp,
      " the probability of the targeted intervention being cost-effective is ",
      100*wtpProb(asthmaICER(p.GA                 = input$p.GA,
                             p.exa.Notr.notGA     = input$p.exa.Notr.notGA,
                             p.exa.withtr.notGA   = input$p.exa.withtr.notGA,
                             wtp                  = input$wtp,
                             p.exa.Notr.GA.mean   = input$p.exa.Notr.GA.mean   ,
                             p.exa.withtr.GA.mean = input$p.exa.withtr.GA.mean ,
                             Cost_treatemnt       = input$Cost_treatemnt       ,
                             cost_exacmean        = input$cost_exacmean        ,
                             Cost_genetic_id      = input$Cost_genetic_id       ),
              wtp=input$wtp), "%")

    })
}



# Run the application
shinyApp(ui = ui, server = server)
