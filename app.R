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
library(ggplot2)
library(ggthemes)

asthmaICER <- function (pGA=0.25,
                        pExacNoTxNoGA=0,
                        pExacTxNoGA=0,
                        pExacNoTxGA =0.55,
                        pExacTxGA = 0.05,
                        c_tx=2.49*2*30,
                        cExacER=575,
                        cExacNoHosp=126,
                        cGeneTest=199*(3796/7212)*(125.9/120.5),
                        uExacModAlpha=0.51,
                        uExacModBeta=0.38,
                        uExacERAlpha=0.36,
                        uExacERBeta=0.45,
                        # uExacHospAlpha=0.15,
                        # uExacHospBeta=0.30,
                        wtp=50000
) {

  p.notGA=1-pGA
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


  p.exa.Notr.GA.SD<-pExacNoTxGA*(0.2/1.96)

  p.exa.withtr.GA.SD<- pExacTxGA*(0.2/1.96)

  #############

  cExacMean<-0.4*cExacER+0.6*cExacNoHosp
  cost_exacSD=(0.2/1.96)* cExacMean



  set.seed(135)
  for (k in 1:1000){
    #########

    p.exa.Notr.GA.Par<-BetaPar(pExacNoTxGA,p.exa.Notr.GA.SD)
    p.exa.Notr.GA<-rbeta(1,p.exa.Notr.GA.Par$alpha,p.exa.Notr.GA.Par$beta)

    ###########


    p.exa.withtr.GA.Par<-BetaPar(pExacTxGA,p.exa.withtr.GA.SD)
    p.exa.withtr.GA<- rbeta(1,p.exa.withtr.GA.Par$alpha,p.exa.withtr.GA.Par$beta)
    ##############


    uExacModerate <-rbeta(1,uExacModAlpha, uExacModBeta)
    uExacER  <-rbeta(1,uExacERAlpha, uExacERBeta)
    qLossExac  <- (1-(uExacModerate*0.6+uExacER*0.4))*15/365
    qLossNoExac  <- 0*15/365

    ############

    cost_exacpar<-gammaPar(cExacMean, cost_exacSD)
    cost_exac<-rgamma(1,cost_exacpar$a,cost_exacpar$b)

    c.1 <- cost_exac
    c.2 <- 0
    c.3 <- cost_exac
    c.4 <-  0
    c.5 <- cost_exac+ c_tx+cGeneTest
    c.6 <-c_tx+cGeneTest
    c.7 <- cost_exac+ cGeneTest
    c.8 <-cGeneTest
    c.9 <-cost_exac+ c_tx
    c.10 <-c_tx
    c.11 <-cost_exac+ c_tx
    c.12 <- c_tx


    #################### ## ########################################

    C.notreatment<-c.1*p.exa.Notr.GA*0.25+c.2 *(1-p.exa.Notr.GA)*0.25+c.3*pExacNoTxNoGA*0.75+c.4*(1-pExacNoTxNoGA)*0.75
    C.OnlyGAs<-c.5*p.exa.withtr.GA*0.25 +c.6*(1-p.exa.withtr.GA)*0.25+c.7*pExacNoTxNoGA*0.75 +c.8*(1-pExacNoTxNoGA)*0.75
    C.all<-c.9*p.exa.withtr.GA*0.25+c.10*(1-p.exa.withtr.GA)*0.25+c.11*pExacTxNoGA*0.75+c.12*(1-pExacTxNoGA)*0.75

    q.notreatment<-qLossExac*p.exa.Notr.GA*0.25+qLossNoExac *(1-p.exa.Notr.GA)*0.25+qLossExac*pExacNoTxNoGA*0.75+qLossNoExac*(1-pExacNoTxNoGA)*0.75
    q.OnlyGAs<-qLossExac*p.exa.withtr.GA*0.25 +qLossNoExac*(1-p.exa.withtr.GA)*0.25+qLossExac*pExacNoTxNoGA*0.75 +qLossNoExac *(1-pExacNoTxNoGA)*0.75
    q.all<-qLossExac*p.exa.withtr.GA*0.25+qLossNoExac*(1-p.exa.withtr.GA)*0.25+qLossExac*pExacTxNoGA*0.75+qLossNoExac*(1-pExacTxNoGA)*0.75

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


wtpPlot <- function(res) {
   df <- tibble(wtp = seq(10000, 200000, by=5000), wtp_met=NA)
   for (i in 1:dim(df)[1]) {

     df[i, 2] <- wtpProb(res, as.numeric(df[i,1]))
   }

  p <- ggplot(df) +
        geom_line(aes(y=100*wtp_met, x=wtp)) +
        ylab("Probability of Being Cost-Effective") +
        xlab("Willingness-to-Pay Threhold") +
        theme_bw()

  return(p)

}

# Define UI for application that draws a histogram
ui <- fluidPage(
#  theme = shinytheme("united"),
    # Application title
    titlePanel("The economics of precision health: preventing air pollution-induced exacerbation in asthma"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          # numericInput("pGA",
          #              "Probability of Genetic Abnormality",
          #              min = 0,
          #              max = 1,
          #              value = 0.25,
          #              step = 0.01),

          numericInput("pExacNoTxNoGA",
                       "Probability of Exacerbations without Tx -  No Genetic Abnormality",
                       min = 0,
                       max = 1,
                       value = 0,
                       step = 0.01) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Assumption"),

          numericInput("pExacTxNoGA",
                       "Probability of Exacerbations with Tx - No Genetic Abnormality",
                       min = 0,
                       max = 1,
                       value = 0,
                       step = 0.01) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Assumption"),
          numericInput ("pExacNoTxGA",
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
          numericInput ("pExacTxGA",
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
          numericInput ("c_tx",
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
          numericInput ("cExacER",
                        "Cost of Exacerbation requiring ER Visit",
                        min = 0,
                        max = 1000,
                        value = 575,
                        step = 10
                        ) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Bielinski SJ, St. Sauver JL, Olson JE, et al. Are patients willing to incur out of pocket costs for pharmacogenomic
testing? Pharmacogenomics J 2017; 17: 1–3."),

          numericInput ("cExacNoHosp",
                        "Cost of Exacerbation Without Hospitalization",
                        min = 0,
                        max = 1000,
                        value = 126,
                        step = 10
          ) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Bielinski SJ, St. Sauver JL, Olson JE, et al. Are patients willing to incur out of pocket costs for pharmacogenomic
          testing? Pharmacogenomics J 2017; 17: 1–3."),

          numericInput ("cGeneTest",
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
          p(strong("Utilities - Exacerbation without hospitalization (Beta Distribution)")),
          fluidRow(
          column(4,
                   numericInput("uExacModAlpha",
                       "Alpha",
                       min = 0,
                       max = 100,
                       value = 0.51,
                       step = 1)),
          column(4,
                 numericInput("uExacModBeta",
                              "Beta",
                              min = 0,
                              max = 100,
                              value = 0.38,
                              step = 1))
          ),

        p(strong("Utilities - Exacerbation requiring ER visit  (Beta Distribution)")),
        fluidRow(
          column(4,
                 numericInput("uExacERAlpha",
                              "Alpha",
                              min = 0,
                              max = 100,
                              value = 0.36,
                              step = 1)),
          column(4,
                 numericInput("uExacERBeta",
                              "Beta",
                              min = 0,
                              max = 100,
                              value = 0.45,
                              step = 1))
        ),

        # p(strong("Utilities - Exacerbation requiring ER visit (Beta Distribution)")),
        # fluidRow(
        #   column(4,
        #          numericInput("uExacHospAlpha",
        #                       "Alpha",
        #                       min = 0,
        #                       max = 100,
        #                       value = 0.15,
        #                       step = 1)),
        #   column(4,
        #          numericInput("uExacHospBeta",
        #                       "Beta",
        #                       min = 0,
        #                       max = 100,
        #                       value = 0.30,
        #                       step = 1))
        # ),


          numericInput("wtp",
                       "Willingness to Pay",
                       min = 0,
                       max = 500000,
                       value = 50000,
                       step = 10000)
        ),



#uExacModAlpha
        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(type="tabs",
                      tabPanel("ICER",
                               tableOutput("ICER"),
                               textOutput("wtpProb"),
                               plotOutput("acceptability")
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
       asthmaICER(pGA               = input$pGA,
                  pExacNoTxNoGA     = input$pExacNoTxNoGA,
                  pExacTxNoGA       = input$pExacTxNoGA,
                  pExacNoTxGA       = input$pExacNoTxGA,
                  pExacTxGA         = input$pExacTxGA,
                  c_tx              = input$c_tx,
                  cExacER           = input$cExacER,
                  cExacNoHosp       = input$cExacNoHosp,
                  cGeneTest         = input$cGeneTest,
                  uExacModAlpha     = input$uExacModAlpha,
                  uExacModBeta      = input$uExacModBeta ,
                  uExacERAlpha      = input$uExacERAlpha ,
                  uExacERBeta       = input$uExacERBeta,
                  wtp               = input$wtp
       ))
    }, digits = 5)

    output$wtpProb <- renderText({

      paste0("At a willingess to pay of $",
      input$wtp,
      " the probability of the targeted intervention being cost-effective is ",
      100*wtpProb(asthmaICER(pGA               = input$pGA,
                             pExacNoTxNoGA     = input$pExacNoTxNoGA,
                             pExacTxNoGA       = input$pExacTxNoGA,
                             pExacNoTxGA       = input$pExacNoTxGA,
                             pExacTxGA         = input$pExacTxGA,
                             c_tx              = input$c_tx,
                             cExacER           = input$cExacER,
                             cExacNoHosp       = input$cExacNoHosp,
                             cGeneTest         = input$cGeneTest,
                             uExacModAlpha     = input$uExacModAlpha,
                             uExacModBeta      = input$uExacModBeta ,
                             uExacERAlpha      = input$uExacERAlpha ,
                             uExacERBeta       = input$uExacERBeta,
                             wtp               = input$wtp
                             ),
              wtp=input$wtp), "%")

    })

    output$acceptability <- renderPlot({

      wtpPlot(res = asthmaICER(pGA               = input$pGA,
                               pExacNoTxNoGA     = input$pExacNoTxNoGA,
                               pExacTxNoGA       = input$pExacTxNoGA,
                               pExacNoTxGA       = input$pExacNoTxGA,
                               pExacTxGA         = input$pExacTxGA,
                               c_tx              = input$c_tx,
                               cExacER           = input$cExacER,
                               cExacNoHosp       = input$cExacNoHosp,
                               cGeneTest         = input$cGeneTest,
                               uExacModAlpha     = input$uExacModAlpha,
                               uExacModBeta      = input$uExacModBeta ,
                               uExacERAlpha      = input$uExacERAlpha ,
                               uExacERBeta       = input$uExacERBeta,
                               wtp               = input$wtp
      ))
    })
}



# Run the application
shinyApp(ui = ui, server = server)
