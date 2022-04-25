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
library(purrr)
library(ftplottools)


asthmaICER <- function (pGA=0.25,
                        pExacNoTxNoGA=0,
                        pExacTxNoGA=0,
                        pExacNoTxGA =0.55,
                        pExacTxGA = 0.05,
                        c_tx=2.49*2*30,
                        cExacER=575,
                        cExacNoHosp=126,
                        cGeneTest=199*(3796/7212)*(125.9/120.5),
                        uExacMedAlpha=0.51,
                        uExacMedBeta=0.38,
                        uExacERAlpha=0.36,
                        uExacERBeta=0.25,
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
    alpha <- ((1 - mu) / (s^2) - 1/mu) * mu^2
    beta <- alpha * (1/mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }

  gammaPar <- function(m, s) {
    a <- (m^2)/(s^2)
    b <- m/(s^2)
    list (a=a,b=b)
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


    uExacMed <-rbeta(1,uExacMedAlpha, uExacMedBeta)

    uExacER  <-rbeta(1,uExacERAlpha, uExacERBeta)

    qLossExac  <- (1-(uExacMed*0.6+uExacER*0.4))*15/365
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

    C.notreatment<-c.1*p.exa.Notr.GA*pGA+c.2 *(1-p.exa.Notr.GA)*pGA+c.3*pExacNoTxNoGA*(1-pGA)+c.4*(1-pExacNoTxNoGA)*(1-pGA)
    C.OnlyGAs<-c.5*p.exa.withtr.GA*pGA +c.6*(1-p.exa.withtr.GA)*pGA+c.7*pExacNoTxNoGA*(1-pGA) +c.8*(1-pExacNoTxNoGA)*(1-pGA)
    C.all<-c.9*p.exa.withtr.GA*pGA+c.10*(1-p.exa.withtr.GA)*pGA+c.11*pExacTxNoGA*(1-pGA)+c.12*(1-pExacTxNoGA)*(1-pGA)

    q.notreatment<-qLossExac*p.exa.Notr.GA*pGA+qLossNoExac *(1-p.exa.Notr.GA)*pGA+qLossExac*pExacNoTxNoGA*(1-pGA)+qLossNoExac*(1-pExacNoTxNoGA)*(1-pGA)
    q.OnlyGAs<-qLossExac*p.exa.withtr.GA*pGA +qLossNoExac*(1-p.exa.withtr.GA)*pGA+qLossExac*pExacNoTxNoGA*(1-pGA) +qLossNoExac *(1-pExacNoTxNoGA)*(1-pGA)
    q.all<-qLossExac*p.exa.withtr.GA*pGA+qLossNoExac*(1-p.exa.withtr.GA)*pGA+qLossExac*pExacTxNoGA*(1-pGA)+qLossNoExac*(1-pExacTxNoGA)*(1-pGA)

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


  LE <- as.character(round(c(q.notreatment,q.OnlyGAs, q.all), 5))

  seqLE <-c("Ref.",
            round(q.OnlyGAs-q.notreatment, 5),
            round(q.all-q.OnlyGAs, 5))


  C  <- c(round(C.notreatment,2), round(C.OnlyGAs,2), round(C.all,2))
  names(LE) <- names(C) <- c("notreatment", "OnlyGAs", "q.all")

  seqC <- c("Ref.",
            round(C.OnlyGAs-C.notreatment, 2),
            round(C.all-C.OnlyGAs, 2))

  seqICER <- c("Ref.",
               round((C.OnlyGAs-C.notreatment)/(-q.OnlyGAs+q.notreatment),0),
               round((C.all-C.OnlyGAs)/(-q.all+q.OnlyGAs),0))

  icer <- tibble (Strategy=c("No Intervention", "Prevention only for GA", "Prevention for all"),
                  `Costs (CAD)`= C,
                  `QALY Loss` = LE,
                  `Sequential incremental costs (CAD)` = seqC,
                  `Sequential incremental QALY loss` = seqLE,
                  `Sequential ICER (CAD/QALY)` = seqICER
  )

  return(icer)

}

wtpProb <- function(wtp, res) {

  res_ICER_prob <- as_tibble(res)
  res_ICER_prob <- res_ICER_prob %>%
    mutate(ICER_OnlyGAs = (`only GA_cost`-notreatment_cost)/(-Q_onlyGAs + Q_notreatment)) %>%
    mutate(wtp_met = ICER_OnlyGAs<=wtp) %>%
    summarise_all(mean, na.rm=TRUE)

  wtp_p <- res_ICER_prob$wtp_met

  return(wtp_p)
}


wtpPlot <- function(res, wtpInput) {

  df <- tibble(wtp = seq(10000, 200000, by=5000), wtp_met=NA) %>%
     mutate (wtp_met=100*map_dbl(wtp, wtpProb, res=res))

  p <- ggplot(df, aes(y=wtp_met, x=wtp)) +
        geom_line(size=1.2, color="coral2") +
        ylab("Probability of Being Cost-Effective") +
        xlab("Willingness-to-Pay Threhold") +
        geom_vline(xintercept = wtpInput,  colour="grey", linetype="dashed") +
        geom_text(aes(x=wtpInput, label="\n Willingness to pay", y=20), check_overlap = TRUE, angle=90) +
        scale_colour_brewer(palette = "Dark2") +
        ft_theme() +
        ggtitle("Cost-Effectiveness Acceptability Curve") +
        theme(text = element_text(size=12)) +
        theme(legend.position = "none") +
        theme(legend.title=element_blank())

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
           sliderInput("pGA",
                        "Prevalence of Genetic Abnormality (%)",
                        min = 0,
                        max = 1,
                        value = 0.25,
                        step = 0.01) %>%
             helper(icon = "question-circle",
                    colour = "black",
                    type = "inline",
                    content = "Genetic abnormality is defined as either GSTT1 null, GSTM1
null or GSTP1 Ile105."),
          p(strong("Risk of additional exacerbations in asthmatics without preventive intervention (%)")),
          fluidRow(
            column(6,tags$head(tags$style(HTML(".not_bold label {font-weight:normal;}"))),
                   div(sliderInput("pExacNoTxNoGA",
                                "No Abnormality",
                                min = 0,
                                max = 1,
                                value = 0,
                                step = 0.01) %>%
                     helper(icon = "question-circle",
                            colour = "black",
                            type = "inline",
                            content = "Genetic abnormality is defined as either GSTT1 null, GSTM1
null or GSTP1 Ile105.")), class="not_bold"),
            column(6,
                   div(sliderInput("pExacNoTxGA",
                                "Genetic Abnormality",
                                min = 0.01,
                                max = 1,
                                value = 0.55,
                                step = 0.01) %>%
                     helper(icon = "question-circle",
                            colour = "black",
                            type = "inline",
                            content = c("Genetic abnormality is defined as either GSTT1 null, GSTM1
null or GSTP1 Ile105. Reference: Orellano P, Quaranta N, Reynoso J, et al. Effect of outdoor air pollution on asthma exacerbations in children and
adults: systematic review and multilevel meta-analysis. PLoS One 2017; 12: e0174050.
23", "Zafari Z, Sadatsafavi M, Marra CA, et al. Cost-effectiveness of bronchial thermoplasty, omalizumab, and standard
therapy for moderate-to-severe allergic asthma. PLoS One 2016; 11: e0146003."))
          )), class="not_bold"),

         sliderInput("TxEffect",
                      "Preventive intervention will cut the risk by a factor of:",
                      min = 1,
                      max = 100,
                      value = 0.55/0.05,
                      step = 1,
                      ticks=F) %>%
  helper(icon = "question-circle",
         colour = "black",
         type = "inline",
         content = c("Reference: Orellano P, Quaranta N, Reynoso J, et al. Effect of outdoor air pollution on asthma exacerbations in children and
adults: systematic review and multilevel meta-analysis. PLoS One 2017; 12: e0174050.
23", "Zafari Z, Sadatsafavi M, Marra CA, et al. Cost-effectiveness of bronchial thermoplasty, omalizumab, and standard
therapy for moderate-to-severe allergic asthma. PLoS One 2016; 11: e0146003.")),


p(strong("Risk of additional exacerbations in asthmatics with preventive intervention (%)")) %>%
  helper(icon = "question-circle",
         colour = "black",
         type = "inline",
         content = c("Probabality of exacerbation without treatment divided by treatment risk ratio provided above.")),
fluidRow(
  column(6,
         p("No Genetic Abnormality"),
         textOutput("pExacTxNoGA")),
  column(6,
         p("Genetic Abnormality"),
         textOutput("pExacTxGA")
  )), p(),

          numericInput ("c_tx",
                        "Cost of Preventive Treatment (2018 CAD)",
                        min = 0,
                        max = 1000,
                        value = 149.4,
                        step = 10
                          ) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Advair 500/50 twice per day, per month. Reference: Minelli C, Granell R, Newson R, et al. Glutathione-S-transferase genes and asthma phenotypes: a Human Genome
Epidemiology (HuGE) systematic review and meta-analysis including unpublished data. Int J Epidemiol 2010; 39:
539–562."),


p(strong("Cost of Exacerbations (2018 CAD)")),
fluidRow(
  column(6,
         numericInput("cExacER",
                      "Requiring ER Visit",
                      min = 0,
                      max = 1000,
                      value = 575,
                      step = 10) %>%
           helper(icon = "question-circle",
                  colour = "black",
                  type = "inline",
                  content = "Reference: Bielinski SJ, St. Sauver JL, Olson JE, et al. Are patients willing to incur out of pocket costs for pharmacogenomic
testing? Pharmacogenomics J 2017; 17: 1–3.")),
  column(6,
         numericInput("cExacNoHosp",
                      "No Hospitalization",
                      min = 0,
                      max = 1000,
                      value = 126,
                      step = 10) %>%
           helper(icon = "question-circle",
                  colour = "black",
                  type = "inline",
                  content = c("Reference: Bielinski SJ, St. Sauver JL, Olson JE, et al. Are patients willing to incur out of pocket costs for pharmacogenomic
testing? Pharmacogenomics J 2017; 17: 1–3."))
  )),



          numericInput ("cGeneTest",
                        "Cost for Genetic Test (2018 CAD)",
                        min = 0,
                        max = 1000,
                        value = 109.44,
                        step = 10) %>%
            helper(icon = "question-circle",
                   colour = "black",
                   type = "inline",
                   content = "Reference: Campbell JD, Spackman DE, Sullivan SD. The costs and consequences of omalizumab in uncontrolled asthma
from a USA payer perspective. Allergy 2010; 65: 1141–1148"),
        #   p(strong("Utilities - Exacerbation without hospitalization (Beta Distribution)")),
        #   fluidRow(
        #   column(4,
        #            numericInput("uExacMedAlpha",
        #                "Alpha",
        #                min = 0,
        #                max = 100,
        #                value = 0.51,
        #                step = 1)),
        #   column(4,
        #          numericInput("uExacMedBeta",
        #                       "Beta",
        #                       min = 0,
        #                       max = 100,
        #                       value = 0.38,
        #                       step = 1))
        #   ),
        #
        # p(strong("Utilities - Exacerbation requiring ER visit (Beta Distribution)")),
        # fluidRow(
        #   column(4,
        #          numericInput("uExacERAlpha",
        #                       "Alpha",
        #                       min = 0,
        #                       max = 100,
        #                       value = 0.36,
        #                       step = 1)),
        #   column(4,
        #          numericInput("uExacERBeta",
        #                       "Beta",
        #                       min = 0,
        #                       max = 100,
        #                       value = 0.45,
        #                       step = 1))
        # ),

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
                       "Willingness to Pay (2018 CAD)",
                       min = 0,
                       max = 500000,
                       value = 50000,
                       step = 10000)
        ),



#uExacMedAlpha
        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(type="tabs",
                      tabPanel("ICER",
                               br(),
                               div(id = "background", includeMarkdown("./background.rmd")),
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

  output$pExacTxNoGA <- renderText({
     as.character(input$pExacNoTxNoGA/input$TxEffect)
    })

    output$pExacTxGA <- renderText({
      as.character(input$pExacNoTxGA/input$TxEffect)
    })

    output$ICER <- renderTable({
      sequentialICER(
       asthmaICER(pGA               = input$pGA,
                  pExacNoTxNoGA     = input$pExacNoTxNoGA,
                  pExacTxNoGA       = input$pExacNoTxNoGA/input$TxEffect,
                  pExacNoTxGA       = input$pExacNoTxGA,
                  pExacTxGA         = input$pExacNoTxGA/input$TxEffect,
                  c_tx              = input$c_tx,
                  cExacER           = input$cExacER,
                  cExacNoHosp       = input$cExacNoHosp,
                  cGeneTest         = input$cGeneTest,
                  # uExacMedAlpha     = input$uExacMedAlpha,
                  # uExacMedBeta      = input$uExacMedBeta ,
                  # uExacERAlpha      = input$uExacERAlpha ,
                  # uExacERBeta       = input$uExacERBeta,
                  wtp               = input$wtp
       ))
    })

    output$wtpProb <- renderText({

      paste0("At a willingess to pay of $",
      input$wtp,
      " the probability of the targeted intervention being cost-effective is ",
      100*wtpProb(res=asthmaICER(pGA               = input$pGA,
                                 pExacNoTxNoGA     = input$pExacNoTxNoGA,
                                 pExacTxNoGA       = input$pExacNoTxNoGA/input$TxEffect,
                                 pExacNoTxGA       = input$pExacNoTxGA,
                                 pExacTxGA         = input$pExacNoTxGA/input$TxEffect,
                                 c_tx              = input$c_tx,
                                 cExacER           = input$cExacER,
                                 cExacNoHosp       = input$cExacNoHosp,
                                 cGeneTest         = input$cGeneTest,
                                 # uExacMedAlpha     = input$uExacMedAlpha,
                                 # uExacMedBeta      = input$uExacMedBeta ,
                                 # uExacERAlpha      = input$uExacERAlpha ,
                                 # uExacERBeta       = input$uExacERBeta,
                                 wtp               = input$wtp
                             ),
              wtp=input$wtp), "%")

    })

    output$acceptability <- renderPlot({


      p <- wtpPlot(res = asthmaICER(pGA               = input$pGA,
                                    pExacNoTxNoGA     = input$pExacNoTxNoGA,
                                    pExacTxNoGA       = input$pExacNoTxNoGA/input$TxEffect,
                                    pExacNoTxGA       = input$pExacNoTxGA,
                                    pExacTxGA         = input$pExacNoTxGA/input$TxEffect,
                                    c_tx              = input$c_tx,
                                    cExacER           = input$cExacER,
                                    cExacNoHosp       = input$cExacNoHosp,
                                    cGeneTest         = input$cGeneTest,
                                    # uExacMedAlpha     = input$uExacMedAlpha,
                                    # uExacMedBeta      = input$uExacMedBeta ,
                                    # uExacERAlpha      = input$uExacERAlpha ,
                                    # uExacERBeta       = input$uExacERBeta,
                                    wtp               = input$wtp
      ),
      wtpInput=input$wtp)
      p
    })
}



# Run the application
shinyApp(ui = ui, server = server)
