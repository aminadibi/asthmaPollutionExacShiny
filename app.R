#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(ggplot2)
library(ggthemes)
library(tibble)

asthmaICER <- function (p.GA=0.25,
                        p.exa.Notr.notGA=0,
                        p.exa.withtr.notGA=0,
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

  p.exa.Notr.GA.mean <- 0.55
  p.exa.Notr.GA.SD<-0.55*(0.2/1.96)

  p.exa.withtr.GA.mean <- 0.05
  p.exa.withtr.GA.SD<-0.05*(0.2/1.96)

  #############
  Cost_treatemnt <-2.49*2*30

  cost_exacmean<-0.4*575+0.6*126
  cost_exacSD=(0.2/1.96)* cost_exacmean

  Cost_genetic_id<-199*(3796/7212)*(125.9/120.5)


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
    LE <- round(c(q.notreatment,q.OnlyGAs, q.all), 5)
    seqLE <-c("Ref.",
                     round(q.OnlyGAs-q.notreatment, 5),
                     round(q.all-q.notreatment, 5))

    C  <- round(c(C.notreatment, C.OnlyGAs, C.all),2)
    names(LE) <- names(C) <- c("notreatment", "OnlyGAs", "q.all")

    seqC <- c("Ref.",
                     round(C.OnlyGAs-C.notreatment, 2),
                     round(C.all-C.notreatment, 2))

    seqICER <- c("Ref.",
                 round((C.OnlyGAs-C.notreatment)/(-q.OnlyGAs+q.notreatment),0),
                 round((C.all-C.notreatment)/(-q.all+q.notreatment),0))

    icer <- tibble (Strategy=c("No Intervention", "Prevention only for GA", "Prevention for all"),
                     Costs = C,
                    `QALY Loss` = LE,
                    `Sequential incremental costs` = seqC,
                    `Sequential incremental QALY loss` = seqLE,
                    `Sequential ICER` = seqICER
                    )


    return(icer)
  }

}


# Define UI for application that draws a histogram
ui <- fluidPage(
  #theme = shinytheme("united"),
    # Application title
    titlePanel("The economics of precision health: preventing air pollution-induced exacerbation in asthma"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          numericInput("p.GA",
                       "Probability of Genetic Abnormality",
                       min = 0,
                       max = 1,
                       value = 0.25),

          numericInput("p.exa.Notr.notGA",
                       "Probability of Exacerbations without Tx -  No Genetic Abnormality",
                       min = 0,
                       max = 1,
                       value = 0),

          numericInput("p.exa.withtr.notGA",
                       "Probability of Exacerbations with Tx - No Genetic Abnormality",
                       min = 0,
                       max = 1,
                       value = 0),

          numericInput("wtp",
                       "Willingness to Pay",
                       min = 0,
                       max = 500000,
                       value = 50000)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tableOutput("LE")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {


    output$LE<- renderTable({
       asthmaICER(p.GA               = input$p.GA,
                  p.exa.Notr.notGA   = input$p.exa.Notr.notGA,
                  p.exa.withtr.notGA = input$p.exa.withtr.notGA,
                  wtp                = input$wtp)
    }, digits = 5)
}

# Run the application
shinyApp(ui = ui, server = server)
