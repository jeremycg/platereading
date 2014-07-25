library(shiny)

startingdir=getwd()#change for your directory
setwd(startingdir)		
strainlist=read.csv("strainlist.csv")

# Define UI for dataset viewer application
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Plate Analysis"),

  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
	sliderInput("cutoff", 
                "cutoff of residuals:", 
                min = 0.0,
                max = 7.0, 
                value = 1.0,step=0.01,ticks=T),
	selectInput("variable", "Choose a variable:", 
                choices = c("mumax","lag","od0","odmax")),
	selectInput("query", "Strain:",
                levels(strainlist$strain))

    
  ),
mainPanel(
    tabsetPanel(
      tabPanel("BoxPlots", plotOutput("Plot1",height="100%")), 
      tabPanel("Table", tableOutput("table1"),tableOutput("table2")), 
      tabPanel("Histograms", plotOutput("Plot3",height="100%")),
	  tabPanel("Fitted Plots", plotOutput("Plotfitted",height="100%"))
    )
  )
))
