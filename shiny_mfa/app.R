#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(ggplot2)
library(MFAg)
library(shiny)
library(png)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(

  # Application title
  titlePanel("Multiple Factor Analysis for Wine Tasting"),

  # Sidebars for input
  sidebarLayout(
    sidebarPanel(

      selectInput("select", label = h3("Plot Options"),
                  choices = list("Eigenvalue Barplot" = 1, "Common factor scores" = 2,
                                 "Partial factors scores" = 3, "Loadings" = 4),
                  selected = 1),

      sliderInput(inputId = "table",
                  label = "Which Wine Critic?",
                  min = 1,
                  max = 10,
                  value = 1)


    ),


    # Show plots
    mainPanel(
      plotOutput("plot")
    )
  )
))

# Define server
server <- shinyServer(function(input, output) {

  wines <- read.csv("wines.csv")

  # Add row names
  row.names(wines) <- wines[,1]

  # Separate the grand table into individual tables (create the set list)
  col_ind = grep("V2", colnames(wines))
  sets = list()
  for (i in 1:(length(col_ind))){
    # First 9 tables
    if (i < 10){
      sets[[i]]= (col_ind[i]-1):(col_ind[i+1]-2)
    }
    # Last table
    else{
      sets[[i]]= (col_ind[i]-1):(col_ind[i]+2)
    }
  }

  col_ind = grep("V15", colnames(wines))
  colnames(wines)[col_ind] <- "Peach"
  col_ind = grep("V14", colnames(wines))
  colnames(wines)[col_ind] <- "Grass"
  col_ind = grep("V13", colnames(wines))
  colnames(wines)[col_ind] <- "Melon"
  col_ind = grep("V12", colnames(wines))
  colnames(wines)[col_ind] <- "Hay"
  col_ind = grep("V11", colnames(wines))
  colnames(wines)[col_ind] <- "Vegetal"
  col_ind = grep("V10", colnames(wines))
  colnames(wines)[col_ind] <- "Flinty"
  col_ind = grep("V9", colnames(wines))
  colnames(wines)[col_ind] <- "Grassy"
  col_ind = grep("V8", colnames(wines))
  colnames(wines)[col_ind] <- "Leafy"
  col_ind = grep("V7", colnames(wines))
  colnames(wines)[col_ind] <- "Tropical"
  col_ind = grep("V6", colnames(wines))
  colnames(wines)[col_ind] <- "Citrus"
  col_ind = grep("V5", colnames(wines))
  colnames(wines)[col_ind] <- "Smoky"
  col_ind = grep("V4", colnames(wines))
  colnames(wines)[col_ind] <- "Mineral"
  col_ind = grep("V3", colnames(wines))
  colnames(wines)[col_ind] <- "Green Pepper"
  col_ind = grep("V2", colnames(wines))
  colnames(wines)[col_ind] <- "Passion Fruit"
  col_ind = grep("V1", colnames(wines))
  colnames(wines)[col_ind] <- "Cat Pee"

  # Defin MFA in S4 Class
  setClass(
    Class = "mfa",
    representation = representation(
      eigenvalues = "numeric",
      cfs = "matrix",
      pfs = "list",
      mol = "list",
      alpha = "numeric"
    ),

    # It seems that we dont need to valid anything

    # validity = function(object){
    #   }

    # A simple prototype
    prototype = prototype(
      eigenvalues = c(0,0,0,0),
      cfs = matrix(0, 10, 10),
      pfs = list(matrix(0, 10, 10)),
      mol = list(matrix(0, 10, 10)),
      alpha = c(0,0,0,0)

    )

  )


  # mfa constructor function
  mfa_const <- function(data, sets, ncomps = NULL, center = TRUE, scale = TRUE){

    # Separate the grand table into small tables, and record the names of table
    table_name = c()
    for (i in 1:(length(sets))){
      assign(paste("x",i ,sep=""), data[,sets[[i]]])
      # record the names of individual tables
      table_name = append(table_name, paste("x",i ,sep=""))
    }


    # Normalize all the data set, and combine the normalized data set
    alpha = c()
    combined_X = c()
    normtable_name = c()
    for (i in table_name){
      # Center the data set
      tmp = scale(eval(as.name(paste(i))) , center = center,  scale = F)

      #Scale the data set (depend on the scale input)
      if (scale){
        # compute the scale array
        scale_array = sqrt(apply(tmp^2,2,sum))
        # Scale the nomalized matrix
        norm_x = scale(tmp, center = F, scale = scale_array)
      }
      else{
        norm_x = tmp
      }


      # alpha weight vector
      alpha = append(alpha, rep(1,length(svd(norm_x)$d)) * (1/(svd(norm_x)$d[1])^2) )

      # record the normalized individual table
      assign(paste("norm",i ,sep=""), norm_x)
      normtable_name = append(normtable_name, paste("norm",i ,sep=""))

      # combined normalized matrx
      combined_X = as.matrix(cbind(combined_X,norm_x))
      # assign(paste("SVD",i ,sep=""), svd(norm_x))
    }

    # Create matrix for A
    A = diag(alpha, length(alpha),length(alpha))
    # Create vector for M
    n = length(data[,1])
    # M = diag(rep(1/n, n),n, n)
    M = rep(1/n, n)

    # GSVD!!
    GSVD = GSVD(combined_X, PLin = M, PCol = alpha)

    # Eigen values
    eigenvalues = (GSVD$d)^2

    # eigen value in matrix form
    d_mat = diag(GSVD$d , length(GSVD$d), length(GSVD$d))

    # Factor scores
    f = GSVD$u %*% d_mat

    # Q
    Q = GSVD$v

    # Partial factor scires
    K = length(table_name)
    tmp = 1
    pfs = list()
    mol = list()
    for (i in 1:K){
      x = eval(as.name(paste(normtable_name[i])))
      Q_tmp = Q[tmp:(tmp+length(x[1,])-1),]
      pfs_tmp = K*alpha[tmp]*x%*%Q_tmp
      tmp = tmp+length(x[1,])
      mol[[i]] = Q_tmp
      pfs[[i]] = pfs_tmp

    }



    # If ncomps is specified
    if (!is.null(ncomps)){
      eigenvalues = eigenvalues[1:ncomps]
      f = f[,1:ncomps]
      # pfs and mol is a list, so for loop it...
      for (i in 1:K){
        pfs[[i]] = pfs[[i]][,1:ncomps]
        mol[[i]] = mol[[i]][,1:ncomps]
      }
    }

    # Set the new class
    new(
      Class = "mfa",
      eigenvalues = eigenvalues,
      cfs = f,
      pfs = pfs ,
      mol = mol,
      alpha = alpha
    )

  }

  # Construct the object first
  wine <- mfa_const(data = wines, sets  = sets, ncomps = 11)

  # Extract features from the object
  eigenv = wine@eigenvalues
  cfs = wine@cfs
  pfs = wine@pfs
  Q = wine@mol
  
  # flag
  
  NZ <- readPNG("nz.png") 
  FR <- readPNG("fr.png")
  CA <- readPNG("ca.png")

  
  plot_pfs_shiny <- function(dimension1,
                             dimension2,
                             rownames_vec = as.character(1:length(dimension1))){

    dat <- data.frame(x = dimension1, y = dimension2, label = rownames_vec)

    ggplot(dat) +
      geom_point(aes(x,y)) +
      labs(title="Partial Factor Scores", x ="1", y = "2") +
      geom_text(data = dat, aes(x,y, label = rownames_vec), vjust = -0.5) +
      xlim(-2,2) +
      ylim(-2,2) + 
      mapply(function(xx, yy) 
        annotation_raster(NZ, xmin=xx-0.07, xmax=xx+0.07, ymin=yy-0.06, ymax=yy+0.06),
        dat$x[1:4], dat$y[1:4]) + 
      mapply(function(xx, yy) 
        annotation_raster(FR, xmin=xx-0.07, xmax=xx+0.07, ymin=yy-0.06, ymax=yy+0.06),
        dat$x[5:8], dat$y[5:8]) + 
      mapply(function(xx, yy) 
        annotation_raster(CA, xmin=xx-0.07, xmax=xx+0.07, ymin=yy-0.06, ymax=yy+0.06),
        dat$x[9:12], dat$y[9:12]
      ) 
      

  }
  
  plot_compromise_shiny <- function(dimension1,
                              dimension2,
                              rownames_vec = as.character(1:length(dimension1))){

    dat <- data.frame(x = dimension1, y = dimension2, label = rownames_vec)
    colnames(dat) <- c("x","y")

    ggplot(dat) +
      geom_point(aes(x,y)) +
      labs(title="Compromise of the tables", x ="1", y = "2") +
      geom_text(data = dat, aes(x,y, label = rownames_vec), vjust = -2) +
      xlim(-2,2) +
      ylim(-2,2) + 
      mapply(function(xx, yy) 
        annotation_raster(NZ, xmin=xx-0.07, xmax=xx+0.07, ymin=yy-0.06, ymax=yy+0.06),
        dat$x[1:4], dat$y[1:4]) + 
      mapply(function(xx, yy) 
        annotation_raster(FR, xmin=xx-0.07, xmax=xx+0.07, ymin=yy-0.06, ymax=yy+0.06),
        dat$x[5:8], dat$y[5:8]) + 
      mapply(function(xx, yy) 
        annotation_raster(CA, xmin=xx-0.07, xmax=xx+0.07, ymin=yy-0.06, ymax=yy+0.06),
        dat$x[9:12], dat$y[9:12]
      ) 
  }

  plot_vload <- function(dimension1,
                         dimension2,
                         rownames_vec = as.character(1:length(dimension1))){

    dat <- data.frame(x = dimension1, y = dimension2, label = rownames_vec)

    ggplot(dat) +
      geom_point(aes(x,y)) +
      labs(title="Variable loadings", x ="1", y = "2") +
      geom_text(data = dat, aes(x,y, label = rownames_vec), vjust = -0.5) +
      xlim(-0.5,0.5) +
      ylim(-0.6,0.6) 
    
       
  }

  # Make the plots

  pt <- reactive({
    input$select
    input$table
    eigenv_name = c(1,2,3,4,5,6,7,8,9,10,11)

    # Eigenvalue?
    if (input$select == 1){
      return(barplot(eigenv,
                     names.arg= eigenv_name,
                     ylab="Magnitude of Eigenvalues",
                     xlab="Number of component"))
    }

    # Compromise?
    else if(input$select == 2){
      cfs_dim1 = cfs[,1]
      cfs_dim2 = cfs[,2]
      return(plot_compromise_shiny(cfs_dim1,cfs_dim2,row.names(wines)))
    }

    # pfs?
    else if(input$select == 3){
      return(plot_pfs_shiny(pfs[[input$table]][,1], pfs[[input$table]][,2], rownames(wines)))
    }

    # Loading?
    else if(input$select == 4){
      return(plot_vload( Q[[input$table]][,1], Q[[input$table]][,2], colnames(wines)[sets[[input$table]]]))
    }

    else{
      return(NULL)
    }


  })

  output$plot = renderPlot({pt()})



})

# Run the application
shinyApp(ui = ui, server = server)

