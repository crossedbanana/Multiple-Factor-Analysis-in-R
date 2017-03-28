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

#private function to check if dataset is matrix or dataframe
check_dataset <- function(dataset){

  if((class(dataset) != "matrix") & (class(dataset) != "data.frame")){

    stop("\n'dataset' must be a of class 'matrix' or of class 'data.frame'")

  }
  TRUE
}

#private function to check if sets is a list and is composed of numeric vectors
check_sets <- function(sets){

  if(class(sets) != "list"){

    stop("\n'sets' must be a of class 'list'")

  }

  for(i in 1:length(sets)){

    if((class(sets[[i]]) != "integer") &
       (class(sets[[i]]) != "numeric") &
       (class(sets[[i]]) != "character")){

      stop("\n'sets' must contain integer, numeric, or character vectors")

    }

  }
  TRUE
}

#private function to check if ncomps is an integer
check_ncomps <- function(ncomps){

  if((class(ncomps) != "integer") & (class(ncomps) != "numeric")){

    stop("\n'ncomps' must be an integer")

  }

  if(ncomps %% 1 != 0){

    stop("\n'ncomps' must be an integer")

  }
  TRUE
}

#private function to check if center is a logical or numeric vector
check_center <- function(center){

  if((class(center) != "logical") & (class(center) != "numeric")){

    stop("\n'center' must be a logical or numeric vector")

  }
  TRUE
}

#private function to check if scale is a logical or numeric vector
check_scale <- function(scale){

  if((class(scale) != "logical") & (class(scale) != "numeric")){

    stop("\n'scale' must be a logical or numeric vector")

  }
  TRUE
}

# A library for GSVD Computation
library(MFAg)
# public constructer function

#' @title mfa_const
#' @description Creates an object of class \code{"mfa"}
#' @param data data set (matrix or data frame)
#' @param sets list of vectors indicating the sets of variables (i.e. the blocks)
#' @param ncomps integer indicating how many number of components (i.e. factors) are to be extracted
#' @param center logical value
#' @param scale logical value
#' @return an object of class mfa
#' @export
#' @examples
#'
#' simulate_data <- matrix(rnorm(200), ncol = 10, byrow = TRUE)
#' test <- mfa_const(data = simulate_data, sets = list(1:3, 4:5, 6:10), ncomps = 2)
#'
mfa_const <- function(data, sets, ncomps = NULL, center = TRUE, scale = TRUE){

  check_dataset(data)
  check_sets(sets)
  check_ncomps(ncomps)
  check_center(center)
  check_scale(scale)

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


#print method
#' @export
setMethod(
  "print",
  signature = "mfa",
  function(x, ...) {
    cat('object "mfa"\n\n')
    cat("eigenvalues: ", x@eigenvalues, "\n\n")
    cat("cfs: \n")
    print(x@cfs)
  }
)

#Compromise of the 10 tables:

library(ggplot2)
library(png)
library(gridGraphics)


#' @title Plot of Compromise of Tables
#' @description Takes two dimensions of compromise of tables and returns a plot of the two dimensions
#' @param dimension1 vector of class \code{"numeric"}
#' @param dimension2 vector of class \code{"numeric"}
#' @param rownames_vec vector of row labels
#' @param img1 optional image of class \code{"png"}
#' @param img2 optional image of class \code{"png"}
#' @param img3 optional image of class \code{"png"}
#' @return plot of compromise of tables
#' @export
#' @examples
#'
#' dimension1 <- rnorm(20)
#' dimension2 <- rnorm(20)
#' plot_compromise(dimension1, dimension2)
#'
plot_compromise <- function(dimension1,
                            dimension2,
                            rownames_vec = as.character(1:length(dimension1)),
                            img1 = "black",
                            img2 = "black",
                            img3 = "black"){

  dat <- data.frame(x = dimension1, y = dimension2, label = rownames_vec)
  colnames(dat) <- c("x","y")

  ggplot(dat)+ geom_point(aes(x,y)) +labs(title="Compromise of the tables",
                                          x ="1", y = "2")+  ylim(-1, 1) + geom_text(data = dat, aes(x,y, label = rownames_vec), vjust = -2)+
    mapply(function(xx, yy)
      annotation_raster(img1, xmin=xx-0.08, xmax=xx+0.08, ymin=yy-0.04, ymax=yy+0.04),
      dat$x[1:4], dat$y[1:4]) +
    mapply(function(xx, yy)
      annotation_raster(img2, xmin=xx-0.08, xmax=xx+0.08, ymin=yy-0.04, ymax=yy+0.04),
      dat$x[5:8], dat$y[5:8]) +
    mapply(function(xx, yy)
      annotation_raster(img3, xmin=xx-0.08, xmax=xx+0.08, ymin=yy-0.04, ymax=yy+0.04),
      dat$x[9:12], dat$y[9:12]
    )

}

#' @title Plot of Partial Factor Scores
#' @description Takes two dimensions of partial factor scores and returns a plot of the two dimensions
#' @param dimension1 vector of class \code{"numeric"}
#' @param dimension2 vector of class \code{"numeric"}
#' @param rownames_vec vector of row labels
#' @param img1 optional image of class \code{"png"}
#' @param img2 optional image of class \code{"png"}
#' @param img3 optional image of class \code{"png"}
#' @return plot of partial factor scores
#' @export
#' @examples
#'
#' dimension1 <- rnorm(20)
#' dimension2 <- rnorm(20)
#' plot_compromise(dimension1, dimension2)
#'
plot_pfs <- function(dimension1,
                     dimension2,
                     rownames_vec = as.character(1:length(dimension1)),
                     img1 = "black",
                     img2 = "black",
                     img3 = "black"){

  dat <- data.frame(x = dimension1, y = dimension2, label = rownames_vec)

  ggplot(dat)+labs(title="Partial Factor Scores",
                               x ="1", y = "2") + geom_text(data = dat, aes(x,y, label = rownames_vec), vjust = -0.5)+ xlim(-1.5,1.8) + ylim(-1.3,1.3) +
    mapply(function(xx, yy)
      annotation_raster(img1, xmin=xx-0.08, xmax=xx+0.08, ymin=yy-0.04, ymax=yy+0.04),
      dat$x[1:4], dat$y[1:4]) +
    mapply(function(xx, yy)
      annotation_raster(img2, xmin=xx-0.08, xmax=xx+0.08, ymin=yy-0.04, ymax=yy+0.04),
      dat$x[5:8], dat$y[5:8]) +
    mapply(function(xx, yy)
      annotation_raster(img3, xmin=xx-0.08, xmax=xx+0.08, ymin=yy-0.04, ymax=yy+0.04),
      dat$x[9:12], dat$y[9:12])

}

#' @title Plot of Variable Loadings
#' @description Takes two dimensions of variable loadings and returns a plot of the two dimensions
#' @param dimension1 vector of class \code{"numeric"}
#' @param dimension2 vector of class \code{"numeric"}
#' @param rownames_vec vector of row labels
#' @return plot of variable loadings
#' @export
#' @examples
#'
#' dimension1 <- rnorm(20)
#' dimension2 <- rnorm(20)
#' plot_compromise(dimension1, dimension2)
#'
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

#' @title eigen_table
#' @description Takes the \code{"mfa"} object and returns a table with the singular values (i.e. square root of eigenvalues), the eigenvalues, cumulative, percentage of intertia, cumulative percentage of inertia, for all the extracted components
#' @param object an object of class \code{"mfa"}
#' @return object of class \code{"data.frame"}
#' @export
#' @examples
#'
#' simulate_data <- matrix(rnorm(200), ncol = 10, byrow = TRUE)
#' test <- mfa_const(simulate_data, sets = list(1:3, 4:5, 6:10), ncomps = 2)
#'
#' eigen_table(test)
#'
# Generic Methods
setGeneric(
  "eigen_table",
  function(object) standardGeneric("eigen_table")
)

# Set the method for S4 class
setMethod(
  "eigen_table",
  signature = "mfa",
  function(object){
    # Everything is based on eigenvalues
    eigenval = object@eigenvalues

    #  Create the data frame with a sample column
    samp = c(0,0,0,0,0)
    # rownames
    rownames = c("Singular value","Eigenvalue","cumulative E","% Inertia","cumulative I")
    df = data.frame(samp,row.names = rownames)
    # for loop for creating different arrays
    for (i in 1:(length(eigenval))){
      tmp_arr = c()

      tmp_sv = sqrt(eigenval[i])
      tmp_arr = append(tmp_arr,tmp_sv)

      tmp_ev = eigenval[i]
      tmp_arr = append(tmp_arr,tmp_ev)

      tmp_cumu = sum(df[2,1:i])+tmp_ev
      tmp_arr = append(tmp_arr,tmp_cumu)

      tmp_inert = (eigenval[i]/sum(eigenval))*100
      tmp_arr = append(tmp_arr,tmp_inert)

      tmp_cumu_inert = sum(df[4,1:i])+tmp_inert
      tmp_arr = append(tmp_arr,tmp_cumu_inert)

      assign(paste("Comp",i ,sep=""), tmp_arr)
      # produce the data frame
      df[,paste("Comp",i ,sep="")] = tmp_arr
    }
    # create the data frame, remove the sample column
    return(subset(df, select = -samp ))

  }

)


#' @title Contribution of an Observation to a Dimension
#' @description Outputs values that help interpret how the observations contribute to a dimension.
#' @param object an object of class \code{"mfa"}
#' @return object of class \code{"matrix"}
#' @export
#' @examples
#'
#' simulate_data <- matrix(rnorm(200), ncol = 10, byrow = TRUE)
#' test <- mfa_const(data = simulate_data, sets = list(1:3, 4:5, 6:10), ncomps = 2)
#' COD(test)
#'
# Generic Methods
setGeneric(
  "COD",
  function(object) standardGeneric("COD")
)
# Set the method for S4 class
setMethod(
  "COD",
  signature = "mfa",
  function(object){
    # Define the mass
    cfs = object@cfs
    n = length(cfs[,1])
    m = 1/n
    cfs_sq = cfs*cfs
    ctr = m*cfs_sq

    #scale by eigenvalues
    eigenval = object@eigenvalues
    for (i in 1:length(eigenval)){
      ctr[,i] = ctr[,i]/eigenval[i]
    }
    return(ctr)
  }
)


#' @title Contributions of a Variable to a Dimension
#' @description Outputs values that help interpret how the variables contribute to a dimension.
#' @param object an object of class \code{"mfa"}
#' @return an object of class \code{"list"}
#' @export
#' @examples
#'
#' simulate_data <- matrix(rnorm(200), ncol = 10, byrow = TRUE)
#' test <- mfa_const(data = simulate_data, sets = list(1:3, 4:5, 6:10), ncomps = 2)
#' CVD(test)
#'
# Generic Methods
setGeneric(
  "CVD",
  function(object) standardGeneric("CVD")
)
# Set the method for S4 class
setMethod(
  "CVD",
  signature = "mfa",
  function(object){
    cvd = list()
    alpha = object@alpha
    tmp = 0
    for (i in 1:length(object@mol)){
      mol = object@mol[[i]]
      mol_sq = mol*mol
      tmp = tmp+length(mol[,1])
      cvd_tmp = alpha[tmp] *mol_sq
      cvd[[i]] = cvd_tmp
    }
    return(cvd)
  }
)


#' @title Contribution of a Table to a Dimension
#' @description Outputs values that help interpret how the tables contribute to a dimension. The contribution of a table reﬂects the proportion of the variance of a dimension that can be attributed to this table.
#' @param object an object of class \code{"mfa"}
#' @return object of class \code{"list"}
#' @export
#' @examples
#'
#' simulate_data <- matrix(rnorm(200), ncol = 10, byrow = TRUE)
#' test <- mfa_const(data = simulate_data, sets = list(1:3, 4:5, 6:10), ncomps = 2)
#' CTD(test)
#'
# Generic Methods
setGeneric(
  "CTD",
  function(object) standardGeneric("CTD")
)
# Set the method for S4 class
setMethod(
  "CTD",
  signature = "mfa",
  function(object){
    # use CVD method to find Contributions of a Variable to a Dimension
    cvd = CVD(object)
    ctd = c()
    # sum columns to find the
    for (i in 1:length(cvd)){
      tmp_ctd = colSums(cvd[[i]])
      ctd = rbind(ctd,tmp_ctd)
    }
    # remove row names
    rownames(ctd) <- c()
    return(ctd)

  }
)

#private function to check if tabe is a matrix
check_table <- function(table){

  if(class(table) != "matrix"){
    stop("\n argument must be a of class 'matrix'")
  }
  TRUE
}

#private function to check if tabe is a matrix
check_rows <- function(table1, table2){

  if(nrow(table1) != nrow(table2)){

    stop("\n The two tables  must have the same number of rows.")

  }
  TRUE
}


#' @title RV Coefficient
#' @description Use the RV coefficient to evaluate the similarity between two tables
#' @param table1 an object of class \code{"matrix"} with same number of rows as table 2
#' @param table2 an object of class \code{"matrix"} with same number of rows as table 1
#' @return an object of class \code{"numeric"}
#' @export
#' @examples
#'
#' table1 <- matrix(rnorm(100), nrow = 10)
#' table2 <- matrix(rnorm(50), nrow = 10)
#' RV(table1, table2)
#'
RV<-function(table1, table2){

  check_table(table1)
  check_table(table2)
  check_rows(table1, table2)

  RV = sum(diag(table1%*%t(table1)%*%table2%*%t(table2)))/
    sqrt(sum(diag(table1%*%t(table1)%*%table1%*%t(table1)))*
           sum(diag(table2%*%t(table2)%*%table2%*%t(table2))))

  return(RV)

}

#' @title RV Coefficient Table
#' @description RV coefficient table contains information of how similar the individual tables are to each other pairwise.
#' @param dataset an object of class \code{"matrix"} or \code{"data.frame"}
#' @param sets object of class \code{"list"}
#' @return an object of class \code{"matrix"}
#' @export
#' @examples
#'
#' simulate_data <- matrix(rnorm(200), ncol = 10, byrow = TRUE)
#' RV_table(simulate_data, sets = list(1:3, 4:5, 6:10))
#'
RV_table <- function(dataset, sets){

  check_dataset(dataset)
  check_sets(sets)

  RV_mat <- matrix(NA,nrow = length(sets), ncol = length(sets))

  for(i in 1:length(sets)){
    for(j in 1:length(sets)){

      table1 = data.matrix(dataset[,sets[[i]]])
      table2 = data.matrix(dataset[,sets[[j]]])
      RV_mat[i,j] <- RV(table1, table2)

    }
  }
  return(RV_mat)
}

#' @title Lg Coefficient
#' @description Computes the Lg coefficient between two tables
#' @param table1 object of class \code{"matrix"} with same number of rows as table 2
#' @param table2 object of class \code{"matrix"} with same number of rows as table 1
#' @return object of class \code{"numeric"}
#' @export
#' @examples
#'
#' table1 <- matrix(rnorm(100), nrow = 10)
#' table2 <- matrix(rnorm(50), nrow = 10)
#' Lg(table1, table2)
#'
Lg<-function(table1, table2){

  check_table(table1)
  check_table(table2)
  check_rows(table1, table2)

  # Center the two tables first
  tmp_1 = scale(table1 , center = T,  scale = F)
  tmp_2 = scale(table2 , center = T,  scale = F)

  #Scale the data set (depend on the scale input)

  scale_array_1 = sqrt(apply(tmp_1^2,2,sum))
  scale_array_2 = sqrt(apply(tmp_2^2,2,sum))

  # Scale to nomalize matrix
  norm_1 = scale(tmp_1, center = F, scale = scale_array_1)
  norm_2 = scale(tmp_2, center = F, scale = scale_array_2)

  # alpha weight vector
  alpha1 = (1/(svd(norm_1)$d[1])^2)
  alpha2 = (1/(svd(norm_2)$d[1])^2)

  # Find the Lg coefficient
  Lg = sum(diag(table1%*%t(table1)%*%table2%*%t(table2)))*alpha1*alpha2
  return(Lg)
}

#' @title Lg Table
#' @description Computes a table of Lg coefficients.
#' @param dataset object of class \code{"matrix"} or \code{"data.frame"}
#' @param sets object of class \code{"list"}
#' @return object of class \code{"matrix"}
#' @export
#' @examples
#'
#' simulate_data <- matrix(rnorm(200), ncol = 10, byrow = TRUE)
#' Lg_table(simulate_data, sets = list(1:3, 4:5, 6:10))
#'
Lg_table <- function(dataset, sets){

  check_dataset(dataset)
  check_sets(sets)

  Lg_mat <- matrix(NA,nrow = length(sets), ncol = length(sets))
  for(i in 1:length(sets)){
    for(j in 1:length(sets)){
      table1 = data.matrix(dataset[,sets[[i]]])
      table2 = data.matrix(dataset[,sets[[j]]])
      Lg_mat[i,j] <- Lg(table1, table2)

    }
  }
  return(Lg_mat)
}




