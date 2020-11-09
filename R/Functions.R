#' @title Create a Factor from a vector
#' @description Use number of tables and observations, or a vector
#' containing the number of observations in every table to create a
#' Factor used un tables splitting for DS-PC analysis.
#'
#' @param n_values When \code{mode = "dim"}, \code{n_values} is a vector containing
#' the number of observations per table in \code{n_values[1]} and the number
#' of tables in \code{n_values[2]}. When \code{mode = "obs"}, n_values is a vector containing
#' the number of observations in each table.
#' @param mode It can be set up to "dimensions mode" when \code{mode = "dim"} or
#' "observations mode" when \code{mode = "obs"}(default).
#' @param name A character string used as prefix in the creation of the Factor.
#' It can be useful to distinguish different set of tables.
#'
#' @return A factor that can be used to splitting a matrix of values into a list of data tables.
#' @export
#'
#' @examples TabFactor( c(5,7) )
#' TabFactor( n_values = c(5,7), mode = "obs" )
#' TabFactor( n_values = c(5,7), mode = "dim" )
#' TabFactor( n_values = c(5,7,3,4) )
#' TabFactor( n_values = c(5,7,3,4), name = "W" )
TabFactor = function( n_values ,mode = "obs",name = "T" ){

  #validation
  if(mode!= "dim" && mode!="obs"){stop("'mode' argument is not valid")}
  if(!is.character(name)){stop("'name' argument must be of character class")}

  #"dim" mode uses equal observation numbers per table
  # n_values[1] is the number of observations per table and
  # n_values[2] is the number of tables

  if((mode=="dim") && (length(n_values)!=2)){stop("n_values must have two elements when 'dim' mode is chosen")}
  if( mode=="dim"){T_Factor = rep(x=1:n_values[2],each=n_values[1])}

  #"Obs" mode uses different number of observations per table
  # n_values is a vector holding the number of observations per table
  if( (mode =="obs") ){ T_Factor = rep(x = 1:length(n_values),times = n_values)}

  T_Factor=paste(rep(x = name,times=length(T_Factor)),T_Factor,sep="")

  return( factor(T_Factor) )
}




#' @title Create sequences over a factor
#' @description Use a factor to create a vector of
#' unique names for tables observations.
#'
#' @param FactorToName The factor from which names are created.
#'
#' @return A character vector of names, matching the length of \code{FactorToName}.
#' @export
#'
#' @examples
#' fac = TabFactor(  c(5,7,3,4) )
#' SeqFactor( fac )
#'
#' SeqFactor( c("A","A","A","B","B","B","C","C","C","C") )
SeqFactor = function( FactorToName ){
  NamesSequency = unsplit(lapply(split(x = FactorToName,f = FactorToName),seq_along),FactorToName)
  T_names = paste(FactorToName,rep(x=".",times=length(FactorToName)),NamesSequency,sep="")
  return(T_names)
}




#' @title  Create a Tables Object
#' @description Uses a matrix or data frame to organize the main parameters
#' needed to perform DS-PC analysis.
#'
#' @param JointTable A matrix or data frame containing the K reference tables
#' stacked downwards one after another. Is preferred that its columns and rows have names.
#' @param TableFactor A factor with the same number of entries as rows are in
#' \code{JointTable}.It is used to define the different tables involved in the analysis.
#' @param ObsNames A vector containing names for each row in \code{JointTable} (optional).
#' @param VarNames A vector containing names for each column in \code{JointTable} (optional).
#'
#' @return
#' A Table Object (a list) containing the separated tables
#' in \code{Xk_data}, the given table with names in \code{Original_data},
#' the \code{TableFactor} used to split the given table, the names
#' assigned to rows and columns in \code{VarNames} and \code{ObsNames},
#' the order of tables extracted from the factor, in \code{TablesOrder},
#' and the number of variables and observations in \code{Nvars} and \code{Nobs}
#' for computation purposes.
#'
#' @export
#'
#' @examples
#' data(iris)
#' TableObject(JointTable = iris[,1:4], TableFactor = iris[,5])
#' TableObject(iris[,1:4], iris[,5], ObsNames = SeqFactor(iris[,5]) )
#' TableObject(iris[,1:4], iris[,5], SeqFactor(iris[,5]), VarNames = paste("V",1:4,sep="") )
TableObject=function(JointTable,TableFactor,ObsNames=NULL,VarNames=NULL){
  Nvars=length(JointTable[1,])
  TableFactor=factor(TableFactor)

  #Validation
  if(length(JointTable[,1])!=length(TableFactor)){stop("number of rows of JointTable does not match with the number of elements in Table Factor")}

  if(is.null(VarNames)){ if(!is.null(colnames(JointTable))){VarNames=colnames(JointTable)} else{VarNames=paste(rep("V",times=Nvars),1:Nvars,sep="") }}
  if(is.null(ObsNames)){ if(!is.null(row.names(JointTable))){ObsNames=row.names(JointTable)} else{ObsNames= SeqFactor(TableFactor)}}

  rownames(JointTable) = ObsNames
  colnames(JointTable)=VarNames
  Original_data=as.matrix(JointTable)
  TablesOrder = unique(TableFactor)

  Xk_data=lapply(split(x = JointTable,f = TableFactor), as.matrix)
  Xk_data=Xk_data[TablesOrder]

  Nobs=sapply(Xk_data,function(M) length(M[,1]))

  Tobj=list(Xk_data=Xk_data,
           Original_data=Original_data,
           TableFactor=TableFactor,
           VarNames=VarNames,
           ObsNames=ObsNames,
           TablesOrder=TablesOrder,
           Nvars=Nvars,
           Nobs=Nobs)
  return(Tobj)}





#' @title Preprocessing of the Table Object
#' @description Performs an optional (but preferred) preprocessing over the \code{Xk_data}
#' in a Table Object.
#'
#' @param Tobj The table object containing the \code{Xk_data} to  be preprocessed
#' @param centering either a logical indicating whether center the data or not (usig global mean), a
#' value to perform the centering or a vector indicating a value for every variable.
#' @param scaling either a logical indicating whether scale the data or not (using global standard deviations), a
#' value to perform the centering or a vector indicating a value for every variable.
#' @param normalizing When TRUE, data is normalized dividing by the corresponding number
#' of rows in each table. If FALSE, the tables are not normalized. Optionally, a function which
#' takes each table as argument can be given, the value obtained will be divided from each
#' corresponding table.
#'
#' @return The Table object (\code{Tobj}) where \code{Xk_data} is preprocessed.
#' @export
#'
#' @examples
#' data(iris)
#' IrisTobj = TableObject(JointTable = iris[,1:4], TableFactor = iris[,5])
PreprocessTobj=function(Tobj,centering=TRUE,scaling=TRUE,normalizing=TRUE){
  Xbar = colMeans(Tobj$Original_data)
  Desv = apply(X=Tobj$Original_data,MARGIN = 2,FUN = sd)

  if( is.logical(centering) &&  isTRUE(centering)){ centering=Xbar }
  if( is.logical(scaling) &&  isTRUE(scaling)){ scaling=Desv }

  if( !is.logical(normalizing) && !is.function(normalizing) ){stop("normalizing must be either logical or a function")}
  if( is.logical(normalizing) &&  isTRUE(normalizing)){ normalizing = function(M) length(M[,1])}
  if( is.logical(normalizing) && !isTRUE(normalizing)){ normalizing = function(M) 1 }


  Tobj$Xk_data = lapply(Tobj$Xk_data,
                       function(M){
                         procM = scale(x=M,center = centering,scale = scaling)/normalizing(M)
                         attributes(procM)[["scaled:center"]]=NULL
                         attributes(procM)[["scaled:scale"]]=NULL
                         return(procM)
                       })
  return(Tobj)
}




#' @title Dual STATIS analysis for reference data
#' @description Performs Dual STATIS on a Table object derived from the \code{TableObject}
#' function. It provides all elements needed for analysis and monitoring
#'
#' @param Tobj The table object containing the \code{Xk_data} to  be preprocessed
#'
#' @return A list containing the Dual STATIS results, coded as \code{DSr}.
#' The projected coordinates are contained on the \code{PC_proj} list. The \code{G} matrix
#' has the projection of tables onto the Interstructure, the \code{V} matrix allows
#' new tables projection onto the Interstructure, the \code{F_scores} matrix
#' contains the principal components projections onto the Intrastructure, the \code{Fk} list
#' is composed by the partial factor scores matrix for every table, the \code{Pk} list
#' considers the projection matrices for every reference table and the \code{COj} list
#' contains the same partial factor scores as \code{Fk} organized by variable. The
#' \code{DSr} list, also has the \code{Xgsvd} list with the information of the generalized
#' singular value decomposition results. Additionally, original data from \code{Tobj} is
#' conserved
#' @export
#'
#' @examples
#' #Iris Data
#' data(iris)
#' IrisTobj = TableObject(JointTable = iris[,1:4], TableFactor = iris[,5])
#' IrisDSr = DualSTATIS( IrisTobj )
#' IrisDSr$PC_proj$G
#' IrisDSr$PC_proj$COj
#'
#' #DS-PC PlasticBags data
#' data( PlasticBags )
#' PBsTobj = TableObject( PlasticBags$Ref$data,  PlasticBags$Ref$factor)
#' PBsDSr = DualSTATIS( PreprocessTobj(PBsTobj)  )
#' PBsDSr$PC_proj$G
#' PBsDSr$PC_proj$COj
DualSTATIS=function(Tobj){

  #Fisrt Step
  Sk = lapply(X=Tobj$Xk_data,FUN = function(a) crossprod(a,a) )
  Z = t( sapply(X = Sk,FUN = as.vector))
  Z_pca = svd(Z)

  nPCS= length(Z_pca$d)
  PCSnames = paste("PC",1:nPCS,sep=" ")

  rownames(Z_pca$u)=Tobj$TablesOrder
  colnames(Z_pca$u)=PCSnames
  rownames(Z_pca$v)=paste("z", rep(1:Tobj$Nvars,times=Tobj$Nvars),rep(1:Tobj$Nvars,each=Tobj$Nvars),sep="")
  colnames(Z_pca$v)=PCSnames

  Z_pca$u=-Z_pca$u; Z_pca$v=-Z_pca$v

  #Alpha Wheights
  Aw = Z_pca$u[,1]/sum(Z_pca$u[,1])

  G = sweep(Z_pca$u,2,Z_pca$d,"*")
  V = Z_pca$v


  #Second Step
  X_stacked=t(matrix(unlist( lapply(Tobj$Xk_data,t)) ,nrow=Tobj$Nvars))
  rownames(X_stacked)=Tobj$ObsNames
  colnames(X_stacked)=Tobj$VarNames

  Aw_expanded = rep( x = Aw, times=Tobj$Nobs)

  X_tilde=sweep(X_stacked,1,sqrt(  Aw_expanded   ),"*" )
  X_pca = svd(X_tilde)

  #Delta
  Delta_vector =(X_pca$d)/(sqrt(Tobj$Nvars))
  CO_PCSnames = paste("PC",1:length(Delta_vector),sep=" ")
  names(Delta_vector)=CO_PCSnames

  #P
  P = sweep(X_pca$u,1,sqrt(  Aw_expanded   ),"/" )
  rownames(P)=Tobj$ObsNames
  colnames(P)=CO_PCSnames

  #Q
  Q = (X_pca$v)*(sqrt(Tobj$Nvars))
  rownames(Q)=Tobj$VarNames
  colnames(Q)=CO_PCSnames

  F_scores = sweep(Q,2,Delta_vector,"*")
  rownames(F_scores)=Tobj$VarNames
  colnames(F_scores)=CO_PCSnames

  #Fk obtaining
  colnames(P)= CO_PCSnames
  Pk = (lapply(split(as.data.frame(P),Tobj$TableFactor ),as.matrix) )[Tobj$TablesOrder]
  Fk = mapply(function(M1,M2){t(M1)%*%M2},Tobj$Xk_data, Pk ,SIMPLIFY = FALSE)

  COj = lapply(1:Tobj$Nvars,function(k){
    TempMat = t( sapply(1:length(Tobj$TablesOrder),function(s) (Fk[[s]])[k,] ) )
    rownames(TempMat) = Tobj$TablesOrder
    return(TempMat)
  }  )

  #Sk_projection_matrix
  Sk_proj= sweep(Q,2,Delta_vector,"/")/Tobj$Nvars

  #Principal components and projection matrices
  PC_proj =  list(G=G,
                  V=V,
                  F_scores=F_scores,
                  Fk=Fk,
                  Pk=Pk,
                  COj=COj)

  Xgsvd=list(P=P,Delta_vector=Delta_vector,Q=Q)

  #Dual Statis Results
  DSr = list(PC_proj=PC_proj,
             Sk_proj=Sk_proj,
             Xgsvd=Xgsvd,
             Xk_data=Tobj$Xk_data,
             Original_data=Tobj$Original_data,
             TableFactor=Tobj$TableFactor,
             VarNames=Tobj$VarNames,
             ObsNames=Tobj$ObsNames,
             TablesOrder=Tobj$TablesOrder,
             Nvars=Tobj$Nvars,
             Nobs=Tobj$Nobs)

  return(DSr)
}




#' @title Dual STATIS projection of new tables
#' @description Projects new tables from a Table object derived from the \code{TableObject}
#' function. It provides all elements needed for monitoring of new tables.
#'
#' @param Dsrs The Dual STATIS results object containing the processsed reference information
#' @param newTobj The table object containing the prepared information of new tables
#'
#' @return A list containing the projection of new tables, coded as \code{nDSr}.
#' The projected coordinates are contained on the \code{PC_proj} list.
#' Additionally, original data from \code{newTobj} is conserved
#' @export
#'
#' @examples
#' #DS-PC PlasticBags data
#' data( PlasticBags )
#' PBsTobj = TableObject( PlasticBags$Ref$ data,  PlasticBags$Ref$factor )
#' Xbar =  colMeans(PBsTobj$Original_data); Desv = apply(PBsTobj$Original_data,2,sd)
#' PBsDSr = DualSTATIS( PreprocessTobj(PBsTobj)  )
#' newTobj = TableObject( PlasticBags$Add$data,  PlasticBags$Add$factor )
#' procTobj = PreprocessTobj(newTobj, centering = Xbar, scaling = Desv)
#' NewTabProj = DUAL_STATIS_projection( PBsDSr, procTobj )
#' NewTabProj$PC_proj$G
DualSTATISprojection = function(Dsrs,newTobj){
  #Validate number of variables?
  if(Dsrs$Nvars!=newTobj$Nvars){stop("Number of variables between results and new table object must match") }

  newSk = lapply(newTobj$Xk_data,function(M) crossprod(M) )
  newZ = t( sapply(newSk,as.vector))
  newG = newZ%*%Dsrs$PC_proj$V

  newPk = lapply(newTobj$Xk_data,function(M) M%*%Dsrs$Sk_proj )
  newFk = lapply(newSk,function(M)  M%*%Dsrs$Sk_proj )

  newCOj = lapply(1:newTobj$Nvars,function(k){
    TempMat = t( sapply(1:length(newTobj$TablesOrder),function(s) (newFk[[s]])[k,] ) )
    rownames(TempMat) = newTobj$TablesOrder
    return(TempMat)
  }  )

  #Principal components and projection matrices
  PC_proj =  list(G=newG,
                  Fk=newFk,
                  Pk=newPk,
                  COj=newCOj)

  nDsr = list(PC_proj=PC_proj,
              Xk_data=newTobj$Xk_data,
              Original_data=newTobj$Original_data,
              TableFactor=newTobj$TableFactor,
              VarNames=newTobj$VarNames,
              ObsNames=newTobj$ObsNames,
              TablesOrder=newTobj$TablesOrder,
              Nvars=newTobj$Nvars,
              Nobs=newTobj$Nobs)
  return(nDsr)
}





#' @title Computing of bagplots for projected tables
#' @description Computing of hull and center of bagplots for Interstructure and
#' variables in the Intrastructure contained in a Dual STATIS results list derived
#' from \code{DualSTATIS} function.
#'
#' @import aplpack
#' @param Dsrs A list containing Dual STATIS results obtained from \code{DualSTATIS} function,
#'
#' @return A list containing coordinates of points conforming the hull and center of
#' the bagplots ordered counterclockwise.
#'
#' @export
#' @examples
#' #DS-PC PlasticBags data
#' data( PlasticBags )
#' PBsTobj = TableObject( PlasticBags$Ref$data,  PlasticBags$Ref$factor)
#' PBsDSr = DualSTATIS( PreprocessTobj(PBsTobj)  )
#' GenBagplots( PBsDSr )
GenBagplots = function(Dsrs){

  modified_bp = function(DataPoints){
    tempBP = compute.bagplot(DataPoints)
    if(is.null(tempBP$hull.loop) ){tempBP$hull.loop=DataPoints[chull(DataPoints),]}
    colnames(tempBP$hull.loop) = c("PC1","PC2")
    rownames(tempBP$hull.loop)= rep("",length(tempBP$hull.loop[,1]))
    names(tempBP$center) = c("PC1","PC2")

    modifiedbagplot = list(loop=tempBP$hull.loop, center=tempBP$center)
    return(modifiedbagplot)
  }

  IS_bp=modified_bp(Dsrs$PC_proj$G[,c(1,2)])
  if( is.null(compute.bagplot(Dsrs$PC_proj$G[,c(1,2)])$hull.loop)){
    warning( "There are few data points, a simple Convex Hull was computed" )}
  CO_bp = lapply(Dsrs$PC_proj$COj, function(M) modified_bp(M[,c(1,2)]) )

  DS_Bags = list(IS_bp=IS_bp, CO_bp = CO_bp)
  return(DS_Bags)
}





#' @title Parallel Coordinates Plot
#' @description A modified version of the Parallel Coordinates plotting function
#' available in the MASS package. If number of rows in the table
#' is greater than 500, a sample is picked.
#'
#' @param x  A multivariate data matrix representing reference batches or
#' a testing batch.
#' @param col A vector of colours, recycled as necessary for each observation.
#' @param lty A vector of line types, recycled as necessary for each observation.
#' @param var.label If TRUE, each variable's axis is labelled with maximum and minimum values.
#' @param liminf a numeric vector conformed by the inferior limit for every variable.
#' @param limsup a numeric vector conformed by the inferior limit for every variable.
#' @param ... Further graphics parameters which are passed to \code{matplot}
#'
#' @details This function creates a parallel coordinates plot using the \code{matplot}
#' function. No assignable output is created.
#' @export
#' @author B. D. Ripley. Enhancements based on ideas and code by Fabian Scheipl.
#' @references Wegman, E. J. (1990) Hyperdimensional data analysis using parallel
#' coordinates. Journal of the American Statistical Association 85, 664 - 675.
#'
#' Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S.
#' Fourth edition. Springer.
#' @examples
#' data( PlasticBags )
#' x = PlasticBags$Add$data[1:100,]
#' colorset = c(rep("black",50),rep("blue",50))
#' Xbar = colMeans(PlasticBags$Ref$data); Desv = apply(PlasticBags$Ref$data,2,sd)
#' inferiorlim = Xbar - 5*Desv
#' superiorlim = Xbar + 5*Desv
#' parcoord2( x, var.label = TRUE , col= colorset ,liminf = inferiorlim, limsup = superiorlim  )
parcoord2=function(x,col = 1,lty = 1,var.label = FALSE,liminf = NULL,limsup=NULL, ...){
  if(is.null(liminf)) liminf=apply(x, 2, min, na.rm = TRUE)
  if(is.null(limsup)) limsup=apply(x, 2, max, na.rm = TRUE)
  if( nrow(x)>500){x = x[sample(1:nrow(x),500),]}
  #ranges
  rx <- rbind(liminf,limsup)

  # Range-Standarization (0 to 1) of axes
  colnames(rx)=colnames(x)
  x <- sapply(1:ncol(x),function(j) (x[,j]-liminf[j])/(limsup[j]-liminf[j])  )
  colnames(x)=colnames(rx)

  matplot(1L:ncol(x), t(x), type = "l", col = col, lty = lty,
          xlab = "", ylab = "", axes = FALSE, ylim=c(0,1), ...)
  axis(1, at = 1L:ncol(x), labels = colnames(x))


  for (i in 1L:ncol(x)) {
    lines(c(i, i), c(0, 1), col = "grey55")
    if (var.label)
      text(c(i, i), c(0, 1), labels = format(rx[, i], digits = 3),
           xpd = NA, offset = 0.3, pos = c(1, 3), cex = 0.7)
  }
  invisible()
}




#' @title Automatic processing of data using Dual STATIS
#' @description Default routine for automatic processing of reference batches and/or
#' testing batches using Dual STATIS for the DS-PC approach
#'
#' @param RefTable A matrix or data frame containing the K reference tables stacked
#' downwards one after another. Is preferred that its columns and rows have names.
#' @param AddTable A matrix or data frame containing the tables of testing batches
#' stacked downwards one after another. Is preferred that its columns and rows have names.
#' @param FactOnTable if \code{TRUE}, \code{RefFactor} and/or \code{AddFactor} are defined
#' as the first column  of \code{RefTable} or \code{AddTable} respectively.
#' @param RefFactor A factor with the same number of entries as rows are in \code{RefTable}.
#' It is used to define the different tables involved in the analysis.
#' @param AddFactor A factor with the same number of entries as rows are in \code{AddTable}.
#' It is used to define the different tables involved in the analysis.
#' @param centering either a logical indicating whether center the data or not (using global mean),
#'  a value to perform the centering or a vector indicating a value for every variable.
#' @param scaling either a logical indicating whether scale the data or not (using global standard
#' deviations), a value to perform the centering or a vector indicating a value for every variable.
#' @param normalizing When TRUE, data is normalized dividing by the corresponding number of rows
#' in each table. If FALSE, the tables are not normalized. Optionally, a function which takes
#' each table as argument can be given, the value obtained will be divided from each corresponding table.
#'
#' @return
#' A list containing the Dual STATIS results of the processing. This includes the bagplots
#' of Interstructure and Intrastructure for monitoring.
#'
#' @export
#'
#' @examples
#' #DS-PC PlasticBags data
#' data( PlasticBags )
#' REFtab = PlasticBags$Ref$data
#' REFfac = PlasticBags$Ref$factor
#' process.results = AutoProcessing(RefTable = REFtab,FactOnTable = FALSE,RefFactor = REFfac )
#' process.results$REF$PC_proj$G
#' process.results$BPS$IS_bp
AutoProcessing = function(RefTable,RefFactor,AddTable=NULL,AddFactor =NULL,
                          centering = TRUE, scaling = TRUE, normalizing =TRUE){

  #Internal validation

  if( is.null(RefTable) || is.null(RefFactor) ){stop("RefTable and RefFactor must be given")}else{
      if(length(RefFactor)!=length(RefTable[,1])){stop("Number of rows of RefFactor must match the number of elements in RefTable")}
    }
  if( !is.null(AddTable) ){
      if( is.null(AddFactor) ){stop("When AddTable is not NULL, AddFactor must be given")}else{
        if(length(AddFactor)!=length(AddTable[,1])){stop("Number of rows of AddFactor must match the number of elements in AddTable")}
      }
    }

  #(scaling options are validated by the scaling function)

  #Processing of Reference
  RefTableObject = TableObject(RefTable,RefFactor,rownames(RefTable),colnames(RefTable))
  PrepRefTableObject = PreprocessTobj(RefTableObject,centering,scaling,normalizing)
  RefDSresults = DualSTATIS(PrepRefTableObject)
  RefBagplots = GenBagplots(RefDSresults)

  #Processing of Additional Tables, if given
  if( !is.null(AddTable) ){
    AddTableObject = TableObject(AddTable,AddFactor,rownames(AddTable),colnames(AddTable))
    PrepAddTableObject = PreprocessTobj(AddTableObject,centering,scaling,normalizing)
    AddDSresults = DualSTATISprojection(RefDSresults,PrepAddTableObject)
  }else{
    AddDSresults=NULL
  }

  return( list(REF = RefDSresults, ADD = AddDSresults, BPS = RefBagplots) )
}





#' @title Automatic plotting of Dual STATIS results
#' @description Default routine for automatic plotting of reference batches and
#' testing batches as part of the DS-PC approach
#'
#' @param REF A Dual STATIS results object from reference tables.
#' @param ADD A Dual STATIS results object from additional tables.
#' @param BPS A list of bagplot hulls and centers for Interstructure and Intrastructure.
#' @param plotIS A logical value indicating whether the Interstructure must be plotted or not.
#' @param plotCO A logical value indicating whether the Intrastructure must be plotted or not.
#' @param plotPC A logical value indicating whether the Parallel Coordinates must be plotted or not.
#' If number of observations in Reftable or Addtable is greater than 250, a sample is picked.
#' @param plotBags A logical value indicating whether the bagplots must be plotted or not.
#'
#' @details This function creates graphs for Interstructure and Intrastructure
#' analysis of reference and/or additional batches
#'
#' @export
#'
#' @examples
#' #DS-PC PlasticBags data
#' data( PlasticBags )
#' REFtab = PlasticBags$Ref$data[1:1500,]
#' ADDtab = PlasticBags$Add$data
#' REFfac = PlasticBags$Ref$factor[1:1500]
#' ADDfac = PlasticBags$Add$factor
#' DSPCres = AutoProcessing(REFtab, REFfac, ADDtab, ADDfac)
#' AutoPlotting( DSPCres$REF, DSPCres$ADD, DSPCres$BPS )
AutoPlotting = function(REF, ADD = NULL, BPS = NULL, plotIS = TRUE,
                        plotCO = TRUE, plotPC = TRUE, plotBags= TRUE){

  #Validations
  if( is.null(REF) ){stop("REF must be given")}
  if(plotBags && is.null(BPS)){stop("When plotBags is TRUE, BPS must be given")}
  if(!is.logical(plotIS) ){stop("plotIS must be a logical value")}
  if(!is.logical(plotCO) ){stop("plotCO must be a logical value")}
  if(!is.logical(plotPC) ){stop("plotPC must be a logical value")}
  if(!is.logical(plotBags) ){stop("plotBags must be a logical value")}
  if( isTRUE(plotPC)){
    if(nrow(REF$Original_data) > 250){
      REF$Original_data = REF$Original_data[sample(1:nrow(REF$Original_data),250),]  }
  }

  if( is.null(ADD)){
    #STATIS analysis mode
    if(plotIS){
      plot(REF$PC_proj$G[,1],REF$PC_proj$G[,2],pch=19,col="black",cex=0.7,
           main = "INTERSTRUCTURE",xlab = "PC1",ylab = "PC2")
      text(REF$PC_proj$G[,1],REF$PC_proj$G[,2],labels = as.character(REF$TablesOrder),pos = 4,cex=0.7)
      if(plotBags){
        polygon(BPS$IS_bp$loop)
        points( t(BPS$IS_bp$center),pch=8,col="red",cex=0.7 )
      }
    }
    if(plotCO){
      invisible(lapply(1:REF$Nvars,function(j){
        plot( (REF$PC_proj$COj[[j]])[,1],(REF$PC_proj$COj[[j]])[,2],pch=19,col="black",cex=0.7,
              main = bquote(paste('CO'[italic(j)]," Control Chart - ", .(REF$VarNames[j]) )) ,
              xlab = "PC1",ylab = "PC2")
        text( (REF$PC_proj$COj[[j]])[,1],(REF$PC_proj$COj[[j]])[,2],
              labels = as.character(REF$TablesOrder),pos=4,cex=0.7)
        if(plotBags){
          polygon((BPS$CO_bp[[j]])$loop)
          points( ((BPS$CO_bp[[j]])$center)[1] ,((BPS$CO_bp[[j]])$center)[2],pch=8,col="red",cex=0.7 )
        }
      }))
    }
    if(plotPC){
      #colors:  black = 1, blue = 4, red = 2
      plotPCtab = rbind( REF$Original_data, colMeans(REF$Original_data) )
      plotPCcol = c( rep(1,nrow(REF$Original_data) ), 2)
      parcoord2(plotPCtab,col = plotPCcol,var.label = TRUE, main = "Parallel Coordinates" )
    }
  }else{

    #STATIS monitoring mode
    if(plotIS){
      plot( rbind(REF$PC_proj$G[,1:2],ADD$PC_proj$G[,1:2]),pch =19,cex=0.7,
            main = "INTERSTRUCTURE",xlab = "PC1",ylab = "PC2",
            col= c(rep(1,length(REF$TablesOrder)),rep(4,length(ADD$TablesOrder)))   )
      text(c(REF$PC_proj$G[,1],ADD$PC_proj$G[,1]),c(REF$PC_proj$G[,2],ADD$PC_proj$G[,2]),
           labels = c(as.character(REF$TablesOrder),as.character(ADD$TablesOrder)),pos = 4,cex=0.7)

      if(plotBags){
        polygon(BPS$IS_bp$loop)
        points( t(BPS$IS_bp$center),pch=8,col="red",cex=0.7 )
      }
    }
    if(plotCO){
      invisible(lapply(1:REF$Nvars,function(j){
        plot( rbind( (REF$PC_proj$COj[[j]])[,1:2],(ADD$PC_proj$COj[[j]])[,1:2] ),pch =19,cex=0.7,
              main = bquote(paste('CO'[italic(j)]," Control Chart - ", .(REF$VarNames[j]) )) ,
              xlab = "PC1",ylab = "PC2",
              col= c(rep(1,length(REF$TablesOrder)),rep(4,length(ADD$TablesOrder)))   )
        text(c( (REF$PC_proj$COj[[j]])[,1],(ADD$PC_proj$COj[[j]])[,1] ),
             c( (REF$PC_proj$COj[[j]])[,2],(ADD$PC_proj$COj[[j]])[,2] ),
             labels = c(as.character(REF$TablesOrder),as.character(ADD$TablesOrder)),pos = 4,cex=0.7)

        if(plotBags){
          polygon((BPS$CO_bp[[j]])$loop)
          points( t((BPS$CO_bp[[j]])$center),pch=8,col="red",cex=0.7 )
        }
      }))
    }
    if(plotPC){
      Add_indices_list = (split(seq_along(ADD$TableFactor),ADD$TableFactor))[ADD$TablesOrder]

      #colors:  black = 1, blue = 4, red = 2
      invisible(lapply( 1:length(ADD$TablesOrder), function(Tab){
        plotPCind = Add_indices_list[[Tab]]
        if( length(plotPCind)>250 ){plotPCind = sample(plotPCind,250)  }
        plotPCtab = rbind( REF$Original_data, ADD$Original_data[ plotPCind, ] ,colMeans(REF$Original_data) )
        plotPCcol = c( rep(1,nrow(REF$Original_data) ), rep(4,length( plotPCind ) ), 2)
        parcoord2(plotPCtab,col = plotPCcol,var.label = TRUE,main = paste("Parallel Coordinates - ",ADD$TablesOrder[Tab] )  )
      }))
    }

  }
}

