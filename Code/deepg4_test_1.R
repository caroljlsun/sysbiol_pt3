library(keras)
library(tensorflow)
library(Biostrings)
library(DeepG4)

#what does tabv look like?
#tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1)
#tabv

#function dna_to_numerical
dna_to_numerical <- function(x,tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1),lower.case=F,seq.size = 1000){
  if(lower.case){ #tabv is a named vector of nucleotides and their assigned numbers
    names(tabv) <- tolower(tabv) #tolower translates upper to lower and vice versa, acts on character vectors. Here, it takes lower cases to upper cases
  }
  x <- Biostrings::as.matrix(x) #making the dna sequence a matrix 
  listMat <- list()                  #set up "listMat"
  for(i in 1:length(tabv)){          #for i in the length of tabv (which is 5), do the following:
    nuc_index <- tabv[[i]]           #assign to "nuc_index" the i'th element, simplified. So, it'll give 2 instead of C;2
    nuc_value <- names(tabv[i])      #assign to "nuc_value" the name of the i'th element. So, it'll give C 
    mat <- matrix(0,nrow(x),ncol(x)) #create an object mat, which is an x by x matrix full of zeros (000000000...)
    mat[x==nuc_value] <- 1           #whenever x is equal to the "nuc_value", e.g. is a C, assign 1 to the position instead of 0
    if(ncol(x)<seq.size){            #if the columns of x(a sequence matrix) is larger than your set "seq.size"...
      mat <- cbind(mat,matrix(0,nrow(x),seq.size-ncol(x))) #add more columns filled with 0s with the rows matching size of x
    }
    listMat[[nuc_index]] <- mat      #in the i'th subset, add the newly created matrix. e.g. in 2nd subset, add the C matrix
  }
  arrayout <- array(unlist(listMat), dim = c(nrow(listMat[[1]]), ncol(listMat[[1]]), length(listMat)))
  return(arrayout)                   #then create object arrayout, which lists out all the matrices you've made, with the appropriate dimensions
  
  
}


#examples
x <- Biostrings::DNAStringSet(c("ACGTNNTG"))
x_onehot <- dna_to_numerical(x)
x_onehot

# Next part of their code
DeepG4 <- function(X = NULL,Y=NULL,model=NULL,lower.case=F,treshold = 0.5,seq.size = 1000,retrain=FALSE,retrain.path="",log_odds=F){
  model.regular.size.accepted <- 1000 #I'm changing it all to 1000 instead of 201
  tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1)#as before
  #Check if X is provided
  if (is.null(X)) {
    stop("X must be provided (see ?DeepG4 for accepted formats).", #An object of class character,list or DNAStringSet/DNAStringSetList with DNA sequences.
         call. = FALSE)
  }
  # Packages check
  if (!requireNamespace("keras", quietly = TRUE)) {
    stop("Package \"keras\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package \"Biostrings\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # Check sequences and convert into one-hot
  ## Check model class and convert into DNAStringSet object if needed
  
  if(!class(X)[[1]] %in%c("DNAString","DNAStringSet","DNAStringSetList")){
    if(class(X) == "character"){
      X <- Biostrings::DNAStringSet(X)
    }else if(class(X) == "list"){
      if(class(X[[1]])[[1]] == "DNAString"){
        X <- as(X,"DNAStringSet")
      }else if(class(X[[1]])[[1]] == "character"){
        X <- Biostrings::DNAStringSet(unlist(X))
      }else{
        stop("X must be a list of DNAString/character class",
             call. = FALSE)
      }
    }else{
      stop("X must be a character, a list or a DNAStringSet/DNAStringSetList class",
           call. = FALSE)
    }
  }else if(class(X)[[1]] =="DNAStringSetList"){
    X <- unlist(Biostrings::DNAStringSetList(X))
  }else if(class(X)[[1]] =="DNAString"){
    X <- Biostrings::DNAStringSet(X)
  }
  ## Check sequences sizes
  message("Check sequences sizes...")
  seqsizes <- Biostrings::nchar(X)
  if(max(seqsizes) > seq.size){
    message(paste0("Warning: Some of your sequences are >",seq.size,", and will be croped."))
    X[(Biostrings::nchar(X)>seq.size & !Biostrings::nchar(X)<seq.size)] <- Biostrings::subseq(X[(Biostrings::nchar(X)>seq.size & !Biostrings::nchar(X)<seq.size)],start = 1,end = seq.size)
  }
  ## Check DNA composition
  message("Check sequences composition...")
  resFreq <- Biostrings::letterFrequency(X,"N",as.prob = T)
  testNFreq <- as.vector(resFreq>0.1)
  if(any(testNFreq)){
    message(paste0("Warning: Some of your sequences have a N frequency > 0.1 and will be removed.\nDeepG4 has difficulty to handle sequences with a N rate > 10%"))
    X <- X[!testNFreq]
    if(length(X)<1){
      stop("Not enough sequences to continue ...",
           call. = FALSE)
    }
  }
  ## One-Hot conversion
  message("One-Hot Conversion...")
  if(length(seqsizes) == 1) {
    X <- DNAToNumerical(X,tabv = tabv,lower.case=lower.case,seq.size = seq.size)
  }else{
    ## Have to apply One-Hot independently because seq sizes are differents
    X_by_size <- lapply(unique(Biostrings::nchar(X)),function(onesize){
      DNAToNumerical(X[Biostrings::nchar(X)==onesize],tabv = tabv,lower.case=lower.case,seq.size = seq.size)
    })
    X <- array(unlist(X_by_size), dim = c(length(X),seq.size,length(tabv)))
  }
  if(retrain){
    # IF RETRAIN = TRUE
    message("retrain == TRUE")
    message("Model will be retrain using user input...")
    # Check Y
    if(is.null(Y)){
      stop("Y must be set if you want retrain our model.",
           call. = FALSE)
    }
    if(class(Y) != "numeric"){
      stop("Y must be a numeric vector of 1 and 0 values",
           call. = FALSE)
    }
    if(FALSE %in% (unique(Y) %in% c(0,1))){
      stop("Y must be a numeric vector of 1 and 0 values",
           call. = FALSE)
    }
    # Build the model
    # Try to load our saved model or custom model if !is.null(model)
    message("Loading model...")
    if(is.null(model)){
      model <-  system.file("extdata", "model.hdf5", package = "DeepG4")
    }else{
      if(class(model) != "character"){
        stop("model must be a path to a keras model in hdf5 format",
             call. = FALSE)
      }
    }
    #Load model with keras (tensorflow must be installed as well)
    model <- keras::load_model_hdf5(model)
    model <- keras::from_config(keras::get_config(model))
    message("Compilation...")
    keras::compile(model,
                   optimizer = 'rmsprop',
                   loss = 'binary_crossentropy',
                   metrics = list('accuracy')
    )
    # Retrain the model
    message("Training...")
    history <- keras::fit(model,
                          X,
                          Y,
                          epochs = 20,
                          batch_size = 128,
                          validation_split = 0.2,
                          verbose= 0)
    message("Done !...")
    res <- stats::predict(model,X)
    if(retrain.path == ""){
      retrain.path <- paste0("DeepG4_retrained_",Sys.Date(),".hdf5")
    }
    keras::save_model_hdf5(model,retrain.path)
  }else{
    #IF RETRAIN = FALSE
    # Try to load our saved model or custom model if !is.null(model)
    message("Loading model...")
    if(is.null(model)){
      model <-  system.file("extdata", "model.hdf5", package = "DeepG4")
      if(seq.size != model.regular.size.accepted){
        message("Please don't manually set seq.size unless you want to use a custom model")
        seq.size <- model.regular.size.accepted
      }
    }else{
      if(class(model) != "character"){
        stop("model must be a path to a keras model in hdf5 format",
             call. = FALSE)
      }
    }
    #Load model with keras (tensorflow must be installed as well)
    model <- keras::load_model_hdf5(model)
    if(log_odds){
      # If log_odds is set to TRUE, return instead a real number computed by the layer before the sigmoid activation (or the last layer without sigmoid)
      model <- keras::keras_model(inputs = model$input,
                                  outputs = keras::get_layer(model, index = 7)$output)
    }
    res <- stats::predict(model,X)
  }
  # If Y is provided, instead of returning prediction, return accuracy / AUC
  if(is.null(Y)){
    return(res)
  }else{
    if(class(Y) != "numeric"){
      stop("Y must be a numeric vector of 1 and 0 values.",
           call. = FALSE)
    }
    if(FALSE %in% (unique(Y) %in% c(0,1))){
      stop("Y must be a numeric vector of 1 and 0 values.",
           call. = FALSE)
    }
    if(length(Y)!= nrow(res)){
      stop("Y must be a vector of same size as X.",
           call. = FALSE)
    }
    # Compute accuracy and AUC
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("yardstick", quietly = TRUE)) {
      stop("Package \"yardstick\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    prediction_table <- data.frame(
      truth = as.factor(Y),
      pred_prob = res[,1]
    )
    prediction_table$estimate <- factor(ifelse(prediction_table$pred_prob<treshold,0,1),levels = c(1,0))
    if(length(levels(prediction_table$truth))==1){
      message("DeepG4: metrics can't be evaluated with no control cases (length(levels(Y))==1), return predictions")
      return(res)
    }
    if(length(levels(prediction_table$estimate))==1){
      prediction_table$estimate <- factor(prediction_table$estimate,levels = c(1,0))
    }
    #Plot AUC
    prediction_table$truth <- factor(prediction_table$truth,levels = c(1,0))
    plot_ROC <- ggplot2::autoplot(yardstick::roc_curve(prediction_table,`truth`,`pred_prob`))
    #Get metrics
    table_metrics <- yardstick::metrics(prediction_table,`truth`,`estimate`,`pred_prob`)
    #Plot confusion matrix
    confusion_matrix <- as.data.frame(yardstick::conf_mat(prediction_table,`truth`, `estimate`)[[1]])
    
    confusion_matrix <- ggplot2::ggplot(confusion_matrix,ggplot2::aes(Prediction, `Truth`, fill = `Freq`)) +
      ggplot2::geom_tile(show.legend = FALSE) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::geom_text(ggplot2::aes(label = Freq), color = "white", alpha = 1, size = 8) +
      ggplot2::labs(
        title = "Confusion matrix"
      ) + ggplot2::theme_minimal(base_size=18)
    return(list(res,plot_ROC,confusion_matrix,table_metrics))
  }
  
  
}