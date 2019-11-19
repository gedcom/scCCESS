###############################
# Written by Thomas A. Geddes #
# Last modified 2019/11/19    #
###############################

library("keras")
library("clue")
library("parallel")

encode = function(dat, seed = 1, max_random_projection = 2048, encoded_dim = 16, hidden_dims = c(128), learning_rate = 0.001, batch_size = 32, epochs = 100, verbose = 2, scale = FALSE, genes_as_rows = FALSE) {
    if (verbose[1] %in% 0:2) {
      verbose = verbose[1]
    } else {
      verbose = 1
    }
  
    set.seed(seed)
    was_data_frame = is.data.frame(dat)
    if (was_data_frame) dat = as.matrix(dat)
    if (class(dat) != "matrix") stop("Input data must be dataframe or matrix.")

    # Transpose
    if (genes_as_rows) dat = t(dat)

    # Strip row and column names
    datrows = rownames(dat)
    rownames(dat) = NULL
    colnames(dat) = NULL

    # Scale columns
    if (scale) {
        dat = apply(dat, 2, function(x) (x - mean(x)) / sd(x))
    }

    # Perform random projection
    num_input_features = ncol(dat)
    final_proj_dim = min(max_random_projection, ceiling(0.8 * num_input_features))
    random_proj_cols = sample(num_input_features, size = final_proj_dim, replace = FALSE)
    dat = dat[, random_proj_cols]

    # Clear deep learning graph
    keras::k_clear_session()
    
    # Construct encoder network
    tns = encoder_input = keras::layer_input(shape = final_proj_dim)
    for (layer_width in hidden_dims) {
        tns = keras::layer_dense(tns, units = layer_width)
        tns = keras::layer_activation_leaky_relu(tns, alpha = 0.01)
    }
    tns = keras::layer_dense(tns, units = encoded_dim)
    encoder = keras::keras_model(inputs = encoder_input, outputs = tns)
    
    # Construct decoder network
    tns = decoder_input = keras::layer_input(shape = encoded_dim)
    
    for (layer_width in rev(hidden_dims)) {
        tns = keras::layer_dense(tns, units = layer_width)
        tns = keras::layer_activation_leaky_relu(tns, alpha = 0.01)
    }
    
    tns = keras::layer_dense(tns, units = final_proj_dim)
    decoder = keras::keras_model(inputs = decoder_input, outputs = tns)

    # Combine networks
    tns = ae_input = keras::layer_input(final_proj_dim)
    tns = decoder(encoder(tns))
    autoencoder = keras::keras_model(inputs = ae_input, outputs = tns)
    keras::compile(autoencoder, optimizer = keras::optimizer_adam(lr = learning_rate), loss = 'mean_squared_error')
    
    # Fit autoencoder model
    keras::fit(autoencoder, dat, dat, batch_size = batch_size, epochs = epochs, verbose = verbose)
    
    # Encode input data, return rownames and return as original data type
    reduced_data = predict(encoder, dat, batch_size = batch_size)
    rownames(reduced_data) = datrows
    # if (genes_as_rows) reduced_data %<>% t
    if (was_data_frame) reduced_data = as.data.frame(reduced_data)

    return(reduced_data)
}


ensemble_cluster = function(dat, seed = 1, cluster_func = function(x) kmeans(x, centers=5), ensemble_sizes = c(1, 5, 10, 20, 50), cores = 1, ...) {
  
  individual_encodings = lapply(1:max(ensemble_sizes), function(i) {
    ae_seed = (seed - 1) * max(ensemble_sizes) + i
    cat("Seed ", ae_seed, "\tEnsemble ", i, "\n")
    encoding = encode(dat, seed = ae_seed, ...)
    cat("\n")
    return(encoding)
  })

  cat("Generating individual clusters...\n")
  individual_clusters = parallel::mclapply(individual_encodings, cluster_func, mc.cores = cores)

  cat("Generating consensus clusters...\n")  
  ensemble_clusters = parallel::mclapply(ensemble_sizes, function(size) {
    if (size == 1) {
      return(sample(individual_clusters, 1)[[1]]$cluster)
    } else {
      consensus = sample(individual_clusters, size, replace = FALSE)
      consensus = clue::cl_consensus(consensus, method = "HE")
      
      consensus = apply(as.matrix(consensus[[1]]), 1, function(row_vals) which(row_vals == 1))
      return(consensus)
    }
  }, mc.cores = cores)
  
  names(ensemble_clusters) = as.character(ensemble_sizes)
  
  return(ensemble_clusters)
}


