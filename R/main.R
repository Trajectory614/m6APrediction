#' Encode Five-Base DNA Contexts
#'
#' Convert each sequence to five categorical predictors used by the fitted model.
#' @param dna_strings Nonempty character vector of uppercase five-base DNA strings
#'   containing only A, T, C and G. Training contexts have A at the central base.
#' @return A data frame with five factor columns nt_pos1 through nt_pos5,
#'   each with levels A, T, C and G, and one row per sequence.
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "TGACC"))
dna_encoding <- function(dna_strings) {
  stopifnot(length(dna_strings) > 0L, !anyNA(dna_strings),
            all(grepl("^[ATCG]{5}$", dna_strings)))
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A for a Feature Table
#'
#' Apply an existing random forest without retraining. Positive calls require a
#' score strictly greater than the threshold; equality is Negative.
#' @param ml_fit A fitted randomForest classification model with a Positive class.
#' @param feature_df Nonempty data frame with gc_content, RNA_type, RNA_region,
#'   exon_length, distance_to_junction, evolutionary_conservation and DNA_5mer.
#'   The four continuous columns must be finite numeric values. Length and
#'   distance are already log2 transformed. RNA_type accepts mRNA, lincRNA,
#'   lncRNA or pseudogene; RNA_region accepts CDS, intron, 3'UTR or 5'UTR.
#'   DNA_5mer contains uppercase five-base DNA contexts (central A in training).
#' @param positive_threshold Finite numeric scalar between zero and one.
#' @return The original data frame with numeric predicted_m6A_prob and character
#'   predicted_m6A_status columns appended (existing columns with these names
#'   are replaced). Probabilities are forest vote fractions for Positive.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                   package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, feature_df[1:3, ], positive_threshold = 0.6)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  feature_names <- c("gc_content", "RNA_type", "RNA_region", "exon_length",
                     "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
  stopifnot(all(feature_names %in% colnames(feature_df)), nrow(feature_df) > 0L,
            length(positive_threshold) == 1L, is.finite(positive_threshold),
            positive_threshold >= 0, positive_threshold <= 1)
  # Keep the original table intact; only the model's working copy uses factors.
  model_data <- feature_df[, feature_names, drop = FALSE]
  model_data$RNA_type <- factor(model_data$RNA_type,
                               levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  model_data$RNA_region <- factor(model_data$RNA_region,
                                 levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  numeric_columns <- c("gc_content", "exon_length", "distance_to_junction",
                       "evolutionary_conservation")
  stopifnot(!anyNA(model_data),
            all(vapply(model_data[numeric_columns], is.numeric, logical(1))),
            all(is.finite(as.matrix(model_data[numeric_columns]))))
  model_data <- cbind(model_data, dna_encoding(as.character(model_data$DNA_5mer)))
  # Select Positive by name, never by an assumed column index.
  feature_df$predicted_m6A_prob <- predict(ml_fit, model_data, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > positive_threshold, "Positive", "Negative")
  return(feature_df)
}

#' Predict m6A for One Site
#'
#' Build a one-row feature table and call prediction_multiple().
#' @param ml_fit A fitted randomForest classification model with a Positive class.
#' @param gc_content Finite numeric scalar for GC content.
#' @param RNA_type One of mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region One of CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length Finite numeric scalar, already log2 transformed.
#' @param distance_to_junction Finite numeric scalar, already log2 transformed.
#' @param evolutionary_conservation Finite numeric scalar for conservation.
#' @param DNA_5mer Uppercase five-base A/T/C/G context; central A in training.
#' @param positive_threshold Finite numeric scalar in [0, 1]; Positive requires
#'   a score strictly greater than this value.
#' @return A named character vector of length two: predicted_m6A_prob (a numeric
#'   score represented as text) and predicted_m6A_status (Positive or Negative).
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'   RNA_region = "CDS", exon_length = 12, distance_to_junction = 5,
#'   evolutionary_conservation = 0.8, DNA_5mer = "GGACA")
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation,
                              DNA_5mer, positive_threshold = 0.5) {
  feature_df <- data.frame(gc_content, RNA_type, RNA_region, exon_length,
                           distance_to_junction, evolutionary_conservation, DNA_5mer)
  stopifnot(nrow(feature_df) == 1L)
  result <- prediction_multiple(ml_fit, feature_df, positive_threshold)
  # Mixed numeric/text atomic vectors are character vectors in R.
  returned_vector <- c(predicted_m6A_prob = as.character(result$predicted_m6A_prob[1]),
                       predicted_m6A_status = result$predicted_m6A_status[1])
  return(returned_vector)
}
