# m6APrediction

m6APrediction is a BIO215 R package for predicting N6-methyladenosine (m6A) sites from nucleotide context and transcript features. It uses a fitted random forest to score individual candidate sites or batches of sites.

## Features

- Score one candidate site with `prediction_single()` or a feature table with `prediction_multiple()`.
- Convert each five-base context to the categorical predictors expected by the model.
- Use the bundled model and example CSV to run predictions without training a model.

## Installation

```r
install.packages("remotes")
remotes::install_github("Trajectory614/m6APrediction")
```

## Single-site prediction

```r
library(m6APrediction)

model <- readRDS(system.file("extdata", "rf_fit.rds",
                            package = "m6APrediction", mustWork = TRUE))
result <- prediction_single(
  model,
  gc_content = 0.6,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 12,
  distance_to_junction = 5,
  evolutionary_conservation = 0.8,
  DNA_5mer = "GGACA"
)
result
```

The bundled model returns a score of `"0.528"` and a status of `"Positive"` for this example.

## Batch prediction

```r
library(m6APrediction)

model <- readRDS(system.file("extdata", "rf_fit.rds",
                            package = "m6APrediction", mustWork = TRUE))
features <- read.csv(system.file("extdata", "m6A_input_example.csv",
                                 package = "m6APrediction", mustWork = TRUE))
prediction_multiple(model, head(features, 3), positive_threshold = 0.6)
```

## Input requirements

Read the CSV into a nonempty data frame with the following case-sensitive column names. Supply the same seven features as individual values for a single-site prediction. Required fields must have no missing values.

| Field | Accepted input |
| --- | --- |
| `gc_content` | Finite numeric value representing GC content. |
| `RNA_type` | `mRNA`, `lincRNA`, `lncRNA`, or `pseudogene`. |
| `RNA_region` | `CDS`, `intron`, `3'UTR`, or `5'UTR`. |
| `exon_length` | Finite numeric value on the log2 scale used in training. |
| `distance_to_junction` | Finite numeric value on the log2 scale used in training. |
| `evolutionary_conservation` | Finite numeric value representing evolutionary conservation. |
| `DNA_5mer` | Exactly five uppercase bases drawn from A, T, C, and G. The training contexts had A at the center. |

Supply length and distance values already transformed as in the training data. The prediction functions apply no further numerical transformation. RNA type and region must match the listed category labels exactly.

## Output interpretation

- `predicted_m6A_prob` is the fitted forest's vote fraction for the `Positive` class. It should be interpreted as a model score, not as an independently calibrated probability.
- `predicted_m6A_status` is `Positive` if the score is strictly greater than `positive_threshold`, and `Negative` otherwise. The threshold defaults to 0.5 and accepts a finite numeric value from 0 to 1 inclusive.

`prediction_multiple()` returns the input data frame with a numeric score column and a character status column. Existing columns with those output names are replaced. `prediction_single()` returns a two-element named character vector. Use `as.numeric(result[["predicted_m6A_prob"]])` to convert the single-site score to a number.

## Model information

The package includes `rf_fit.rds`, the 500-tree random forest fitted in Practical 5. The model was trained on all 48,055 rows of the course dataset and uses five position-specific nucleotide features, including `nt_pos3`. Load it with `readRDS(system.file(...))` as shown above and pass it to either prediction function.

## Historical ROC and precision–recall curves

The original figure below comes from Practical 3, Assignment Task 2. That random forest was fitted on 38,444 rows and evaluated with `precrec` on 9,611 held-out rows (2,273 positive and 7,338 negative). The Practical 3 assignment reported an AUROC of **0.8864** and an AUPRC of **0.7042**.

![ROC and precision–recall curves for the Practical 3 random forest on its held-out test set](figures/practical3-holdout-roc-prc.png)

The Practical 5 model bundled here was fitted on all 48,055 rows, including the 9,611 rows held out in Practical 3, and uses the additional feature `nt_pos3`. The displayed metrics therefore apply only to the Practical 3 model and provide no independent test of the bundled model.

## Limitations

The 100-row example CSV demonstrates the input format and prediction workflow. It was sampled from the training data. This repository provides no independent test results for the bundled model, and performance on new datasets remains unestablished here.
