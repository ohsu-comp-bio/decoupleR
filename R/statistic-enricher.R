#' Enricher Python module used to predict the activity of regulatory proteins
#'
#' @description
#' Calculates regulatory activities using Enricher
#'
#' @details
#' This function is a wrapper for the python method [enricher.enrich()].
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param cohort Character, cohort to be processed.
#' @param regulon_size Integer indicating the minimum number of targets allowed per
#' regulon.
#' @param thresh.filter Float, Prior to normalization remove features that have a standard deviation per feature less
#' than {thresh_filter}
#' @param scaler_type Character indicating whether to scale and by what method to scale dataset
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#'
#' @import dplyr
#' @import tibble
#' @import purrr
#' @import tidyr
#' @import reticulate
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_enrich(mat, network, .source='tf', verbose = FALSE)

run_enrich <- function(mat,
                      network,
                      .source = .data$source,
                      .target = .data$target,
                      .mor = .data$mor,
                      .likelihood = .data$likelihood,
                      minsize = 5,
                      scaler_type = NULL,
                      ...) {
    # Check for NAs/Infs in mat
    check_nas_infs(mat)
    use_condaenv('enrich_env')
    enr <- import('enricher.enrich')

    # Before to start ---------------------------------------------------------
    network <- network %>%
        rename_net({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})
    network <- filt_minsize(rownames(mat), network, minsize)
    network <- network %>% convert_f_defaults( UpGene = source, DownGene = target) %>% dplyr::mutate(Type = 'controls-expresssion-of') %>% select(UpGene,Type, DownGene) 
    pd_network <- enr$pd$DataFrame(network)
    pd_mat <- enr$pd$DataFrame(mat, index=rownames(mat), columns=colnames(mat))

    # Analysis ----------------------------------------------------------------
    enr_scores <- enr$Enrichment(cohort='decoupler',expr=pd_mat,regulon=pd_network,regulon_size=minsize)
    enr_scores$scale(scaler_type=scaler_type)
    enr_scores$assign_weights()
    enr_scores$calculate_enrichment()
    t(enr_scores$total_enrichment) %>%
    as.data.frame() %>%
    rownames_to_column("source") %>%
    pivot_longer(-.data$source, names_to = "condition", values_to = "score") %>%
    add_column(statistic = "enricher", .before = 1) %>%
    mutate(p_value = 2*stats::pnorm(-abs(.data$score)))
}
