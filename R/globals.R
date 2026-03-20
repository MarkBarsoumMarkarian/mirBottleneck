# Suppress "no visible binding for global variable" NOTEs from dplyr pipelines
utils::globalVariables(c(
  # classify.R + composite_score.R
  "mirna", "vss", "coherence_score", "bottleneck_index", "class",
  "vss_norm", "coherence_norm", "coef", "hr", "p", "p_adj",
  "OS_days", "OS_status", "age", "stage_clean", "mirna_expr",
  "score_group", "patient", "bottleneck_score",
  # build_network.R (multiMiR result columns)
  "database", "mature_mirna_id", "target_symbol",
  "mirna_norm", "target", "n_targets"
))
