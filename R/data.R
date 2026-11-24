#' Energy Efficiency data
#'
#' Heating and cooling loads of 768 simulated residential building configurations, from the
#' UCI Machine Learning Repository. Surface area equals wall area plus twice the roof area,
#' so the paper's benchmark drops \code{roof_area} and predicts \code{heating_load} from the
#' remaining seven inputs.
#'
#' @format A data frame with 768 rows and 10 columns: \code{relative_compactness},
#'   \code{surface_area}, \code{wall_area}, \code{roof_area}, \code{overall_height},
#'   \code{orientation}, \code{glazing_area}, \code{glazing_area_distribution},
#'   \code{heating_load} and \code{cooling_load}.
#' @source Tsanas, A. and Xifara, A. (2012). Energy Efficiency [Dataset]. UCI Machine Learning
#'   Repository, \doi{10.24432/C51307}. Licensed under CC BY 4.0.
#' @references Tsanas, A. and Xifara, A. (2012). Accurate quantitative estimation of energy
#'   performance of residential buildings using statistical machine learning tools. Energy and
#'   Buildings, 49, 560--567.
"energy"
