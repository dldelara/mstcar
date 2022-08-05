#' Heart Attack Mortality Data in North Carolina
#'
#' A dataset containing counts for Myocardial Infarction (ICD-9 code 410.0)
#' deaths in North Carolina counties for individuals in 4 age groups from 1979-1988. This
#' dataset also contains the corresponding population counts and county FIPS codes.
#' The variables are as follows:
#'
#' @format a data frame with 100 rows and 6 variables:
#'
#' \describe{
#'   \item{year}{year of death count (1979)}
#'   \item{agegrp}{age group of individuals counted (65-74)}
#'   \item{county}{name of county for event count}
#'   \item{FIPS}{FIPS code for corresponding county name}
#'   \item{deaths}{death count (1-129)}
#'   \item{population}{population count (400-21061)}
#' }
#' @source \url{https://wonder.cdc.gov/cmf-icd9.html}
"ncheart"

#' North Carolina Neighbor Data
#'
#' A dataset containing U.S. Census TIGER shape data for North Carolina
#'
#' @format a list with two objects:
#'
#' \describe{
#'   \item{neigh}{List of neighbor adjacency vectors}
#'   \item{num}{Vector of number of neighbors for each county}
#' }
"ncnb"

#' North Carolina Shapefile
#'
#' A dataset containing U.S. Census TIGER shape data for North Carolina
#'
#' @format a data frame with 100 rows and 6 variables:
#'
#' \describe{
#'   \item{STATEFP}{State FIPS code}
#'   \item{COUNTYFP}{County FIPS code}
#'   \item{COUNTYNS}{County GNIS code}
#'   \item{AFFGEOID}{Census Unique Identifier}
#'   \item{GEOID}{Census Unique Identifier, truncated}
#'   \item{NAME}{County Name}
#'   \item{LSAD}{County LSAD code}
#'   \item{ALAND}{amount of land in square meters}
#'   \item{AWATER}{amount of water in square meters}
#'   \item{geometry}{shape data for each county}
#' }
#' @source \url{https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html}
"ncshp"
