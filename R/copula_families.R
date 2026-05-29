# Input:  marginal  - families-object
# Output: character "continuous", "discrete" or "binary"
get_marginal_type <- function(marginal){
  fam_name <- tolower(attr(marginal, "name"))
  
  continuous <- c("gaussian", "gamma", "log-normal", "beta", "student t",
                  "log-log", "weibull")
  discrete <- c("zip", "negative binomial", "zinbi")
  binary <- c(0)
  
  if (fam_name %in% continuous) return("continuous")
  if (fam_name %in% discrete) return("discrete")
  if (fam_name %in% binary) return("binary")
  
  stop("Unknown marginal type for family: '", fam_name, "'")
}

# Input:  marginal1, marginal2   - families-object
# Output: character, e.g. "continuous_continuous", "binary_continuous",...
get_marginal_case <- function(marginal1, marginal2){
  type1 <- get_marginal_type(marginal1)
  type2 <- get_marginal_type(marginal2)
  
  paste0(type1, "_", type2)
}