
#' This script defines functions for calculating contrapositive consistency 
#' (c-consistency) and contrapositive coverage (c-coverage).
#' 
#' The functions are provided by Mathias Ambühl.

# load cna library
library(cna)


# ===== Description of the functions =====

#' The functions condTbl(), asf() and csf() will be most useful for cna 
#' users.
#' These functions take the same main arguments as the functions condTbl(), 
#' asf() and csf() from the regular cna library on CRAN.
#' They return the same results as their counterparts from the regular cna 
#' library, but with two additional columns containing contrapositive 
#' consistency (c-consistency) and contrapositive coverage (c-coverage).
#' 
#' Examples on how to use condTbl(), asf() and csf() as defined in this script 
#' are provided at the bottom of this script.


#' The functions contrapositives() and add_contrapositives() are helper 
#' functions: they will usually not be used directly by cna users but need to 
#' be defined because the functions condTbl(), asf() and csf() as defined in 
#' this script rely on them.
#' 
#' contrapositives() takes two arguments:
#' * cond : a character vector containing conditions. Only asfs or csfs in 
#' disjunctive normal form can be processed by contrapositives()
#' * d : a data.frame or configTable
#' The return value is a data.frame with two columns named "c_consistency" and 
#' "c_coverage".
#' 
#' add_contrapositives() takes two arguments, too:
#' * tbl : a condTbl, the output of functions contTbl(), asf(), csf()
#' * d : the data that was used to create tbl
#' It adds the columns c_consistency and c_coverage to the table and
#' rearranges the columns such that c_consistency and c_coverage are displayed 
#' before consistency and coverage.

#' The function fourfold will usually not be used directly by cna users, but 
#' may be helpful to some.
#' 
#' fourfold() takes a condition or asf and data (crisp-set or multi-value) 
#' arguments and prints a fourfold table (confusion matrix) for each condition.
#' 


# ===== Defining of the functions =====

# Contrapositives()
contrapositives <- function(cond, d){
  if (!inherits(d, "configTable")) d <- configTable(d)
  cti <- ctInfo(d)
  condTypes <- as.vector(cna:::.qcondType(cond, values = colnames(cti$scores), cti$type))
  is_atomic <- !is.na(condTypes) & condTypes == "stdAtomic"
  is_complex <- !is.na(condTypes) & condTypes == "stdComplex"
  ncond <- length(cond)
  out <- data.frame(c_consistency = rep(NA_real_, ncond), 
                    c_coverage = rep(NA_real_, ncond))
  # apply negation&invertion for asf components of complex conditions
  if (any(is_complex)){
    asfs <- extract_asf(cond[is_complex])
    # complex with only one asf are treated as asf:
    len2plus <- lengths(asfs) > 1
    unlistedAsfs <- unlist(asfs)
    negUnlisted <- paste0("!(", rhs(unlistedAsfs), ")->!(", lhs(unlistedAsfs), ")")
    negComplex <- split.default(negUnlisted, rep(seq_along(asfs), lengths(asfs)))
    negComplex <- paste0(ifelse(len2plus, "(", ""),
                         cna:::C_mconcat(negComplex, sep = ")*("), 
                         ifelse(len2plus, ")", ""))
    out[is_complex, ] <- cna::condTbl(negComplex, d)[c("consistency", "coverage")]
  }
  # negate both sides of atomic conditions and invert them
  if (any(is_atomic)){
    negAtomic <- paste0("!(", rhs(cond[is_atomic]), ")->!(", lhs(cond[is_atomic]), ")")
    out[is_atomic, ] <- cna::condTbl(negAtomic, d)[c("consistency", "coverage")]
  }
  out  
}

# add_contrapositives()
add_contrapositives <- function(tbl, d){
  if (!inherits(d, "configTable")) d <- configTable(d)
  cond <- tbl$condition
  tbl[c("c_consistency", "c_coverage")] <- contrapositives(cond, d)
  tbl[union(intersect(c("outcome", "condition", "c_consistency", "c_coverage"), 
                      names(tbl)), 
            names(tbl))]
}
# condTbl()
condTbl <- function(cond, d){
  co <- condition(cond, d)  
  std <- attr(co, "info")$condType %in% c("stdAtomic", "stdComplex")
  out <- as.condTbl(co)
  out$consistency[!std] <- out$coverage[!std] <- NA_real_
  add_contrapositives(out, d)
}
# asf()
asf <- function(x, ...){
  stopifnot(inherits(x, "cna"))
  add_contrapositives(cna::asf(x, ...), x$configTable)
}
# csf()
csf <- function(x, ...){
  stopifnot(inherits(x, "cna"))
  add_contrapositives(cna::csf(x, ...), x$configTable)
}

# fourfold()
fourfold <- function(cond, d){
  if (!inherits(d, "configTable")) d <- configTable(d)
  condTypes <- as.vector(getCondType(cond, ctInfo(d)))
  stopifnot(attr(d, "type") %in% c("cs", "mv"))
  if (!all(grepl("atomic", condTypes, ignore.case = TRUE))){
    stop("fourfold() accepts atomic conditions only.")
  }
  out <- lapply(condition(cond, d), .fourfold1)
  out <- array(unlist(out), dim = c(2, 2, length(out)), 
               dimnames = list(c("TRUE", "FALSE"), c("TRUE", "FALSE"), names(out)))
  out[is.na(out)] <- 0L
  out
}
# .fourfold1() is an auxiliary function used in fourfold
.fourfold1 <- function(cond1){
  nn <- attr(cond1, "n")
  cond1 <- lapply(cond1, as.logical)
  out <- tapply(nn, cond1, sum)[2:1, 2:1]
  out  
}


# ===== examples =====

#' Before running the examples below, make sure to execute all previous 
#' lines in this script to correctly define the functions
#' condTbl(), asf() and csf().


### crisp-set ###

#' Run a crisp-set CNA analysis.
ana <- cna(d.women, con = 0.7, cov = 0.7)

#' Apply the new asf() function to display the returned atomic solution 
#' formulas (i.e. single-outcome models) with additional columns for 
#' c-consistency and c-coverage.
asfs <- asf(ana)
asfs

#' Apply the new csf() function to display the returned complex solution 
#' formulas (i.e. multiple-outcome models) with additional columns for 
#' c-consistency and c-coverage.
csfs <- csf(ana)
csfs

#' Apply the new condTbl() function to display the c-consistency and 
#' c-coverage of the minimally sufficient conditions.
mscs <- msc(ana)
condTbl(mscs, d.women)


### multi-value ###

#' Run a multi-value CNA analysis.
ana <- cna(d.pban, con = 0.85, cov = 0.85)

#' Apply the new asf() function to display the returned atomic solution 
#' formulas (i.e. single-outcome models) with additional columns for 
#' c-consistency and c-coverage.
asfs <- asf(ana)
asfs

#' Apply the new csf() function to display the returned complex solution 
#' formulas (i.e. multiple-outcome models) with additional columns for 
#' c-consistency and c-coverage.
csfs <- csf(ana)
csfs

#' Apply the new condTbl() function to display the c-consistency and 
#' c-coverage of the minimally sufficient conditions.
mscs <- msc(ana)
print(condTbl(mscs, d.pban), 100)


### fuzzy-set ###

#' Run a fuzzy-set CNA analysis.
ana <- cna(d.pacts, con = 0.75, cov = 0.75)

#' Apply the new asf() function to display the returned atomic solution 
#' formulas (i.e. single-outcome models) with additional columns for 
#' c-consistency and c-coverage.
asfs <- asf(ana)
asfs

#' Apply the new csf() function to display the returned complex solution 
#' formulas (i.e. multiple-outcome models) with additional columns for 
#' c-consistency and c-coverage.
csfs <- csf(ana)
csfs

#' Apply the new condTbl() function to display the c-consistency and 
#' c-coverage of the minimally sufficient conditions.
mscs <- msc(ana)
condTbl(mscs, d.pacts)
