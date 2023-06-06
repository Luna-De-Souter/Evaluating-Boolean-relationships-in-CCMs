

# This script defines functions for calculating contrapositive consistency
# and contrapositive coverage based on a cna output and the dataset from which 
# that output was inferred.
# The dataset should be a crisp-set dataset in the form of a dataframe.
# Examples on how to use the functions are provided at the bottom of
# the script.



#################
### functions ###
#################

# asf_contrapositives(ana, df) calculates contrapositive consistency and
# contrapositive coverage for atomic solution formulas (asfs).
# The first argument, ana, is a cna output (as returned by the function cna()).
# The second argument, df, is the crisp-set dataset from which ana was 
# inferred using the cna function.
# The function returns a dataframe of asfs, as is standardly returned by the
# asf() function in the cna package, with two additional columns, one 
# containing the contrapositive consistency and one containing the 
# contrapositive coverage of the asfs in the dataframe.

asf_contrapositives <- function(ana, df){
  # get the asfs of the cna output
  asfs <- asf(ana)
  # make a copy
  asfs2 <- asfs
  # get the total number of cases in the dataset
  N <- nrow(df)
  # get the number of cases with a positive outcome value: |Y|
  asfs2$Y <- sapply(asfs2$outcome, function(y) sum(df[,y]))
  # get the number of cases with a positive antecedent value and positive 
  # outcome value: |X*Y| (coverage = |X*Y|/|Y|, so |X*Y| = coverage * |Y|)
  asfs2$XY <- asfs2$coverage * asfs2$Y
  # get the number of cases with a negative antecedent value and positive 
  # outcome value: |x*Y| (|Y| = |X*Y| + |x*Y|, so |x*Y| = |Y| - |X*Y|)
  asfs2$xY <- asfs2$Y - asfs2$XY
  # get the number of cases with a positive antecedent value: |X|
  # (coverage/consistency = |X|/|Y|, so coverage/consistency*|Y| = |X|)
  asfs2$X <- asfs2$coverage/asfs2$consistency*asfs2$Y
  # get the number of cases with a negative antecedent value: |x| (this is 
  # equal to the total number of cases in the dataset which do not have a 
  # positive antecedent value, so |x| = N - |X|)
  asfs2$x_small <- N - asfs2$X
  # get the number of cases with a negative outcome value: |y| (this is 
  # equal to the total number of cases in the dataset which do not have a 
  # positive outcome value, so |y| = N - |Y|)
  asfs2$y_small <- N - asfs2$Y
  # get the number of cases with a negative antecedent value and a negative 
  # outcome value: |x*y| (|x| = |x*Y| + |x*y|, so |x*y| = |x| - |x*Y|)
  asfs2$xy <- asfs2$x_small - asfs2$xY
  # get contrapositive consistency: |x*y|/|y|
  asfs$contrapositive_consistency <- asfs2$xy / asfs2$y_small
  # get contrapositive coverage: |x*y|/|x|
  asfs$contrapositive_coverage <- asfs2$xy / asfs2$x_small
  return(asfs)
}


# csf_contrapositives(ana, df) calculates contrapositive consistency and
# contrapositive coverage for complex solution formulas (csfs).
# The first argument, ana, is a cna output (as returned by the function cna()).
# The second argument, df, is the crisp-set dataset from which ana was 
# inferred using the cna function.
# The function returns a dataframe of csfs, as is standardly returned by the
# csf() function in the cna package, with two additional columns, one 
# containing the contrapositive consistency and one containing the 
# contrapositive coverage of the csfs in the dataframe.
#
# csf_contrapositives() depends on asf_contrapositives(). So, make sure to 
# execute the definition of both asf_contrapositives() and 
# csf_contrapositives() before using csf_contrapositives().

csf_contrapositives <- function(ana, df){
  # get the contrapositives of the asfs of the cna output
  asfs_contr <- asf_contrapositives(ana, df)
  # get the csfs of the cna output
  csfs <- csf(ana)
  # for each csf, get its contrapositive consistency by retrieving the 
  # contrapositive consistencies of all asfs in the csf (from asfs_contr) and 
  # then taking the minimum of those contrapositive consistencies
  csfs$contrapositive_consistency <- sapply(csfs$condition, function(x) {
    if (grepl("\\*", x)){
      parts <- strsplit(x, "\\)\\*\\(")[[1]]
      parts <- trimws(parts)
      parts <- gsub("\\(", "", parts)
      parts <- gsub("\\)", "", parts)
      parts <- sapply(parts, function(y) asfs_contr$contrapositive_consistency[asfs_contr$condition == y])
      return(min(parts))
    } else {
      return(asfs_contr$contrapositive_consistency[asfs_contr$condition == x])
      }})
  # for each csf, get its contrapositive coverage by retrieving the 
  # contrapositive coverages of all asfs in the csf and then taking the 
  # minimum of those contrapositive coverages
  csfs$contrapositive_coverage <- sapply(csfs$condition, function(x) {
    if (grepl("\\*", x)){
      parts <- strsplit(x, "\\)\\*\\(")[[1]]
      parts <- trimws(parts)
      parts <- gsub("\\(", "", parts)
      parts <- gsub("\\)", "", parts)
      parts <- sapply(parts, function(y) asfs_contr$contrapositive_coverage[asfs_contr$condition == y])
      return(min(parts))
      } else {
      return(asfs_contr$contrapositive_coverage[asfs_contr$condition == x])
        }
    })
  return(csfs)
}


# msc_contrapositives(ana, df) calculates contrapositive consistency and
# contrapositive coverage for minimally sufficient conditions (mscs).
# The first argument, ana, is a cna output (as returned by the function cna()).
# The second argument, df, is the crisp-set dataset from which ana was 
# inferred using the cna function.
# The function returns a dataframe of mscs, as is standardly returned by the
# msc() function in the cna package, with two additional columns, one 
# containing the contrapositive consistency and one containing the 
# contrapositive coverage of the mscs in the dataframe.

msc_contrapositives <- function(ana, df){
  # get the mscs of the cna output
  mscs <- msc(ana)
  # make a copy
  mscs2 <- mscs
  # get the total number of cases in the dataset
  N <- nrow(df)
  # get the number of cases with a positive outcome value: |Y|
  mscs2$Y <- sapply(mscs2$outcome, function(y) sum(df[,y]))
  # get the number of cases with a positive antecedent value and positive 
  # outcome value: |X*Y| (coverage = |X*Y|/|Y|, so |X*Y| = coverage * |Y|)
  mscs2$XY <- mscs2$coverage * mscs2$Y
  # get the number of cases with a negative antecedent value and positive 
  # outcome value: |x*Y| (|Y| = |X*Y| + |x*Y|, so |x*Y| = |Y| - |X*Y|)
  mscs2$xY <- mscs2$Y - mscs2$XY
  # get the number of cases with a positive antecedent value: |X|
  # (coverage/consistency = |X|/|Y|, so coverage/consistency*|Y| = |X|)
  mscs2$X <- mscs2$coverage/mscs2$consistency*mscs2$Y
  # get the number of cases with a negative antecedent value: |x| (this is 
  # equal to the total number of cases in the dataset which do not have a 
  # positive antecedent value, so |x| = N - |X|)
  mscs2$x_small <- N - mscs2$X
  # get the number of cases with a negative outcome value: |y| (this is 
  # equal to the total number of cases in the dataset which do not have a 
  # positive outcome value, so |y| = N - |Y|)
  mscs2$y_small <- N - mscs2$Y
  # get the number of cases with a negative antecedent value and a negative 
  # outcome value: |x*y| (|x| = |x*Y| + |x*y|, so |x*y| = |x| - |x*Y|)
  mscs2$xy <- mscs2$x_small - mscs2$xY
  # get contrapositive consistency: |x*y|/|y|
  mscs$contrapositive_consistency <- mscs2$xy / mscs2$y_small
  # get contrapositive coverage: |x*y|/|x|
  mscs$contrapositive_coverage <- mscs2$xy / mscs2$x_small
  return(mscs)
}



################
### examples ###
################

library(cna)

# Run a cna analysis.
ana <- cna(d.irrigate, con = 0.8, cov = 0.8)

# Calculate contrapositive consistency and contrapositive coverage of asfs,
# using the asf_contrapositives function.
asfs_contr <- asf_contrapositives(ana, d.irrigate)
asfs_contr

# Calculate contrapositive consistency and contrapositive coverage of csfs,
# using the csf_contrapositives function.
csfs_contr <- csf_contrapositives(ana, d.irrigate)
csfs_contr

# Calculate contrapositive consistency and contrapositive coverage of mscs,
# using the msc_contrapositives function.
mscs_contr <- msc_contrapositives(ana, d.irrigate)
mscs_contr


