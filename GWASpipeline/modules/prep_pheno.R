# Extract parameters from snakemake object.
phenotable <- snakemake@config$phenotable
genoprefix <- snakemake@config$genotypes
phenotype <- snakemake@params$phenotype
transformation <- snakemake@params$transformation
method <- snakemake@params$method
covariates <- snakemake@params$covariates
indvid <- snakemake@params$indvid
out_indvs <- snakemake@output$indvs
out_pheno <- snakemake@output$pheno
out_covar <- snakemake@output$covar

# Rank-based inverse normal transformation
qtrans <- function(x) {
  k <- !is.na(x)
  ran <- rank(x[k])
  y <- qnorm((1:sum(k) - 0.5) / sum(k))
  x[k] <- y[ran]
  x
}

# Load .fam order of individuals.
indvs <- data.frame(
    IID=as.character(data.table::fread(paste0(genoprefix, '.fam'), data.table = F)[,2]),
    intercept=1
)
colnames(indvs)[1] <- indvid

# Add sex if not covariate for QT.
colstokeep <- c(indvid, phenotype, covariates)
if(!'sex' %in% colstokeep){colstokeep <- c(colstokeep, 'sex')}

# Load phenotypes
pt <- data.table::fread(phenotable, data.table = F, select = colstokeep)

# Remove duplicated rows (indvs with multiple identical rows).
pt <- unique(pt)

# Left join phenotypes and covariants to indvs
indvs <- dplyr::left_join(indvs, pt)

# Only keep indvs with phenotype and covariates defined.
indvs <- indvs[rowSums(is.na(indvs[,c(phenotype, covariates)]))==0, ]

# Transformation
if(transformation=='QT'){
    if(!'sex' %in% colnames(indvs)){
        stop(paste0('Quantile transformation needs a sex covariate with 1s and 2s'))
    }
    indvs[indvs$sex==1, phenotype] <- qtrans(indvs[indvs$sex==1, phenotype])
    indvs[indvs$sex==2, phenotype] <- qtrans(indvs[indvs$sex==2, phenotype])
} else if(transformation=='UT'){
    dummy <- 'dummy' # for avoiding returning NULL
} else if(transformation=='BT') {
    # For binary transform, we assume that the value with the largest
    # number of counts is the not-measuremable threshold. All equal or
    # below is set to 0 and all above is set to 1.
    val_with_maxcount <- as.numeric(names(sort(table(indvs[,phenotype]), decreasing = T)[1]))
    print(paste0('Binary transform lower threshold: ', val_with_maxcount))
    indvs[,phenotype] <- ifelse(indvs[,phenotype]<=val_with_maxcount, 0, 1)
} else {
    stop(paste0("Tranformation '", transformation, "' not implemented yet."))
}

# One-hot encode the cohort variables, using model.matrix.
possible_cohort_vars <- c('cohort') # Manually add cohort-variable names used here.
if (any(possible_cohort_vars %in% covariates)){
    cohort_varname <- covariates[covariates %in% possible_cohort_vars]
    if (length(unique(indvs[,cohort_varname]))>1){ # If more than one cohort represented.
        cohort_onehot <- model.matrix(lm(as.formula(paste(phenotype, '~', cohort_varname)), data=indvs))
        indvs <- cbind(indvs, cohort_onehot)
        cohort_names <- colnames(cohort_onehot)[-1]
        covariates <- c(covariates[!covariates==cohort_varname], cohort_names)
    } else { # If only one cohort represented.
        covariates <- covariates[!covariates==cohort_varname]
    }
}

# Add needed columns to ordered individuals.
colstokeep <- c(indvid, phenotype, covariates)

# Write outputs.
write.table(indvs[,c(indvid, indvid)], out_indvs, sep='\t',
    row.names = F, quote = F, col.names = F)

if(method %in% c('gemma', 'gemmaloco')){
    write.table(indvs[,paste(phenotype)], out_pheno, row.names=F, quote = F, col.names = F)
    write.table(indvs[,c('intercept', covariates)], out_covar, row.names=F, quote = F, col.names = F)
} else {
    stop(paste0("Method '", method, "' not implemented."))
}