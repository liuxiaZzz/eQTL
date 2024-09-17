# .libPaths("/home/users/nus/e1124313/project/e1124313/env/PEER/lib/R/library")
# 一直memory error, 直到从github上找到了一个人分享的yaml，用conda env update -f xxx.yaml来安装，然后就可以了

library(peer)
# library(data.table)
# library(dplyr)

set.seed(1234)

compute_PEER <- function(cell_type){

    # Get expression file --------------------------------
    expr = read.csv(sprintf("/home/users/nus/e1124313/scratch/eqtl/0716_input/%s/%s.txt", cell_type, cell_type), sep="\t", header=T)
    dim(expr)
    rnames = as.character(expr$Sample_ID)
    expr2 = expr[, -which(names(expr) %in% c('Sample_ID', 'Batch_ID', 'Age', 'Sex', 'Ethnicity_1', 'Ethnicity_2', 'Status', 'Disease_status', 'Assay_method', 'CT_1', 'CT_2', 'geno_PC1', 'geno_PC2', 'geno_PC3', 'geno_PC4', 'geno_PC5', 'geno_PC6'))]
    expr2[1:5,1:5]
    cnames = colnames(expr2)

    # Get covariate files ---------------------------------
    covs = expr[c('Age', 'Sex', 'geno_PC1', 'geno_PC2', 'geno_PC3', 'geno_PC4', 'geno_PC5', 'geno_PC6')]

    # Set PEER paramaters based on the instructions from PEER package website

    model = PEER()

    PEER_setPhenoMean(model, as.matrix(expr2))
    print(sprintf("%s Pheno mean set successfully", cell_type))

    PEER_setCovariates(model, as.matrix(covs))
    print(sprintf("%s Covariates set successfully", cell_type))

    dim(PEER_getPhenoMean(model))

    # PEER_setAdd_mean(model, TRUE)

    # PEER_setNmax_iterations(model, 100)

    PEER_setNk(model,10) # Set to generate 10 PEER factors
    print(sprintf("%s PEER factors set successfully", cell_type))

    PEER_getNk(model)

    print(sprintf("%s Starting to run PEER", cell_type))
    PEER_update(model)

    # Set a directory for each cell type to save outputs
    dir.create(sprintf("", cellLabel))

    # Calculate and save the PEER factors
    factors = PEER_getX(model)
    dim(factors)
    factors_df = data.frame(factors)
    factors_df$sampleid = rnames
    factors_df = factors_df[c(19,1:18)]
    colnames(factors_df) = c("Sample_Id", "Age", "Sex", 'geno_PC1', 'geno_PC2', 'geno_PC3', 'geno_PC4', 'geno_PC5', 'geno_PC6', 
        "pf1", "pf2","pf3", "pf4", "pf5",
        "pf6", "pf7","pf8", "pf9", "pf10")
    combined_1 = cbind(factors_df, expr2)
    write.table(combined_1, file=sprintf("/home/users/nus/e1124313/scratch/eqtl/0716_input/%s/%s_add_peer.txt", cell_type, cell_type),
        sep="\t", col.names=T, row.names=F, quote=F)    

    # # Calculate and save the weights for each factor
    # weights = PEER_getW(model)
    # dim(weights)
    # weights_df = data.frame(weights)
    # weights_df$geneid = cnames
    # weights_df = weights_df[c(19,1:18)]
    # colnames(weights_df) = c("geneid", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age",
    #     "pf1", "pf2","pf3", "pf4", "pf5",
    #     "pf6", "pf7","pf8", "pf9", "pf10")
    # write.table(weights_df, file=sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/peer_factors/%s/%s_peer_factor_weights.tsv", cellLabel, cellLabel),
    #     sep="\t", col.names=T, quote=F, row.names=F)   

    # # Calculate and save the precision values
    # precision = PEER_getAlpha(model)
    # dim(precision)
    # precision_df = data.frame(precision)
    # precision_df$covariate = c("sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age",
    #     "pf1", "pf2","pf3", "pf4", "pf5",
    #     "pf6", "pf7","pf8", "pf9", "pf10")
    # write.table(precision_df, file=sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/peer_factors/%s/%s_peer_factor_precision.tsv", cellLabel, cellLabel),
    #     sep="\t", col.names=T, quote=F, row.names=F)   

    # Calculate and save the residuals
    residuals = PEER_getResiduals(model)
    dim(residuals)
    residuals_df = data.frame(residuals)
    colnames(residuals_df) = cnames
    combined_2 = cbind(factors_df, residuals_df)
    write.table(combined_2, file=sprintf("/home/users/nus/e1124313/scratch/eqtl/0716_input/%s/%s_residual_peer.txt", cell_type, cell_type),
        sep="\t", col.names=T, quote=F, row.names=F)
}

# # Session info ------------------------------------------
# print_session(here(output))
# cell_type_df = read.csv("/home/users/nus/e1124313/scratch/eqtl/input/cell_types.csv", sep="\t", header=F)
# for (i in 1:nrow(cell_type_df)){
#     cell_type = as.character(cell_type_df[i,1])
#     print(cell_type)
#     compute_PEER(cell_type)
# }
compute_PEER("CytoT_GZMH+")
compute_PEER("B_atypical")
compute_PEER('B_mem')
compute_PEER('B_naive')
compute_PEER('CytoT_GZMK+')
compute_PEER('NK_bright')
compute_PEER('NK_dim')
compute_PEER('Progen')
compute_PEER('T_mait')
compute_PEER('T4_em')
compute_PEER('T4_naive')
compute_PEER('T4_reg')
compute_PEER('T8_naive')