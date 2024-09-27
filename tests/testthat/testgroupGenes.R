# Load test dataset
#library(testthat)
#library(alakazam)

data <- data.frame(sequence_id = c(unlist(lapply(1:30, function(x) paste0("seq", x)))),
                   subject_id = "S1",
                   v_call = c("IGHV1-1*01","IGKV1-1*01","IGHV1-1*01","IGKV1-1*01",
                              "IGHV1-1*01","IGKV2-1*01","IGHV1-1*01","IGKV1-1*01",
                              "IGKV2-1*01","IGHV1-1*01","IGHV1-1*01","IGHV1-1*01",
                              "IGKV1-1*01,IGKV2-1*01","IGHV1-1*01","IGKV3-1*01",
                              "IGHV2-1*01","IGKV1-1*01","IGKV1-1*01","IGKV1-1*01",
                              "IGHV1-1*01","IGHV1-1*01","IGHV1-1*01","IGKV4-1*01",
                              "IGHV1-1*01","IGKV4-1*01","IGHV1-1*01","IGHV1-1*01",
                              "IGKV1-1*01","IGHV1-1*01","IGKV1-1*01"),
                   j_call = c("IGHJ2*01","IGKJ1*01","IGHJ2*01","IGKJ1*01","IGHJ2*01",
                              "IGKJ1*01","IGHJ2*01","IGKJ1*01","IGKJ1*01","IGHJ2*01",
                              "IGHJ2*01","IGHJ2*01","IGKJ1*01","IGHJ2*01","IGKJ1*01",
                              "IGHJ1*01","IGKJ1*01","IGKJ1*01","IGKJ1*01","IGHJ2*01",
                              "IGHJ2*01","IGHJ2*01","IGKJ1*01","IGHJ2*01","IGKJ1*01",
                              "IGHJ2*01","IGHJ2*01","IGKJ1*01","IGHJ2*01","IGKJ1*01"),
                   junction = c("TGTAAAAAATGG","TGTCCCCCCTGG","TGTAAAAAATGG","TGTCCCCCCTGG",
                                "TGTAAAAAATGT","TGTCCCCCCTGG","TGTAAAAAATGG","TGTCCCCCCTGG",
                                "TGTCCCCCCTGG","TGTAAAAAATGG","TGTAAAAAATGG","TGTAAAAAATGG",
                                "TGTCCCCCCTGG","TGTAAAAAATCG","TGTCCCCCCTGG","TGTAAAACCTGG",
                                "TGTCCCCCCTGG","TGTCCCCCCTGG","TGTCCCCCCTGG","TGTAAAAAATGT",
                                "TGTAAAAAATCG","TGTAAAAAATGG","TGTCCCCCCTGG","TGTAAAAAATGG",
                                "TGTCCCCCCTGG","TGTAGAAAATGG","TGTAAAAAATGG","TGTCCCCCCTGG",
                                "TGTAAAAAATGG","TGTCCCCCCTGG"),
                   locus = c("IGH","IGK","IGH","IGK","IGH","IGK","IGH","IGK","IGK",
                             "IGH","IGH","IGH","IGK","IGH","IGK","IGH","IGK","IGK",
                             "IGK","IGH","IGH","IGH","IGK","IGH","IGK","IGH","IGH",
                             "IGK","IGH","IGK"),
                   cell_id = c(1,1,2,2,3,3,4,4,4,6,NA,8,8,5,5,7,7,9,NA,10,11,12,
                               12,13,13,14,15,15,16,16),
                   junction_length = 12,
                   clone_id = c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,NA,NA,1,2,4,4,
                                4,4,1,1,1,1,1),
                   expected_clone_subgroup = c(1,1,1,1,2,2,1,1,NA,1,1,1,1,1,1,1,
                                               1,NA,NA,2,1,1,1,1,1,1,1,1,1,1))

# run group genes and make sure that it works as expected

# # HEAVY only
# results <- groupGenes(data)
# expect_equal(results$vj_group, results$expected_vj_group_heavy_only)
# # with cell_id
# results <- groupGenes(data, cell_id = "cell_id")
# expect_equal(results$vj_group, results$expected_vj_group_heavy_only_cell_id)
# # should have dropped the ones without a cell id (seqs 11 and 19)
# expect_equal("seq11" %in% results$sequence_id, FALSE)
# expect_equal("seq19" %in% results$sequence_id, FALSE)
# 
# # Heavy only == FALSE
# # with cell_id
# results <- groupGenes(data, only_heavy = FALSE, cell_id = "cell_id")
# expect_equal(results$vj_group, results$expected_vj_group_heavy_false_cell_id)
# # should have dropped the ones without a cell id (seqs 11 and 19)
# expect_equal("seq11" %in% results$sequence_id, FALSE)
# expect_equal("seq19" %in% results$sequence_id, FALSE)
# # no cell_id
# results <- groupGenes(data, only_heavy = FALSE)
# expect_equal(results$vj_group, results$expected_vj_group_heavy_false_no_cell_id)
