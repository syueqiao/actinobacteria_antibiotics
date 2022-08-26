new_ssmd_spline <- function(turi2_df_bac_DMSO, turi2_df_library_DMSO){
  
  turi2_df_bac_DMSO_splines <- select(turi2_df_bac_DMSO_spread, c(6:13))
  turi2_bac_time_DMSO <- c(0:7)
  turi2_bac_lm_DMSO <- apply(turi2_df_bac_DMSO_splines, 1, function(x) lm(x ~ lspline(turi2_bac_time_DMSO, knots = knots)))
  turi2_df_bac_DMSO_splines_coeff <- as.data.frame(unlist(lapply(turi2_bac_lm_DMSO, function(x) coef(x)[2])))
  
  turi2_df_bac_DMSO_spread$spline <- turi2_df_bac_DMSO_splines_coeff
  turi2_df_bac_DMSO_spline_vector <- turi2_df_bac_DMSO_spread$spline
  
  turi2_df_bac_in_DMSO_avg_rep1 <- turi2_df_bac_DMSO_spread %>%
    filter(rep == "rep1") %>%
    .$spline %>%
    mean
  
  turi2_df_bac_in_DMSO_avg_rep2 <- turi2_df_bac_DMSO_spread %>%
    filter(rep == "rep2") %>%
    .$spline %>%
    mean
  
  
  #combine
  turi2_df_bac_DMSO_spread$spline <- turi2_df_bac_DMSO_splines_coeff
  turi2_df_bac_DMSO_spline_vector <- turi2_df_bac_DMSO_spread$spline
  
  turi2_df_bac_in_DMSO_avg_rep3 <- turi2_df_bac_DMSO_spread %>%
    filter(rep == "rep3") %>%
    .$spline %>%
    mean
  
  turi2_df_library_DMSO_splines <- select(turi2_df_library_DMSO_spread, c(6:13))
  turi2_library_time_DMSO <- c(0:7)
  turi2_library_lm_DMSO <- apply(turi2_df_library_DMSO_splines, 1, function(x) lm(x ~ lspline(turi2_library_time_DMSO, knots = knots)))
  turi2_df_library_DMSO_splines_coeff <- as.data.frame(unlist(lapply(turi2_library_lm_DMSO, function(x) coef(x)[2])))
  
  turi2_df_library_DMSO_spread$spline <- turi2_df_library_DMSO_splines_coeff
  turi2_df_library_DMSO_spline_vector <- turi2_df_library_DMSO_spread$spline  
 
   turi2_df_lib_in_DMSO_avg_rep1 <- turi2_df_library_DMSO_spread %>%
    filter(rep == "rep1")
  turi2_df_lib_in_DMSO_avg_rep1$d_1 <- turi2_df_lib_in_DMSO_avg_rep1$spline -  turi2_df_bac_in_DMSO_avg_rep1
  
  turi2_df_lib_in_DMSO_avg_rep2 <- turi2_df_library_DMSO_spread %>%
    filter(rep == "rep2")
  turi2_df_lib_in_DMSO_avg_rep2$d_1 <- turi2_df_lib_in_DMSO_avg_rep2$spline -  turi2_df_bac_in_DMSO_avg_rep2
  
  turi2_df_lib_in_DMSO_avg_rep3 <- turi2_df_library_DMSO_spread %>%
    filter(rep == "rep3")
  turi2_df_lib_in_DMSO_avg_rep3$d_1 <- turi2_df_lib_in_DMSO_avg_rep3$spline -  turi2_df_bac_in_DMSO_avg_rep3
  
  turi2_df_lib_in_DMSO_avg_repall <- rbind (turi2_df_lib_in_DMSO_avg_rep1, turi2_df_lib_in_DMSO_avg_rep2, turi2_df_lib_in_DMSO_avg_rep3)
  turi2_df_lib_in_DMSO_avg_repall_avg <- group_by(turi2_df_lib_in_DMSO_avg_repall, well) %>%
    summarize(d_calc = mean(d_1)) %>%
    left_join(., turi2_df_lib_in_DMSO_avg_repall, by = "well")
  
  turi2_df_lib_in_DMSO_avg_repall_avg$s_1 <- (turi2_df_lib_in_DMSO_avg_repall_avg$d_1 - turi2_df_lib_in_DMSO_avg_repall_avg$d_calc)^2
  
  turi2_df_lib_in_DMSO_avg_repall_s_2 <- group_by(turi2_df_lib_in_DMSO_avg_repall_avg, well)%>%
    summarize(s_3 = sum(s_1))
  
  turi2_df_lib_in_DMSO_avg_repall_s_2$s_calc <- sqrt(turi2_df_lib_in_DMSO_avg_repall_s_2$s_3/2)
  turi2_df_lib_in_DMSO_avg_repall_avg <- left_join(turi2_df_lib_in_DMSO_avg_repall_s_2, turi2_df_lib_in_DMSO_avg_repall_avg, by = "well")
  turi2_df_lib_in_DMSO_avg_repall_avg$mm <- turi2_df_lib_in_DMSO_avg_repall_avg$d_calc/turi2_df_lib_in_DMSO_avg_repall_avg$s_calc
  turi2_df_lib_in_DMSO_avg_repall_avg <- turi2_df_lib_in_DMSO_avg_repall_avg[!duplicated(turi2_df_lib_in_DMSO_avg_repall_avg$well),]
  return(turi2_df_lib_in_DMSO_avg_repall_avg)
}

turi2_DMSO_ssmd_spline <- new_ssmd_spline(turi2_df_bac_DMSO, turi2_df_library_DMSO)
turi2_MeOH_ssmd_spline <- new_ssmd_spline(turi2_df_bac_MeOH, turi2_df_library_MeOH)
turi2_H2O_ssmd_spline<- new_ssmd_spline(turi2_df_bac_H2O, turi2_df_library_H2O)

############
turi2_DMSO_ssmd_spline$color_2 <- ifelse(turi2_DMSO_ssmd_spline$mm <= -1.645, "inhib", 
                                          ifelse(turi2_DMSO_ssmd_spline$mm >=1.645, "enh",
                                                 ifelse(turi2_DMSO_ssmd_spline$mm > -1.645 & turi2_DMSO_ssmd_spline$mm < 1.645, "noef", "NA")))
turi2_H2O_ssmd_spline$color_2 <- ifelse(turi2_H2O_ssmd_spline$mm <= -1.645, "inhib", 
                                         ifelse(turi2_H2O_ssmd_spline$mm >=1.645, "enh",
                                                ifelse(turi2_H2O_ssmd_spline$mm > -1.645 & turi2_H2O_ssmd_spline$mm < 1.645, "noef", "NA")))
turi2_MeOH_ssmd_spline$color_2 <- ifelse(turi2_MeOH_ssmd_spline$mm <= -1.645, "inhib", 
                                          ifelse(turi2_MeOH_ssmd_spline$mm >=1.645, "enh",
                                                 ifelse(turi2_MeOH_ssmd_spline$mm > -1.645 & turi2_MeOH_ssmd_spline$mm < 1.645, "noef", "NA")))
#######
turi2_DMSO_ssmd_spline <- left_join(turi2_DMSO_ssmd_spline, platemap_DMSO_merge, by = "well")
turi2_DMSO_ssmd_spline <- select(turi2_DMSO_ssmd_spline, well, mm, color_2, solv)
colnames(turi2_DMSO_ssmd_spline) <- c("well", "mm", "color_2", "solv")

turi2_MeOH_ssmd_spline <- left_join(turi2_MeOH_ssmd_spline, platemap_MeOH_merge, by = "well")
turi2_MeOH_ssmd_spline <- select(turi2_MeOH_ssmd_spline, well, mm, color_2, solv)

turi2_H2O_ssmd_spline <- left_join(turi2_H2O_ssmd_spline, platemap_H2O_merge_form, by = "well")
turi2_H2O_ssmd_spline <- select(turi2_H2O_ssmd_spline, well, mm, color_2, solv)

#combine
test_turi2_ssmd_spline <- rbind(turi2_DMSO_ssmd_spline, turi2_H2O_ssmd_spline, turi2_MeOH_ssmd_spline)
table(test_turi2_ssmd_spline$color_2)
colnames(test_turi2_ssmd_spline) <- c("well", "mm", "spline_cat", "solv")


#combine with the OD differences

turi2_DMSO_ssmd_diff_cat <- turi2_DMSO_ssmd
turi2_MeOH_ssmd_diff_cat <- turi2_MeOH_ssmd
turi2_H2O_ssmd_diff_cat <- turi2_H2O_ssmd


turi2_DMSO_ssmd_diff_cat$diff_cat <- ifelse(turi2_DMSO_ssmd_diff_cat$mm <= -1.645, "inhib", 
                                    ifelse(turi2_DMSO_ssmd_diff_cat$mm >=1.645, "enh",
                                           ifelse(turi2_DMSO_ssmd_diff_cat$mm > -1.645 & turi2_DMSO_ssmd_diff_cat$mm < 1.645, "noef", "NA")))
turi2_H2O_ssmd_diff_cat$diff_cat <- ifelse(turi2_H2O_ssmd_diff_cat$mm <= -1.645, "inhib", 
                                         ifelse(turi2_H2O_ssmd_diff_cat$mm >=1.645, "enh",
                                                ifelse(turi2_H2O_ssmd_diff_cat$mm > -1.645 & turi2_H2O_ssmd_diff_cat$mm < 1.645, "noef", "NA")))
turi2_MeOH_ssmd_diff_cat$diff_cat <- ifelse(turi2_MeOH_ssmd_diff_cat$mm <= -1.645, "inhib", 
                                          ifelse(turi2_MeOH_ssmd_diff_cat$mm >=1.645, "enh",
                                                 ifelse(turi2_MeOH_ssmd_diff_cat$mm > -1.645 & turi2_MeOH_ssmd_diff_cat$mm < 1.645, "noef", "NA")))


#correlate to chemical names from regex
turi2_DMSO_ssmd_diff_cat <- left_join(turi2_DMSO_ssmd_diff_cat, platemap_DMSO_merge, by = "well")
turi2_DMSO_ssmd_diff_cat <- select(turi2_DMSO_ssmd_diff_cat, well, mm, diff_cat, solv)
colnames(turi2_DMSO_ssmd_diff_cat) <- c("well", "mm", "diff_cat", "solv")

turi2_MeOH_ssmd_diff_cat <- left_join(turi2_MeOH_ssmd_diff_cat, platemap_MeOH_merge, by = "well")
turi2_MeOH_ssmd_diff_cat <- select(turi2_MeOH_ssmd_diff_cat, well, mm, diff_cat, solv)

turi2_H2O_ssmd_diff_cat <- left_join(turi2_H2O_ssmd_diff_cat, platemap_H2O_merge_form, by = "well")
turi2_H2O_ssmd_diff_cat <- select(turi2_H2O_ssmd_diff_cat, well, mm, diff_cat, solv)

#combine
test_turi2_diff_cat <- rbind(turi2_DMSO_ssmd_diff_cat, turi2_H2O_ssmd_diff_cat, turi2_MeOH_ssmd_diff_cat)
table(test_turi2_diff_cat$diff_cat)
colnames(test_turi2_diff_cat) <- c("well", "mm", "diff_cat", "solv")

#combine categories

test_turi2_ssmd_cat <- left_join(test_turi2_diff_cat, test_turi2_ssmd_spline, by = "solv")

test_turi2_ssmd_cat$category <- ifelse(test_turi2_ssmd_cat$diff_cat == "noef" & test_turi2_ssmd_cat$spline_cat == "noef", 5, 
                                                      ifelse(test_turi2_ssmd_cat$diff_cat == "inhib" & test_turi2_ssmd_cat$spline_cat == "inhib", 9, 
                                                             ifelse(test_turi2_ssmd_cat$diff_cat == "enh" & test_turi2_ssmd_cat$spline_cat == "enh", 1, 
                                                                    ifelse(test_turi2_ssmd_cat$diff_cat == "enh" & test_turi2_ssmd_cat$spline_cat == "noef", 2, 
                                                                           ifelse(test_turi2_ssmd_cat$diff_cat == "enh" & test_turi2_ssmd_cat$spline_cat == "inhib", 3, 
                                                                                  ifelse(test_turi2_ssmd_cat$diff_cat == "noef" & test_turi2_ssmd_cat$spline_cat == "enh", 4,
                                                                                         ifelse(test_turi2_ssmd_cat$diff_cat == "noef" & test_turi2_ssmd_cat$spline_cat == "inhib", 6, 
                                                                                                ifelse(test_turi2_ssmd_cat$diff_cat == "inhib" & test_turi2_ssmd_cat$spline_cat == "enh", 7,
                                                                                                       ifelse(test_turi2_ssmd_cat$diff_cat == "inhib" & test_turi2_ssmd_cat$spline_cat == "noef",9, 0)))))))))
table(test_turi2_ssmd_cat$category)
#quick histogram visualization/comparision
histo_turi2_DMSO_ssmd_cat <- left_join(turi2_DMSO_ssmd_diff_cat, turi2_DMSO_ssmd_spline, by = "solv")
histo_turi2_MeOH_ssmd_cat <- left_join(turi2_MeOH_ssmd_diff_cat, turi2_MeOH_ssmd_spline, by = "solv")
histo_turi2_H2O_ssmd_cat <- left_join(turi2_H2O_ssmd_diff_cat, turi2_H2O_ssmd_spline, by = "solv")
colnames(histo_turi2_DMSO_ssmd_cat) <- c("well.x", "mm.x", "diff_cat", "solv", "well.y", "mm.y", "spline_cat")
colnames(histo_turi2_MeOH_ssmd_cat) <- c("well.x", "mm.x", "diff_cat", "solv", "well.y", "mm.y", "spline_cat")
colnames(histo_turi2_H2O_ssmd_cat) <- c("well.x", "mm.x", "diff_cat", "solv", "well.y", "mm.y", "spline_cat")


###Category function because why not

category_func <- function(df){
  df$category <- ifelse(df$diff_cat == "noef" & df$spline_cat == "noef", 5, 
                                         ifelse(df$diff_cat == "inhib" & df$spline_cat == "inhib", 9, 
                                                ifelse(df$diff_cat == "enh" & df$spline_cat == "enh", 1, 
                                                       ifelse(df$diff_cat == "enh" & df$spline_cat == "noef", 2, 
                                                              ifelse(df$diff_cat == "enh" & df$spline_cat == "inhib", 3, 
                                                                     ifelse(df$diff_cat == "noef" & df$spline_cat == "enh", 4,
                                                                            ifelse(df$diff_cat == "noef" & df$spline_cat == "inhib", 6, 
                                                                                   ifelse(df$diff_cat == "inhib" & df$spline_cat == "enh", 7,
                                                                                          ifelse(df$diff_cat == "inhib" & df$spline_cat == "noef", 8, 0)))))))))
hist(df$category)
return(df)
  }

histograph_turi2_DMSO_ssmd_cat <- category_func(histo_turi2_DMSO_ssmd_cat)
histograph_turi2_MeOH_ssmd_cat <- category_func(histo_turi2_MeOH_ssmd_cat)
histograph_turi2_H2O_ssmd_cat <- category_func(histo_turi2_H2O_ssmd_cat)

histograph_turi2_all_ssmd_cat <- rbind(histograph_turi2_DMSO_ssmd_cat, histograph_turi2_MeOH_ssmd_cat, histograph_turi2_H2O_ssmd_cat)
setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics")

write.csv(histograph_turi2_all_ssmd_cat, "histograph_turi2_all_ssmd_cat.csv")
