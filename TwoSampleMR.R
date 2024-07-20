# Set Working Rirectory
setwd("/MR")

library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(mr.raps)
library(ieugwasr)
library(dplyr)
library(MendelianRandomization)

# exposure
FileNames <-list.files(paste0(getwd()),pattern=".tsv")  
exp_dat_ids <- FileNames
exps <- FileNames

# outcome
out<-fread("outcome.txt",header = T)
out$phenotype <- 'phenotypename'  
outcomeid <- out
rm(out)
head(outcomeid)

dir.create(path = "mendelian test") #Create a new folder to store the results

# Performing Mendelian randomization analyses
for (qaq in 1:length(exp_dat_ids)) { 
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  
  d3<- try(fread(paste0(getwd(),"/",FileNames[qaq]),fill=TRUE),silent = F)
  d3<-subset(d3,d3$pval<5e-8)  # Selecting exposure-related SNPs

  if (nrow(d3) > 0) {
    d3<-as.data.frame(d3)
    d3<-format_data(d3,
                    type="exposure")
    
    # Clumping
    library(ieugwasr)     
    
    d4<- ld_clump(
      clump_kb = 10000,
      clump_r2 = 0.001,    
      pop = "EUR",
      dplyr::tibble(rsid=d3$SNP, pval=d3$pval.exposure, id=d3$id.exposure),
      # set the path to the plink executable file
      plink_bin = "C:/Users/Wenda/Desktop/plink.exe",
      # set the path to the binary file for Plink
      bfile = "C:/Users/Wenda/Desktop/EUR/EUR"
    )
    exp_data<-subset(d3,SNP %in% d4$rsid) 

    # phenoscanner- Removing confounders
    confounder_list <- lapply(seq(1, nrow(exp_data), 100), function(start_index) {
      end_index <- min(start_index + 99, nrow(exp_data))
      batch_data <- exp_data[start_index:end_index, ]
      
      MendelianRandomization::phenoscanner(
        snpquery = unique(batch_data$SNP),
        catalogue = "GWAS",
        pvalue = 5e-8,
        proxies = "None",
        r2 = 0.8,
        build = 37
      )
    })
    
    confounder <- bind_rows(lapply(confounder_list, `[[`, "results"))
    if (nrow(confounder) > 0) {
      a=confounder[,c("snp","trait","p")];a
      a$trait <- tolower(a$trait)
      
      confounders <- "education|smoking|BMI|autoimmune disease|cortisol" 
      pattern <- paste0("^(", paste(confounders, collapse = "|"), ")$")
      
      c <- a %>%
        filter(grepl(pattern, trait))
      
      if (nrow(c) > 0) {
        c$phenotype <- FileNames[qaq]
        
        exp_data <- exp_data[!(exp_data$SNP %in% c$snp), ]
        openxlsx::write.xlsx(c,
                             file = paste0("mendelian test/", exp, "-confounders.xlsx"),
                             row.names = FALSE)
      }} else {
        exp_data <- exp_data
      } 
    
    if(length(exp_data[,1])>0){
      outcome_dat<-merge(exp_data,outcomeid,by.x = "SNP",by.y = "SNP")
      if(length(outcome_dat[,1])>0){  
        write.csv(outcome_dat,file = "d.csv")
        out_data <- read_outcome_data(
          snps = exp_data$SNP, 
          filename = "d.csv",
          sep = ",")
        
        # Removing outcome-related SNPs 
        out_data <- subset(out_data,pval.outcome>5e-8)
        
        # Harmonization
        dat <- TwoSampleMR::harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)
        
        # Removing the palindromic SNPs
        dat <-subset(dat,mr_keep==TRUE)
        
        # Calculating the F-statistic
        dat <- get_f(dat, F_value = 10)
        
        # Horizontal_pleiotropy by using MR_PRESSO    
        if (nrow(dat) > 3) {
        mr_presso_res <- mr_Presso(dat, num = 1000)
        mr_presso_main <- mr_presso_pval(mr_presso_res)
        dat <- mr_presso_snp(mr_presso_res, mr_presso_main, dat, type = "data")
        resMRPRESSO=mr_presso_res[["Main MR results"]]
        resMRPRESSO
        global_test_p <- mr_presso_res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
        se1=sqrt(((resMRPRESSO[1,3])^2)/qchisq(Pval_raw <- resMRPRESSO[1,6],1,lower.tail=F))
        se2=sqrt(((beta_cor <- resMRPRESSO[2,3])^2)/qchisq(Pval_cor <- resMRPRESSO[2,6],1,lower.tail=F))
        resMRPRESSO <- resMRPRESSO %>%
          dplyr::mutate(se = c(se1,se2))
        outliers <- dat$SNP[mr_presso_res[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]]
        outliers = as.data.frame(outliers)
        outliers <- paste(outliers$outliers, collapse = ", ")
        global_test_p = as.data.frame(global_test_p)
        resMRPRESSO
        TTT <- as.data.frame(resMRPRESSO)
        TTT
        openxlsx::write.xlsx(TTT, paste0("mendelian test/",exp,"-MR-PRESSO.xlsx"))
        } 
        dat <-subset(dat,mr_keep==TRUE)
        res <-choose_MR(dat=dat)
        res_hete <-purrr::map(.x=seq_along(res),
                              .f=~res[[.x]][[1]])
        res_hete <-do.call(rbind,res_hete)
        res_hete
        res1 <- generate_odds_ratios(res[[1]][[2]])
        res1
        res1$estimate <- paste0(
          format(round(res1$or, 2), nsmall = 2), " (", 
          format(round(res1$or_lci95, 2), nsmall = 2), "-",
          format(round(res1$or_uci95, 2), nsmall = 2), ")")
        res1
        print(paste0(exp,"_SNP_",res1$nsnp[1]))
        resdata <- dat
        openxlsx::write.xlsx(dat,file = paste0("mendelian test/",exp,"-dat.xlsx"), rowNames = FALSE)
        
        openxlsx::write.xlsx(res1,paste0("mendelian test/",exp,"-res.xlsx"))
        
        # MR Steiger test
        res_steiger <- steiger_test(dat) 
        
        # statistical power
        N <- dat$samplesize.outcome[1]
        alpha <- 0.05
        R2xz <- sum(dat$R2)
        K <- (dat$ncase.outcome[1] / dat$ncontrol.outcome[1])
        if (length(dat[, 1]) == 1) {
          OR <- res1 %>% filter(method == "Wald ratio") %>% pull(or)
        } else if (length(dat[, 1]) > 1) {
          OR <-
            res1 %>%filter(grepl("Inverse variance weighted", method)) %>%pull(or)
        }
        epower = NA
        power <- results_binary(N, alpha, R2xz, K, OR, epower)
        
        library(magrittr)
        # Main result 
        res3 <- res1#[1:5,]
        res3 <- res3[,-c(10:14)]
        res4 <- tidyr::pivot_wider(
          res3,names_from ="method",names_vary = "slowest",
          values_from = c("b","se","pval","estimate") )
        
        col_names <- colnames(res4)
        
        new_col_names <- gsub("\\(.*\\)", "", col_names)
        
        colnames(res4) <- new_col_names
        
        ##steiger
        res_steiger2 <- dplyr::select(res_steiger,
                                      correct_causal_direction,steiger_pval)
        #Power
        power2 <- tidyr::pivot_wider(
          power,names_from ="Parameter",names_vary = "slowest",
          values_from = c("Value","Description") )
        power2 <- power2[,1]
        
        # Merge
        res_ALL <- cbind(res4, res_steiger2,power2
        )
        res_ALL$F <- mean(dat$F,na.rm=TRUE)
        res_ALL$R2 <- sum(dat$R2)
       
        if (length(dat[, 1]) > 0 && length(dat[, 1]) <= 2) {
          write.csv(res_ALL, file = paste0("mendelian test/", exp, "1.csv"), row.names = FALSE)
        } else {
          
          # Cochran's Q Test and LOO
          res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
          res_leaveone <- mr_leaveoneout(dat)
          
          p1 <- mr_scatter_plot(res[[1]][[2]], dat)
          pdf(paste0("mendelian test/", exp, "_scatter.pdf"))
          print(p1[[1]])
          dev.off()
          
          res_single <- mr_singlesnp(dat,all_method)
          p2 <- mr_forest_plot(res_single)
          pdf(paste0("mendelian test/", exp, "_forest.pdf"))
          print(p2[[1]])
          dev.off()
          
          p3 <- mr_funnel_plot(res_single)
          pdf(paste0("mendelian test/", exp, "_funnel.pdf"))
          print(p3[[1]])
          dev.off()
          
          res_loo <- mr_leaveoneout(dat)
          pdf(paste0("mendelian test/", exp, "_leave_one_out.pdf"))
          print(mr_leaveoneout_plot(res_loo))
          dev.off()
          
          res_hete2 <- tidyr::pivot_wider(
            res_hete, names_from = "method", names_vary = "slowest",
            values_from = c("Q", "Q_df", "Q_pval")
          ) %>% 
            dplyr::select(-id.exposure, -id.outcome, -outcome, -exposure)
          res_hete2 <- res_hete2[, 4:6]
          
          # MR-Egger
          res_plei2 <- dplyr::select(res_plei, egger_intercept, se, pval)
          
          res_ALL <- cbind(res_ALL, res_hete2, res_plei2)
          
          if (nrow(dat) > 3) {
  
          res_ALL$outliers <- outliers
          res_ALL <- cbind(res_ALL, global_test_p)
          }
          write.csv(res_ALL, file = paste0("mendelian test/", exp, ".csv"), row.names = FALSE)
        }
      }
    }
  }
}

# Merge the results
csv_files <- list.files("/MR/mendelian test", pattern = "tsv.csv", full.names = TRUE)      
combined_df <- read.csv(csv_files[1], stringsAsFactors = FALSE)      

# Convert the 'exposure' column to character    
if (!is.null(combined_df$exposure)) {      
  combined_df$exposure <- as.character(combined_df$exposure)      
}   

# Convert the 'global_test_p' column to character  
if (!is.null(combined_df$global_test_p)) {      
  combined_df$global_test_p <- as.character(combined_df$global_test_p)      
}  

for (i in 1:length(csv_files)) {      
  temp_df <- read.csv(csv_files[i], stringsAsFactors = FALSE)      
  
  # Convert the 'exposure' column to character    
  if (!is.null(temp_df$exposure)) {      
    temp_df$exposure <- as.character(temp_df$exposure)      
  }    
  
  # Convert the 'global_test_p' column to character  
  if (!is.null(temp_df$global_test_p)) {      
    temp_df$global_test_p <- as.character(temp_df$global_test_p)      
    temp_df$global_test_p[temp_df$global_test_p == "<"] <- ""      
  }      
  
  combined_df <- bind_rows(combined_df, temp_df)      
}
write.csv(combined_df,"exposure-outcome.csv") 