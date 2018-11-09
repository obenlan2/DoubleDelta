Read_Combine_Quality <- function(File1, File2 = FALSE, File3 = FALSE, New_name.csv){
  if (File3 != FALSE){
    Data1 <- read.csv(File1)
    Data2 <- read.csv(File2)
    Data3 <- read.csv(File3)
    All_data <- rbind(Data1, Data2, Data3)
    Remove_undetermined <- All_data[which(All_data$Ct != "Undetermined"),]
    if(any(Remove_undetermined$NAW == 'TRUE')){
      warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
      write.csv(Remove_undetermined, New_name.csv)
      return(Remove_undetermined)
    }
    
  } else {
    if (File2 != FALSE) {
      Data1 <- read.csv(File1)
      Data2 <- read.csv(File2)
      All_data <- rbind(Data1, Data2)
      Remove_undetermined <- All_data[which(All_data$Ct != "Undetermined"),]
      if(any(Remove_undetermined$NAW == 'TRUE')){
        warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
        write.csv(Remove_undetermined, New_name.csv)
        return(Remove_undetermined)
      }
    }
  }
  
  Data <- read.csv(File1)
  Remove_undetermined <- Data[which(Data$Ct != "Undetermined"),]
  if(any(Remove_undetermined$NAW == 'TRUE')){
    warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
    write.csv(Remove_undetermined, New_name.csv)
    return(Remove_undetermined)
  }
}

deltas <- function(data, gene_name, time_point, treatment1, treatment2, control_gene, control_gene2 = FALSE, control_gene3 = FALSE){
  #outcon argument prob go here outcon <- file(NAME_of_File_argument, open = "w") -default vaulue should be false for arugment
  Find_All_Gene_Samples <- grep(gene_name, data$Sample.Name)
  All_Gene_Samples <- data[Find_All_Gene_Samples,]
  Find_Specific_Time <- grep(time_point, All_Gene_Samples$Sample.Name)
  All_Gene_Samples_At_Specific_Time <- All_Gene_Samples[Find_Specific_Time,]
  Find_All_treatment1 <- grep(treatment1, All_Gene_Samples_At_Specific_Time$Sample.Name)
  All_treatment1_Samples <- All_Gene_Samples_At_Specific_Time[Find_All_treatment1,]
  Find_All_treatment2 <-grep(treatment2, All_Gene_Samples_At_Specific_Time$Sample.Name)
  All_treatment2_Samples <- All_Gene_Samples_At_Specific_Time[Find_All_treatment2,]
  Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), NA == "Undetermined"),
                                 na.rm = TRUE)
  Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), NA == "Undetermined"),
                                 na.rm = TRUE)
  Delta1 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)
  print("This is the delta of treatments for your gene of interest")
  print(Delta1)
  if (control_gene3 != FALSE) {
    Find_All_Control_Gene_Samples1 <- grep(control_gene, data$Sample.Name)
    All_Control_Gene_Samples1 <- data[Find_All_Control_Gene_Samples1,]
    Find_All_Control_Gene_Samples2 <- grep(control_gene2, data$Sample.Name)
    All_Control_Gene_Samples2 <- data[Find_All_Control_Gene_Samples2,]
    Find_All_Control_Gene_Samples3 <- grep(control_gene3, data$Sample.Name)
    All_Control_Gene_Samples3 <- data[Find_All_Control_Gene_Samples3,]
    All_Control_Gene_Samples <- rbind(All_Control_Gene_Samples1, All_Control_Gene_Samples2, 
                                      All_Control_Gene_Samples3)
    Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)
    All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time,]
    Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
    All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1,]
    Find_All_treatment2 <-grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
    All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2,]
    Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), NA == "Undetermined"),
                                   na.rm = TRUE)
    Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), NA == "Undetermined"),
                                   na.rm = TRUE)
    Delta2 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)
    print("This is the delta of treatments for your control gene")
    print(Delta2)
    Double_Delta <- abs(Delta1 - Delta2)
    print("this is your double delta Ct", Double_Delta)
    print(Double_Delta)
    Fold_induction <- 2^Double_Delta
    print("This is the fold induction for your gene of interest")
    print(Fold_induction)
    return(Fold_induction)
  } else {
    if (control_gene2 != FALSE) {
      Find_All_Control_Gene_Samples1 <- grep(control_gene, data$Sample.Name)
      All_Control_Gene_Samples1 <- data[Find_All_Control_Gene_Samples1,]
      Find_All_Control_Gene_Samples2 <- grep(control_gene2, data$Sample.Name)
      All_Control_Gene_Samples2 <- data[Find_All_Control_Gene_Samples2,]
      All_Control_Gene_Samples <- rbind(All_Control_Gene_Samples1, All_Control_Gene_Samples2)
      Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)
      All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time,]
      Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
      All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1,]
      Find_All_treatment2 <-grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
      All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2,]
      Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), NA == "Undetermined"),
                                     na.rm = TRUE)
      Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), NA == "Undetermined"),
                                     na.rm = TRUE)
      Delta2 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)
      print("This is the delta of treatments for your control gene")
      print(Delta2)
      Double_Delta <- abs(Delta1 - Delta2)
      print("this is your double delta Ct", Double_Delta)
      print(Double_Delta)
      Fold_induction <- 2^Double_Delta
      print("This is the fold induction for your gene of interest")
      print(Fold_induction)
      return(Fold_induction) #return will end your function
    }
    Find_All_Control_Gene_Samples <- grep(control_gene, data$Sample.Name)
    All_Control_Gene_Samples <- data[Find_All_Control_Gene_Samples,]
    Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)
    All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time,]
    Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
    All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1,]
    Find_All_treatment2 <-grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
    All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2,]
    Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), NA == "Undetermined"),
                                   na.rm = TRUE)
    Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), NA == "Undetermined"),
                                   na.rm = TRUE)
    Delta2 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)
    print("This is the delta of treatments for your control gene(s)")
    print(Delta2)
    Double_Delta <- abs(Delta1 - Delta2)
    print("this is your double delta Ct", Double_Delta)
    print(Double_Delta)
    Fold_induction <- 2^Double_Delta
    print("This is the fold induction for your gene of interest")
    print(Fold_induction)
    return(Fold_induction)
  }
}

combine_folds <- function(f1, new_file, gene_names){ #new_file and gene_names need to be in quotes. 
  outcon <- file(new_file, open = "w")
  Fold_Induction <- as.character(f1)
  writeLines(Fold_Induction, outcon)
  close(outcon)
  read <- read.csv(new_file, header = FALSE)
  read$Gene_Name <- c(gene_names)
  read$Fold_change <- read$V1
  read2 <- read[,-1]
  write.csv(read2,new_file, row.names = FALSE)
  return(read2)
}