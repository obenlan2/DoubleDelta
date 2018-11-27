Read_Combine_Quality <- function(File1, File2 = FALSE, File3 = FALSE, New_name.csv){
  if (File3 != FALSE){ #lines 3-11 will only be executed if the user defines the File1, File2 and File3 arguments 
    Data1 <- read.csv(File1) #reads in the first CSV file
    Data2 <- read.csv(File2) #reads in the second CSV file
    Data3 <- read.csv(File3 )#reads in the third CSV file
    All_data <- rbind(Data1, Data2, Data3) #combines the rows of Data1, Data2, and Data3
    Remove_undetermined <- All_data[which(All_data$Ct != "Undetermined"),] #indexes the results of All_data by Ct values that don't equal "Undetermined"
    write.csv(Remove_undetermined, New_name.csv) #creates a new csv file with a name given by the user 
    if(any(Remove_undetermined$NAW == 'TRUE')){#lines 8-9 send a warning message if NAWs were detected
      warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
    }
    return(Remove_undetermined) #returns a data frame
  } else {
    if (File2 != FALSE) {#lines 16-23 will only be executed if the user defines the File1 and File2 arguments
      Data1 <- read.csv(File1)#reads in the first CSV file
      Data2 <- read.csv(File2)#reads in the second CSV file
      All_data <- rbind(Data1, Data2)#combines the rows of Data1 and Data2 
      Remove_undetermined <- All_data[which(All_data$Ct != "Undetermined"),]#indexes the results of All_data by Ct values that don't equal "Undetermined"
      write.csv(Remove_undetermined, New_name.csv)#creates a new csv file with a name given by the user
      if(any(Remove_undetermined$NAW == 'TRUE')){#lines 20-21 send a warning message if NAWs were detected
        warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
      }
      return(Remove_undetermined)#returns a data frame
    }
  }
  #lines 28-33 will only be executed if the user defines onle the File1 argument
  Data <- read.csv(File1)#reads in the first CSV file
  Remove_undetermined <- Data[which(Data$Ct != "Undetermined"),]#indexes the results of All_data by Ct values that don't equal "Undetermined"
  write.csv(Remove_undetermined, New_name.csv)#creates a new csv file with a name given by the user
  if(any(Remove_undetermined$NAW == 'TRUE')){#lines 30-31 send a warning message if NAWs were detected
    warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
  }
  return(Remove_undetermined)#returns a data frame
}

Data <- read.csv("C:/Users/obenlan2/Documents/R Project Functions/qRT-PCR_Results.csv")#reads in the first CSV file
Remove_undetermined <- Data[which(Data$Ct != "Undetermined"),]#indexes the results of All_data by Ct values that don't equal "Undetermined"
write.csv(Remove_undetermined, "New_name.csv")#creates a new csv file with a name given by the user
return(Remove_undetermined)#returns a data frame
if(any(Remove_undetermined$NAW == 'TRUE')){#lines 30-31 send a warning message if NAWs were detected
  warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
}
return(Remove_undetermined)#returns a data frame

deltas <- function(data, gene_name, time_point, treatment1, treatment2, control_gene, control_gene2 = FALSE, control_gene3 = FALSE){
  Find_All_Gene_Samples <- grep(gene_name, data$Sample.Name) #finds all samples matching gene_name
  All_Gene_Samples <- data[Find_All_Gene_Samples,] #indexes by samples matching gene_name
  Find_Specific_Time <- grep(time_point, All_Gene_Samples$Sample.Name) #finds data for a specified time point in All_Gene_Samples
  All_Gene_Samples_At_Specific_Time <- All_Gene_Samples[Find_Specific_Time,] #indexes by time in All_Gene_Samples
  Find_All_treatment1 <- grep(treatment1, All_Gene_Samples_At_Specific_Time$Sample.Name)#finds samples that match treatment1 
  All_treatment1_Samples <- All_Gene_Samples_At_Specific_Time[Find_All_treatment1,]#indexes by samples that match treatment1
  Find_All_treatment2 <-grep(treatment2, All_Gene_Samples_At_Specific_Time$Sample.Name)#finds samples that match treatment2 
  All_treatment2_Samples <- All_Gene_Samples_At_Specific_Time[Find_All_treatment2,]#indexes by samples that match treatment2
  Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), NA == "Undetermined"),
                                 na.rm = TRUE)#finds the mean Ct for all samples that had treatment1
  Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), NA == "Undetermined"),
                                 na.rm = TRUE)#finds the mean Ct for all samples that had treatment2
  Delta1 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)#finds the absolute value of the difference between Mean_CT_for_treatment1 and Mean_CT_for_treatment2
  print("This is the delta of treatments for your gene of interest")#prints text
  print(Delta1)#prints the delta for your gene of interest
  if (control_gene3 != FALSE) {#lines 54-81 will only be executed if the user defines the control_gene, control_gene2 and control_gene3 arguments
    Find_All_Control_Gene_Samples1 <- grep(control_gene, data$Sample.Name)#finds all samples matching control_gene
    All_Control_Gene_Samples1 <- data[Find_All_Control_Gene_Samples1,]#indexes by samples matching control_gene
    Find_All_Control_Gene_Samples2 <- grep(control_gene2, data$Sample.Name)#finds all samples matching control_gene2
    All_Control_Gene_Samples2 <- data[Find_All_Control_Gene_Samples2,]#indexes by samples matching control_gene2
    Find_All_Control_Gene_Samples3 <- grep(control_gene3, data$Sample.Name)#indexes by samples matching control_gene3
    All_Control_Gene_Samples3 <- data[Find_All_Control_Gene_Samples3,]#indexes by samples matching control_gene3
    All_Control_Gene_Samples <- rbind(All_Control_Gene_Samples1, All_Control_Gene_Samples2, 
                                      All_Control_Gene_Samples3)#combines the rows of All_Control_Gene_Samples1,All_Control_Gene_Samples2, and All_Control_Gene_Samples3
    Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)#finds data for a specified time point in All_Control_Gene_Samples
    All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time,]#indexes by time in All_Control_Gene_Samples
    Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)#finds samples that match treatment1 in All_Control_Gene_Samples
    All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1,]#indexes by treatment1 in All_Control_Gene_Samples
    Find_All_treatment2 <-grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)#finds samples that match treatment2 in All_Control_Gene_Samples
    All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2,]#indexes by treatment2 in All_Control_Gene_Samples
    Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), NA == "Undetermined"),
                                   na.rm = TRUE)#finds mean Ct for all samples that had treatement1
    Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), NA == "Undetermined"),
                                   na.rm = TRUE)#finds the mean Ct for all samples that had treatment2
    Delta2 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)#finds the absolute value of the difference between Mean_CT_for_treatment1 and Mean_CT_for_treatment2
    print("This is the delta of treatments for your control gene")#prints text
    print(Delta2)#prints the delta for your control genes
    Double_Delta <- abs(Delta1 - Delta2)#finds the absolute value of the difference between Delta1 and Delta2
    print("this is your double delta Ct")#prints text
    print(Double_Delta)#prints Double_Delta
    Fold_induction <- 2^Double_Delta #Finds the fold change
    print("This is the fold change for your gene of interest") #prints text
    print(Fold_induction)#prints Fold_induction
    return(Fold_induction)#returns Fold_induction
  } else {
    if (control_gene2 != FALSE) {#lines 84-108 will only be executed if the user only defines the control_gene and control_gene2 arguments
      Find_All_Control_Gene_Samples1 <- grep(control_gene, data$Sample.Name)#finds all samples matching control_gene
      All_Control_Gene_Samples1 <- data[Find_All_Control_Gene_Samples1,]#indexes by samples matching control_gene
      Find_All_Control_Gene_Samples2 <- grep(control_gene2, data$Sample.Name)#finds all samples matching control_gene2
      All_Control_Gene_Samples2 <- data[Find_All_Control_Gene_Samples2,]#indexes by samples matching control_gene2
      All_Control_Gene_Samples <- rbind(All_Control_Gene_Samples1, All_Control_Gene_Samples2)#combines the rows of All_Control_Gene_Samples1 and All_Control_Gene_Samples2,
      Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)#finds data for a specified time point in All_Control_Gene_Samples
      All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time,]#indexes by time in All_Control_Gene_Samples
      Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)#finds samples that match treatment1 in All_Control_Gene_Samples
      All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1,]#indexes by treatment1 in All_Control_Gene_Samples
      Find_All_treatment2 <-grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)#finds samples that match treatment2 in All_Control_Gene_Samples
      All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2,]#indexes by treatment2 in All_Control_Gene_Samples
      Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), NA == "Undetermined"),
                                     na.rm = TRUE)#finds mean Ct for all samples that had treatement1
      Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), NA == "Undetermined"),
                                     na.rm = TRUE)#finds mean Ct for all samples that had treatement2
      Delta2 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)#finds the absolute value of the difference between Mean_CT_for_treatment1 and Mean_CT_for_treatment2
      print("This is the delta of treatments for your control gene")#prints text
      print(Delta2)#prints the delta for your control genes
      Double_Delta <- abs(Delta1 - Delta2)#finds the absolute value of the difference between Delta1 and Delta2
      print("this is your double delta Ct")#prints text
      print(Double_Delta)#prints Double_Delta
      Fold_induction <- 2^Double_Delta#Finds the fold change
      print("This is the fold change for your gene of interest")#prints text
      print(Fold_induction)#prints Fold_induction
      return(Fold_induction)#returns Fold_induction
    }#lines 110-131 will only be executed if the user only defines the control_gene
    Find_All_Control_Gene_Samples <- grep(control_gene, data$Sample.Name)#finds all samples matching control_gene
    All_Control_Gene_Samples <- data[Find_All_Control_Gene_Samples,]#indexes by samples matching control_gene
    Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)#finds data for a specified time point in All_Control_Gene_Samples
    All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time,]#indexes by time in All_Control_Gene_Samples
    Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)#finds samples that match treatment1 in All_Control_Gene_Samples
    All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1,]#indexes by treatment1 in All_Control_Gene_Samples
    Find_All_treatment2 <-grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)#finds samples that match treatment2 in All_Control_Gene_Samples
    All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2,]#indexes by treatment2 in All_Control_Gene_Samples
    Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), NA == "Undetermined"),
                                   na.rm = TRUE)#finds mean Ct for all samples that had treatement1
    Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), NA == "Undetermined"),
                                   na.rm = TRUE)#finds mean Ct for all samples that had treatement2
    Delta2 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)#finds the absolute value of the difference between Mean_CT_for_treatment1 and Mean_CT_for_treatment2
    print("This is the delta of treatments for your control gene(s)")#prints text
    print(Delta2)#prints the delta for your control gene
    Double_Delta <- abs(Delta1 - Delta2)#finds the absolute value of the difference between Delta1 and Delta2
    print("this is your double delta Ct")#prints text
    print(Double_Delta)#prints Double_Delta
    Fold_induction <- 2^Double_Delta#Finds the fold change
    print("This is the fold induction for your gene of interest")#prints text
    print(Fold_induction)#prints Fold_induction
    return(Fold_induction)#returns Fold_induction
  }
}

combine_folds <- function(Fold, new_file, gene_names){  
  outcon <- file(new_file, open = "w")#creates a new file with a name matching new_file and opens a writing connection to it
  Fold_Induction <- as.character(Fold)#converts Fold from numeric class to the character class
  writeLines(Fold_Induction, outcon)#added fold data to new_file 
  close(outcon)#closes connection to new_file
  read <- read.csv(new_file, header = FALSE)#new_file is read in
  read$Gene_Name <- c(gene_names)#New column with gene_names is added
  read$Fold_change <- read$V1#New column with with the fold change data is added so data is to the right of Gene_Name
  read2 <- read[,-1] #removes the first column
  write.csv(read2,new_file, row.names = FALSE)#creates a csv file in the containing all descriptions associatied with the fold changes
  return(read2)#a data frame returned
}