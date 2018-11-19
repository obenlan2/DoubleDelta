\name{deltas}
\alias{deltas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Fold Changes for RT-qPCR Data
}
\description{
The Ct values from a specified gene of interest (GOI) for two specified treatments will be extracted, the mean Ct for both treatments will be calculated, and the delta Ct value between the treatments will be calculated. The same will be done for a specified reference gene(s). This function can normalize the GOI delta Ct with up to three reference genes. The delta Ct of the reference gene(s) will be used to normalize the delta Ct of the GOI and generate the double delta Ct. The double delta Ct will be used to calculate the fold change. Delta and double delta Cts for the GOI and the reference gene(s) will be printed along with the fold change which will be returned at the end of the function. 
}
\usage{
deltas(data, gene_name, time_point, treatment1, treatment2, control_gene, control_gene2 = FALSE, control_gene3 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A CSV file containing RT-qPCR results generated from SDS 2.3 software or a data frame generated from the Read_Combine_Quality function.
}
  \item{gene_name}{
The name of the user's GOI. Needs to be in quotes.
}
  \item{time_point}{
The name of the user's time-point of interest. Needs to be in quotes.
}
  \item{treatment1}{
The name of the user's treatment. Needs to be in quotes.
}
  \item{treatment2}{
The name of the user's treatment. Needs to be in quotes.
}
  \item{control_gene}{
The name of one control gene used in the experiment. Needs to be in quotes.
}
  \item{control_gene2}{
The name of the second control gene used in the experiment. Needs to be in quotes.
}
  \item{control_gene3}{
The name of the third control gene used in the experiment. Needs to be in quotes.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Delta and double delta Cts for the GOI and the reference gene(s) will be printed along with the fold change in order to view results of the calculations. However, the fold change (which is in the numeric class) is the only value returned with the return function and qwhen pointing the function to a name, the fold change will be the only value that will be associated with the given name.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
mypath <- system.file("extdata", "example1", package = "DoubleDelta_0.0.tar.gz")


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, gene_name, time_point, treatment1, treatment2, 
    control_gene, control_gene2 = FALSE, control_gene3 = FALSE) 
{
    Find_All_Gene_Samples <- grep(gene_name, data$Sample.Name)
    All_Gene_Samples <- data[Find_All_Gene_Samples, ]
    Find_Specific_Time <- grep(time_point, All_Gene_Samples$Sample.Name)
    All_Gene_Samples_At_Specific_Time <- All_Gene_Samples[Find_Specific_Time, 
        ]
    Find_All_treatment1 <- grep(treatment1, All_Gene_Samples_At_Specific_Time$Sample.Name)
    All_treatment1_Samples <- All_Gene_Samples_At_Specific_Time[Find_All_treatment1, 
        ]
    Find_All_treatment2 <- grep(treatment2, All_Gene_Samples_At_Specific_Time$Sample.Name)
    All_treatment2_Samples <- All_Gene_Samples_At_Specific_Time[Find_All_treatment2, 
        ]
    Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), 
        NA == "Undetermined"), na.rm = TRUE)
    Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), 
        NA == "Undetermined"), na.rm = TRUE)
    Delta1 <- abs(Mean_CT_for_treatment2 - Mean_CT_for_treatment1)
    print("This is the delta of treatments for your gene of interest")
    print(Delta1)
    if (control_gene3 != FALSE) {
        Find_All_Control_Gene_Samples1 <- grep(control_gene, 
            data$Sample.Name)
        All_Control_Gene_Samples1 <- data[Find_All_Control_Gene_Samples1, 
            ]
        Find_All_Control_Gene_Samples2 <- grep(control_gene2, 
            data$Sample.Name)
        All_Control_Gene_Samples2 <- data[Find_All_Control_Gene_Samples2, 
            ]
        Find_All_Control_Gene_Samples3 <- grep(control_gene3, 
            data$Sample.Name)
        All_Control_Gene_Samples3 <- data[Find_All_Control_Gene_Samples3, 
            ]
        All_Control_Gene_Samples <- rbind(All_Control_Gene_Samples1, 
            All_Control_Gene_Samples2, All_Control_Gene_Samples3)
        Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)
        All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time, 
            ]
        Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
        All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1, 
            ]
        Find_All_treatment2 <- grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
        All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2, 
            ]
        Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), 
            NA == "Undetermined"), na.rm = TRUE)
        Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), 
            NA == "Undetermined"), na.rm = TRUE)
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
    }
    else {
        if (control_gene2 != FALSE) {
            Find_All_Control_Gene_Samples1 <- grep(control_gene, 
                data$Sample.Name)
            All_Control_Gene_Samples1 <- data[Find_All_Control_Gene_Samples1, 
                ]
            Find_All_Control_Gene_Samples2 <- grep(control_gene2, 
                data$Sample.Name)
            All_Control_Gene_Samples2 <- data[Find_All_Control_Gene_Samples2, 
                ]
            All_Control_Gene_Samples <- rbind(All_Control_Gene_Samples1, 
                All_Control_Gene_Samples2)
            Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)
            All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time, 
                ]
            Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
            All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1, 
                ]
            Find_All_treatment2 <- grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
            All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2, 
                ]
            Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), 
                NA == "Undetermined"), na.rm = TRUE)
            Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), 
                NA == "Undetermined"), na.rm = TRUE)
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
        }
        Find_All_Control_Gene_Samples <- grep(control_gene, data$Sample.Name)
        All_Control_Gene_Samples <- data[Find_All_Control_Gene_Samples, 
            ]
        Find_Specific_Time <- grep(time_point, All_Control_Gene_Samples$Sample.Name)
        All_Control_Gene_Samples_At_Specific_Time <- All_Control_Gene_Samples[Find_Specific_Time, 
            ]
        Find_All_treatment1 <- grep(treatment1, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
        All_treatment1_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment1, 
            ]
        Find_All_treatment2 <- grep(treatment2, All_Control_Gene_Samples_At_Specific_Time$Sample.Name)
        All_treatment2_Samples <- All_Control_Gene_Samples_At_Specific_Time[Find_All_treatment2, 
            ]
        Mean_CT_for_treatment1 <- mean(as.numeric(as.character(All_treatment1_Samples$Ct), 
            NA == "Undetermined"), na.rm = TRUE)
        Mean_CT_for_treatment2 <- mean(as.numeric(as.character(All_treatment2_Samples$Ct), 
            NA == "Undetermined"), na.rm = TRUE)
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
