\name{deltas}
\alias{deltas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
deltas(data, gene_name, time_point, treatment1, treatment2, control_gene, control_gene2 = FALSE, control_gene3 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{gene_name}{
%%     ~~Describe \code{gene_name} here~~
}
  \item{time_point}{
%%     ~~Describe \code{time_point} here~~
}
  \item{treatment1}{
%%     ~~Describe \code{treatment1} here~~
}
  \item{treatment2}{
%%     ~~Describe \code{treatment2} here~~
}
  \item{control_gene}{
%%     ~~Describe \code{control_gene} here~~
}
  \item{control_gene2}{
%%     ~~Describe \code{control_gene2} here~~
}
  \item{control_gene3}{
%%     ~~Describe \code{control_gene3} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
