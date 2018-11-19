\name{Read_Combine_Quality}
\alias{Read_Combine_Quality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine RT-qPCR Data in CSV files.
}
\description{
Function will combine data from up to three CSV files, remove data that resulted in an undetermined Ct, and warn users if any of the data was flagged as a non-amplified well (NAW). 
}
\usage{
Read_Combine_Quality(File1, File2 = FALSE, File3 = FALSE, New_name.csv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{File1}{
A CSV file containing RT-qPCR results generated from SDS 2.3 software. 
}
  \item{File2}{
An additional CSV file containing RT-qPCR results generated from SDS 2.3 software that the user would like to combine with File1. Default value is FALSE.
}
  \item{File3}{
An additional CSV file containing RT-qPCR results generated from SDS 2.3 software that the user would like to combine with File1 and File2. Default value is FALSE. 
}
  \item{New_name.csv}{
The name the user would like the new file to be called. This should be in quoutes and have the .csv file extension.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns a data frame.
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
function (File1, File2 = FALSE, File3 = FALSE, New_name.csv) 
{
    if (File3 != FALSE) {
        Data1 <- read.csv(File1)
        Data2 <- read.csv(File2)
        Data3 <- read.csv(File3)
        All_data <- rbind(Data1, Data2, Data3)
        Remove_undetermined <- All_data[which(All_data$Ct != 
            "Undetermined"), ]
        if (any(Remove_undetermined$NAW == "TRUE")) {
            warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
            write.csv(Remove_undetermined, New_name.csv)
            return(Remove_undetermined)
        }
    }
    else {
        if (File2 != FALSE) {
            Data1 <- read.csv(File1)
            Data2 <- read.csv(File2)
            All_data <- rbind(Data1, Data2)
            Remove_undetermined <- All_data[which(All_data$Ct != 
                "Undetermined"), ]
            if (any(Remove_undetermined$NAW == "TRUE")) {
                warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
                write.csv(Remove_undetermined, New_name.csv)
                return(Remove_undetermined)
            }
        }
    }
    Data <- read.csv(File1)
    Remove_undetermined <- Data[which(Data$Ct != "Undetermined"), 
        ]
    if (any(Remove_undetermined$NAW == "TRUE")) {
        warning("Non-Amplified Well (NAW) was TRUE for values. Consider removing them.")
        write.csv(Remove_undetermined, New_name.csv)
        return(Remove_undetermined)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
