\name{Read_Combine_Quality}
\alias{Read_Combine_Quality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Read_Combine_Quality(File1, File2 = FALSE, File3 = FALSE, New_name.csv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{File1}{
%%     ~~Describe \code{File1} here~~
}
  \item{File2}{
%%     ~~Describe \code{File2} here~~
}
  \item{File3}{
%%     ~~Describe \code{File3} here~~
}
  \item{New_name.csv}{
%%     ~~Describe \code{New_name.csv} here~~
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
