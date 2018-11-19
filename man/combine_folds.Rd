\name{combine_folds}
\alias{combine_folds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine Results from Deltas Function into One CSV File
}
\description{
Output taken from the deltas function and adds it to a new CSV file that is named by the user. The user will can also to add labels for each of the fold changes.
}
\usage{
combine_folds(Fold, new_file, gene_names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fold}{
Fold change results from the Deltas function. Users can use a single result or a vector of results in this argument.
}
  \item{new_file}{
Name of the new CSV file. Must be in quotes and have the .csv extension.
}
  \item{gene_names}{
A character string containing descriptive information about the fold changes.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Fold, new_file, gene_names) 
{
    outcon <- file(new_file, open = "w")
    Fold_Induction <- as.character(Fold)
    writeLines(Fold_Induction, outcon)
    close(outcon)
    read <- read.csv(new_file, header = FALSE)
    read$Gene_Name <- c(gene_names)
    read$Fold_change <- read$V1
    read2 <- read[, -1]
    write.csv(read2, new_file, row.names = FALSE)
    return(read2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
