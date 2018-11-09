\name{combine_folds}
\alias{combine_folds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
combine_folds(f1, new_file, gene_names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f1}{
%%     ~~Describe \code{f1} here~~
}
  \item{new_file}{
%%     ~~Describe \code{new_file} here~~
}
  \item{gene_names}{
%%     ~~Describe \code{gene_names} here~~
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
function (f1, new_file, gene_names) 
{
    outcon <- file(new_file, open = "w")
    Fold_Induction <- as.character(f1)
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
