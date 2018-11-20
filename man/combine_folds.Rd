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
data <- system.file("extdata", "example1", package = "DoubleDelta_0.0.tar.gz")
mydata1 <- read.table(data)
data3 <- system.file("extdata", "example3", package = "DoubleDelta_0.0.tar.gz")
mydata2 <- read.table(data2)

a <- deltas(mydata1, "Gene1", "6HAT", "UT", "ST", "BTUB", "GTPB", "GAPDH")
b <- deltas(mypath2, "Gene1", "12HAT", "UT", "ST", "BTUB", "GTPB", "GAPDH")

combine_folds(c(a,b), "All_Data.csv", c("Gene1 6HAT", "Gene2 12HAT"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
