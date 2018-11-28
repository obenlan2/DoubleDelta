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

\value{
This function returns a data frame.
}

\seealso{
\code{\link{deltas}}
\code{\link{Read_Combine_Quality}}
}
\examples{
mypath1 <- system.file("extdata", "example1.csv", package = "DoubleDelta")
mypath3 <- system.file("extdata", "example3.csv", package = "DoubleDelta")
data1 <- read.csv(mypath1)
data3 <- read.csv(mypath3)

a <- deltas(data1, "Gene1", "6HAT", "UT", "ST", "BTUB", "GTPB", "GAPDH")
b <- deltas(data3, "Gene1", "12HAT", "UT", "ST", "BTUB", "GTPB", "GAPDH")

combine_folds(c(a,b), "All_Data.csv", c("Gene1 6HAT", "Gene2 12HAT"))
}

