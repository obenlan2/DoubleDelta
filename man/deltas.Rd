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
#mypath <- system.file("extdata", "example1.csv", package = "DoubleDelta")

#deltas(mypath, "Gene1", "12HAT", "UT", "ST", "BTUB", "GTPB", "GAPDH")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
