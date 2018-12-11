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
\value{
This function returns a data frame. Column and row names will be identical to the inputted files and any rows that had "Undetermined" results for the Ct value will be removed. Empty cells will have "NA".
}

\seealso{
\code{\link{deltas}}
\code{\link{combine_folds}}
}
\examples{
data <- system.file("extdata", "example1.csv", package = "DoubleDelta")
data2 <- system.file("extdata", "example2.csv", package = "DoubleDelta")
data3 <- system.file("extdata", "example3.csv", package = "DoubleDelta")

#Making Temporary Files Located in Temp Directory
outfile1 <- tempfile(fileext = "NewFile.csv")
outfile2 <- tempfile(fileext = "NewFile2.csv")
outfile3 <- tempfile(fileext = "NewFile3.csv")


#Run the following code to put "NewFile.csv", "NewFile2.csv", or "NewFile3.csv" in the Temp Directory.

##combining three files
Read_Combine_Quality(File1 = data, File2 = data2, File3 = data3, outfile1)
##combining two files
Read_Combine_Quality(File1 = data, File2 = data2, New_name.csv = outfile2)
##inputing one file
Read_Combine_Quality(File1 = data, New_name.csv = outfile3)

\dontrun{
#Run the following code to view "NewFile.csv", "NewFile2.csv", or "NewFile3.csv" in your designated directory.

##combining three files
Read_Combine_Quality(File1 = data, File2 = data2, File3 = data3, "NewFile.csv")
##combining two files
Read_Combine_Quality(File1 = data, File2 = data2, New_name.csv = "NewFile2.csv")
##inputing one file
Read_Combine_Quality(File1 = data, New_name.csv = "NewFile3.csv")
}

}


