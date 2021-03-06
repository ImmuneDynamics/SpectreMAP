% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.spatial.files.R
\name{read.spatial.files}
\alias{read.spatial.files}
\title{Read TIFF files into R and and create spatial data object}
\usage{
read.spatial.files()
}
\arguments{
\item{rois}{vector of ROI names (directory names)}

\item{roi.loc}{working directory for ROIs}

\item{multi.tiff}{Default FALSE}

\item{correct.extent}{Default TRUE}

\item{flip.y}{Default TRUE}

\item{value.modifier}{Default 65535}

\item{ext}{Default = ".tiff"}
}
\value{
Returns a spatial data object.
}
\description{
Read TIFF files into R and and create spatial data object
}
\examples{

}
