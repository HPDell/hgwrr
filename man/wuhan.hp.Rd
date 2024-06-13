\name{wuhan.hp}
\alias{wuhan.hp}
\docType{data}
\title{Wuhan Second-hand House Price and POI Data (DataFrame)}
\description{
A data set of second-hand house price in Wuhan, China collected in 2018.
}
\usage{data(multisampling)}
\format{
A list of two items called "data" and "coords".
Item "data" is a data frame with 13862 second-hand properties at 779 neighbourhoods
and the following 22 variables.
\describe{
    \item{\code{Price}}{House price per square metre.}
    \item{\code{Floor.High}}{1 if a property is on a high floor, otherwise 0.}
    \item{\code{Floor.Low}}{1 if a property is on a low floor, otherwise 0.}
    \item{\code{Decoration.Fine}}{1 if a property is well decorated, otherwise 0.}
    \item{\code{PlateTower}}{1 if a property is of the plate-tower type, otherwise 0.}
    \item{\code{Steel}}{1 if a property is of 'steel' structure, otherwise 0.}
    \item{\code{BuildingArea}}{Building area in square metres.}
    \item{\code{Fee}}{Management fee per square meter per month.}
    \item{\code{d.Commercial}}{Distance to the nearest commercial area.}
    \item{\code{d.Greenland}}{Distance to the nearest green land.}
    \item{\code{d.Water}}{Distance to the nearest river or lake.}
    \item{\code{d.University}}{Distance to the nearest university.}
    \item{\code{d.HighSchool}}{Distance to the nearest high school.}
    \item{\code{d.MiddleSchool}}{Distance to the nearest middle school.}
    \item{\code{d.PrimarySchool}}{Distance to the nearest primary school.}
    \item{\code{d.Kindergarten}}{Distance to the nearest kindergarten.}
    \item{\code{d.SubwayStation}}{Distance to the nearest subway station.}
    \item{\code{d.Supermarket}}{Distance to the nearest supermarket.}
    \item{\code{d.ShoppingMall}}{Distance to the nearest shopping mall.}
    \item{\code{lon}}{Longitude coordinates (Projected CRS: EPSG 3857).}
    \item{\code{lat}}{Latitude coordinates (Projected CRS: EPSE 3857).}
    \item{\code{group}}{Group id of each sample.}
}

The following variables are group level:

- \code{Fee}
- \code{d.Commercial}
- \code{d.Greenland}
- \code{d.Water}
- \code{d.University}
- \code{d.HighSchool}
- \code{d.MiddleSchool}
- \code{d.PrimarySchool}
- \code{d.Kindergarten}
- \code{d.SubwayStation}
- \code{d.Supermarket}
- \code{d.ShoppingMall}

The follwoing variables are sample level:

- \code{Price}
- \code{Floor.High}
- \code{Floor.Low}
- \code{Decoration.Fine}
- \code{PlateTower}
- \code{Steel}
- \code{BuildingArea}

Item "coords" is a 779-by-2 matrix of coordinates of all neighbourhoods.
}
\author{Yigong Hu \email{yigong.hu@bristol.ac.uk}}
\examples{
\dontrun{
data(wuhan.hp)
hgwr(
  formula = Price ~ L(d.Water + d.Commercial + d.PrimarySchool +
            d.Kindergarten + Fee) + BuildingArea + (Floor.High | group),
  data = wuhan.hp$data,
  coords = wuhan.hp$coords, bw = 50, kernel = "bisquared")
}
}
