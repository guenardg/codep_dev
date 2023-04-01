\encoding{utf8}
\name{Mite}
\docType{data}
\alias{Mite}
\alias{mite.species}
\alias{mite.env}
\alias{mite.geo}
\title{Lac Geai oribatid mites community data}
\description{
  A data set containing abundances of 35 species of Oribatid mites identified 
  in 70 soil cores (sites) extracted from a small area (2.6 m x 10 m) of the 
  peat mat around Lac Geai (+45.9954; -73.9936), 
  a bog lake located on the territory of the Station de Biologie des 
  Laurentides of Université de Montréal, St-Hippolyte, QC, Canada.
}
\usage{Mite}
\format{ Contains three matrices:

  \code{mite.species} The abundance of 35 morpho-species of oribatid
  mites (Acari).

  \code{mite.env} 14 environmental variables (quantitative and binary).

  \code{mite.geo} The coordinates of the sites (cores) in m.
}
\details{
  Values in \code{mite.species} are counts of individuals of each of the
  35 morpho-species obtained from 5-cm diameter cores going from the
  surface of the peat down to a depth of 7 cm. See Bordard & Legendre
  (1994) and reference therein for details about sample treatment and
  species identification.

  \code{mite.env} contains two quantitative variables, namely the
  substratum density (g/L) and water content (percent wet mass over dry
  mass), and 12 binary (or dummy) variables coding for three qualitative variables 	
  (factors). The first seven dummy variables represent
  the substratum composition: \emph{Sphagnum magellacinum} (with
  dominance of \emph{S. rubellum}), \emph{S. rubellum},
  \emph{S. nemorum}, (with a majority of \emph{S. augustifollium}),
  \emph{S. rubellum} + \emph{S. magellicum} (in equal proportions),
  lignous litter, bare peat, and interface between \emph{Sphagnum}
  species. The next three dummy variables represent the presence and
  abundance of the shrubs (\emph{Kalmia polifolia}, \emph{K. angustifolia},
  and \emph{Rhododentron groenlandicum}): none, few, and many. The last
  two dummy variables represent the microtopography of the peat: blanket
  (flat) or hummock (raised).

  Matrix \code{mite.geo} contains the geographic coordinates of the soil cores 
  (sites), in meters, on the peat mat around the lake. The origin of the reference 
  axes is in the lower-left corner of the plot, for an observer looking from the 
  shore towards the open water. The \code{x} coordinate is the distance to the 
  origin along the shore (from left to right) while the \code{y} coordinate is 
  the distance from the shore towards the open water. See Borcard & Legendre 
  (1994, Figs. 1 and 2) for details on the study area.
}
\source{ Daniel Borcard, Département de sciences biologiques, Université
  de Montréal, Montréal, Québec, Canada. }
\references{
  Borcard, D. & Legendre, P. 1994. Environmental control and spatial
  structure in ecological communities: an example using Oribatid mites
  (Acari, Oribatei). Environmental and Ecological Statistics 1: 37-61
}
\seealso{
  Borcard, D.; P. Legendre & P. Drapeau. 1992. Partialling out the
  spatial component of ecological variation. Ecology 73: 1045-1055

  Legendre, P. 2005. Species associations: the Kendall coefficient of
  concordance revisited. Journal of Agricultural, Biological and
  Environmental Statistics 10: 226-245

  Borcard, D.; Gillet, F. & Legendre, P. 2011. Numerical Ecology with
  R. Springer, New-York, NY, USA.
}
\examples{
data(Mite)
summary(mite.species)
summary(mite.env)
summary(mite.geo)
}
\keyword{ Mite }
