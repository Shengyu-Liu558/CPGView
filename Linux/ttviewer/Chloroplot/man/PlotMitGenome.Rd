% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genome.R
\name{PlotMitGenome}
\alias{PlotMitGenome}
\title{Generate Mitochondrion Genome Plot}
\usage{
PlotMitGenome(
  plot.tables,
  save = TRUE,
  file.type = "pdf",
  text.size = 1,
  height = c(0.1, 0.2, 0.07),
  file.name = NULL,
  gc.per.gene = TRUE,
  pseudo = TRUE,
  legend = TRUE,
  genome.length = TRUE,
  total.gc = TRUE,
  gene.no = TRUE,
  rrn.no = TRUE,
  trn.no = TRUE,
  organelle_type = TRUE,
  background = "grey90",
  gc.color = "grey30",
  gc.background = "grey70",
  info.background = "black",
  nad.color = "#2A6332",
  sdh.color = "#4C8805",
  cob.color = "#7F992C",
  cox.color = "#FEEE50",
  atp.color = "#9FBB3D",
  ccmF.color = "#4D9E3F",
  mmt.color = "#AE2D29",
  rps.color = "#D6AD7C",
  rpl.color = "#9C7A4B",
  mat.color = "#D9662D",
  orf.color = "#71B8A9",
  trn.color = "#172C7F",
  rrn.color = "#D1382A",
  other_gene.color = "#7D7D7D",
  show.genes = c("nad|nd", "sdh", "cob", "cox|cytb", "atp", "ccmF", "mtt", "rps",
    "rpl", "mat", "orf", "trn", "rrn", "OTHER"),
  cu.bias = TRUE,
  customize.ring1 = NULL,
  customize.ring1.color = "grey30",
  customize.ring1.type = "line",
  customize.ring2.color = "grey30",
  customize.ring2 = NULL,
  customize.ring2.type = "line",
  customize.ring3 = NULL,
  customize.ring3.color = "grey30"
)
}
\arguments{
\item{plot.tables}{a list contains information of IR region, gene, and gc
count of the genome. It can be generated by function \code{\link{PlotTab}}.}

\item{save}{A logical value. If it is \code{TRUE}, the plot will be saved in
work directory.}

\item{file.type}{A charactor. It indicates the format of the file in which
the plot will be saved. Options are: "pdf", "png", "jpeg","bmp", "tiff".
By default, it is set as "pdf".}

\item{text.size}{A numeric value. It indicates the size of all texts in the
plot. For exmple, \code{text.size = 1.5} means all the texts in the plot will
be enlarged as 1.5 times of their original size.}

\item{height}{A vector of numeric value. The elements of it indicat the
height of gene layer, GC count layer and IR region layer, respectively.
Default setting is "0.1, 0.2, 0.07". The taltal circle plot region always
has a radius of 1, so a height of 0.1 means 10\% of the circle radius.}

\item{file.name}{A charactor. It indicates the name of the file in which
the plot will be saved. By default, it is set as the specie's name.}

\item{gc.per.gene}{A logical value. If it is \code{TRUE}, the GC content of
each gene will be show by a darker part in gene rectangles.}

\item{pseudo}{A logical value. If it is \code{TRUE}, the pseudo genes (if
there is some in the species) will be marked with a "*" at the end of the
labels.}

\item{legend}{A logical value. If it is \code{TRUE}, the legend for gene
colors will be shown.}

\item{genome.length}{A logical value. If it is \code{TRUE}, the length of
genome will be shown in the center of the plot.}

\item{total.gc}{A logical value. If it is \code{TRUE}, the GC content of
whole genome will be shown in the center of the plot.}

\item{gene.no}{A logical value. If it is \code{TRUE}, the number of genes in
the genome will be shown in the center of the plot.}

\item{rrn.no}{A logical value. If it is \code{TRUE}, the number of rRNAs in
the genome will be shown in the center of the plot.}

\item{trn.no}{A logical value. If it is \code{TRUE}, the number of tRNAs in
the genome will be shown in the center of the plot.}

\item{background}{An R color object. It indicates the color for the background of
entire plot area.}

\item{gc.color}{An R color object. It indicates the color for the lines in gc count
plot.}

\item{gc.background}{An R color object. It indicates the color for the background
of gc count plot area.}

\item{info.background}{An R color object. It indicates the color for the background
of central area where species' information was shown.}

\item{nad.color}{An R color string. It indicates the color for genes of
complex I (NADH dehydrogenase)}

\item{sdh.color}{An R color string. It indicates the color for genes of
complex II (succinate dehydrogenase)}

\item{cob.color}{An R color string. It indicates the color for genes of
complex III (ubichinol cytochrome reductase)}

\item{cox.color}{An R color string. It indicates the color for genes of
complex IV (cytochrome c oxidase)}

\item{atp.color}{An R color string. It indicates the color for genes of
ATP synthase}

\item{ccmF.color}{An R color string. It indicates the color for genes of
cytochrome c biogenesi}

\item{mmt.color}{An R color string. It indicates the color for genes of
transport membrane protein}

\item{rps.color}{An R color string. It indicates the color for genes of
ribosomal proteins (SSU)}

\item{rpl.color}{An R color string. It indicates the color for genes of
ribosomal proteins (LSU)}

\item{mat.color}{An R color string. It indicates the color for genes of
maturases}

\item{orf.color}{An R color string. It indicates the color for ORFs}

\item{trn.color}{An R color object. It indicates the color for genes of
transfer RNA}

\item{rrn.color}{An R color object. It indicates the color for genes of
ribosomal RNA}

\item{other_gene.color}{An R color object. It indicates the color for other
genes}

\item{cu.bias}{A logical value. If it is \code{TRUE}, the condon usage bias
metric for each gene will be shown in the labels. The metric ???Measure
Independent of Length and Composition (MILC)??? was used for evaluate the bias.
You can find more details in the reference paper.}

\item{customize.ring1}{A data frame. It must contain 2 columns:
\itemize{
  \item \strong{position}: 1-base genomic coordinate for the features.
  \item \strong{value}: the values for the features.
}}

\item{customize.ring1.color}{An R color object. It indicates the color for
the plots in customized ring 1.}

\item{customize.ring1.type}{A character. It indicate the plot type in
customize.ring1. Available values are "line", "line + filling", "line + dot",
"line + dot + filling", "step line", "step line + filling", "vertical line"}

\item{customize.ring2.color}{An R color object. It indicates the color for
the plots in customized ring 2.}

\item{customize.ring2}{A data frame. It must contain 2 columns:
\itemize{
  \item \strong{position}: 1-base genome coordinate for the features.
  \item \strong{value}: the values for the features.
}}

\item{customize.ring2.type}{A character. It indicate the plot type in
customize.ring2. Available values are "line", "line + filling", "line + dot",
"line + dot + filling", "step line", "step line + filling", "vertical line"}

\item{customize.ring3}{A data frame. It must contain 2 columns:
\itemize{
  \item \strong{start}: 1-base genome coordinate for the start point of the
  features.
  \item \strong{end}: 1-base genome coordinate for the end point of the
  features.
  \item \strong{value}: the values for the features.
}}

\item{customize.ring3.color}{An R color object. It indicates the color for
the plots in customized ring 3.}

\item{organelle}{A logical value. If it is \code{TRUE}, the organelle type of
the genome will be shown in the center of the plot.}

\item{show.gene}{A vector of characters. It indicates which classes of genes
will be shown on the plot. Avaliable values are "nad|nd","sdh","cob",
"cox|cytb","atp", "ccmF","mtt","rps","rpl", "mat","orf","trn","rrn", "OTHER"}
}
\value{
A plot for chloroplast genome.
}
\description{
Generate Mitochondrion Genome Plot
}
\references{
Supek, Fran, and Kristian Vlahovicek. ???Comparison of codon usage
measures and their applicability in prediction of microbial gene
expressivity.??? BMC bioinformatics vol. 6 182. 19 Jul.
2005, doi:10.1186/1471-2105-6-182
}
