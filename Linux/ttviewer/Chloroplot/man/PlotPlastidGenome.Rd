% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genome.R
\name{PlotPlastidGenome}
\alias{PlotPlastidGenome}
\title{Generate Plasit Genome Plot}
\usage{
PlotPlastidGenome(
  plot.tables,
  save = TRUE,
  file.type = "pdf",
  text.size = 1,
  height = c(0.1, 0.2, 0.07),
  file.name = NULL,
  shadow = TRUE,
  ir.gc = TRUE,
  gc.per.gene = TRUE,
  pseudo = TRUE,
  legend = TRUE,
  ssc.converse = FALSE,
  lsc.converse = FALSE,
  ira.converse = FALSE,
  irb.converse = FALSE,
  genome.length = TRUE,
  total.gc = TRUE,
  show.indel = TRUE,
  gene.no = TRUE,
  rrn.no = TRUE,
  trn.no = TRUE,
  organelle_type = TRUE,
  background = "grey90",
  gc.color = "grey30",
  gc.background = "grey70",
  info.background = "black",
  ir.color = "#2F3941",
  ssc.color = "#82B6E2",
  lsc.color = "#299E96",
  shadow.color = "#0000FF20",
  psa.color = "#2A6332",
  psb.color = "#4C8805",
  pet.color = "#7F992C",
  atp.color = "#9FBB3D",
  ndh.color = "#FEEE50",
  rbc.color = "#4D9E3F",
  rpo.color = "#AE2D29",
  rsp.color = "#D6AD7C",
  rpl.color = "#9C7A4B",
  clp_mat_inf.color = "#D9662D",
  ycf.color = "#71B8A9",
  trn.color = "#172C7F",
  rrn.color = "#D1382A",
  other_gene.color = "#7D7D7D",
  show.genes = c("psa", "psb", "pet", "atp", "ndh", "rbc", "rpo", "rps", "rpl",
    "clp|mat|inf", "ycf", "trn", "rrn", "OTHER"),
  gene_axis_ir.color = ir.color,
  gene_axis_ssc.color = ssc.color,
  gene_axis_lsc.color = lsc.color,
  cu.bias = TRUE,
  customize.ring1 = NULL,
  customize.ring1.type = "line",
  customize.ring1.color = gc.color,
  customize.ring2 = NULL,
  customize.ring2.type = "line",
  customize.ring2.color = gc.color,
  customize.ring3 = NULL,
  customize.ring3.color = gc.color
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

\item{shadow}{A logical value. If it is \code{TRUE}, the section of IR
regions will be highlighted by shadows.}

\item{ir.gc}{A logical value. If it is \code{TRUE}, the gc content for
IR, LSC, and SSC regions will be plot with deep colors.}

\item{gc.per.gene}{A logical value. If it is \code{TRUE}, the GC content of
each gene will be show by a darker part in gene rectangles.}

\item{pseudo}{A logical value. If it is \code{TRUE}, the pseudo genes (if
there is some in the species) will be marked with a "*" at the end of the
labels.}

\item{legend}{A logical value. If it is \code{TRUE}, the legend for gene
colors will be shown.}

\item{ssc.converse}{A logical value. If it is \code{TRUE}, the SSC region
will be converted to its reverse complementary version}

\item{lsc.converse}{A logical value. If it is \code{TRUE}, the LSC region
will be converted to its reverse complementary version}

\item{ira.converse}{A logical value. If it is \code{TRUE}, the IRA region
will be converted to its reverse complementary version}

\item{irb.converse}{A logical value. If it is \code{TRUE}, the IRB region
will be converted to its reverse complementary version}

\item{genome.length}{A logical value. If it is \code{TRUE}, the length of
genome will be shown in the center of the plot.}

\item{total.gc}{A logical value. If it is \code{TRUE}, the GC content of
whole genome will be shown in the center of the plot.}

\item{show.indel}{A logical value. If it is \code{TRUE}, the SNP, insertion
and deletion areas in IR regions will be highlighted on the plot.}

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

\item{ir.color}{An R color object. It indicates the color for the background
of IR sectors.}

\item{ssc.color}{An R color object. It indicates the color for the background
of SSC sectors.}

\item{lsc.color}{An R color object. It indicates the color for the background
of LSC sectors.}

\item{shadow.color}{An R color object. It indicates the color for the shadow
casted from IR sectors.}

\item{psa.color}{An R color object. It indicates the color for genes of
photosystem I}

\item{psb.color}{An R color object. It indicates the color for genes of
photosystem II}

\item{pet.color}{An R color object. It indicates the color for genes of
cytochrome b/f complex}

\item{atp.color}{An R color object. It indicates the color for genes of
ATP synthesis}

\item{ndh.color}{An R color object. It indicates the color for genes of
NADH dehydrogenase}

\item{rbc.color}{An R color object. It indicates the color for genes of
RubisCO larg subunit}

\item{rpo.color}{An R color object. It indicates the color for genes of
RNA polymerase}

\item{rsp.color}{An R color object. It indicates the color for genes of
small ribosomal protein}

\item{rpl.color}{An R color object. It indicates the color for genes of
large ribosomal protein}

\item{clp_mat_inf.color}{An R color object. It indicates the color for genes
of clpP, matK, or infA}

\item{ycf.color}{An R color object. It indicates the color for genes of
hypothetical reading frame}

\item{trn.color}{An R color object. It indicates the color for genes of
transfer RNA}

\item{rrn.color}{An R color object. It indicates the color for genes of
ribosomal RNA}

\item{other_gene.color}{An R color object. It indicates the color for other
genes}

\item{gene_axis_ir.color}{An R color object. It indicates the color for other
genes.}

\item{gene_axis_ssc.color}{An R color object. It indicates the color for other
genes.}

\item{gene_axis_lsc.color}{An R color object. It indicates the color for other
genes.}

\item{cu.bias}{A logical value. If it is \code{TRUE}, the condon usage bias
metric for each gene will be shown in the labels. The metric ???Measure
Independent of Length and Composition (MILC)??? was used for evaluate the bias.
You can find more details in the reference paper.}

\item{customize.ring1}{A data frame. It must contain 2 columns:
\itemize{
  \item \strong{position}: 1-base genomic coordinate for the features.
  \item \strong{value}: the values for the features.
}}

\item{customize.ring1.type}{A charactor. It indicate the plot type in
customize.ring1. Avaliable values are "line", "line + filling", "line + dot",
"line + dot + filling", "step line", "step line + filling", "vertical line"}

\item{customize.ring1.color}{An R color object. It indicates the color for
the plots in customized ring 1.}

\item{customize.ring2}{A data frame. It must contain 2 columns:
\itemize{
  \item \strong{position}: 1-base genomic coordinate for the features.
  \item \strong{value}: the values for the features.
}}

\item{customize.ring2.type}{A charactor. It indicate the plot type in
customize.ring2. Avaliable values are "line", "line + filling", "line + dot",
"line + dot + filling", "step line", "step line + filling", "vertical line"}

\item{customize.ring2.color}{An R color object. It indicates the color for
the plots in customized ring 2.}

\item{customize.ring3}{A data frame. It must contain 2 columns:
\itemize{
  \item \strong{start}: 1-base genomic coordinate for the start point of the
  features.
  \item \strong{end}: 1-base genomic coordinate for the end point of the
  features.
  \item \strong{value}: the values for the features.
}}

\item{customize.ring3.color}{An R color object. It indicates the color for
the plots in customized ring 3.}

\item{organelle}{A logical value. If it is \code{TRUE}, the organelle type of
the genome will be shown in the center of the plot.}

\item{show.gene}{A vector of characters. It indicates which classes of genes
will be shown on the plot. A valiable values are "psa","psb","pet","atp",
"ndh","rbc", "rpo","rps","rpl", "clp|mat|inf","ycf", "trn","rrn", "OTHER"}
}
\value{
A plot for chloroplast genome.
}
\description{
Generate Plasit Genome Plot
}
\references{
Supek, Fran, and Kristian Vlahovicek. ???Comparison of codon usage
measures and their applicability in prediction of microbial gene
expressivity.??? BMC bioinformatics vol. 6 182. 19 Jul.
2005, doi:10.1186/1471-2105-6-182
}
