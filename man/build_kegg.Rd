% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_kegg.R
\name{build_kegg}
\alias{build_kegg}
\title{Constructs a metabolic network from KEGG data for specified species.}
\usage{
build_kegg(
  Merged_data,
  info,
  species = "rno",
  gene_id,
  met_id,
  size = "MetMN_Expand",
  write_output = FALSE,
  output_dir = "outputs"
)
}
\arguments{
\item{Merged_data}{Input dataset containing merged metabolic and transcriptomic data}

\item{info}{Metadata for column standardization}

\item{species}{Species identifier ("rno", "mmu", or "hsa")}

\item{gene_id}{Column name for gene identifiers}

\item{met_id}{Column name for metabolite identifiers}

\item{size}{Network size type (default: "MetMN_Expand")}

\item{write_output}{Should output files be written? (default: FALSE)}

\item{output_dir}{Output directory path (default: "outputs")}
}
\value{
List containing node attributes, edge attributes, and solid edges
}
\description{
Constructs a metabolic network from KEGG data for specified species.
}
