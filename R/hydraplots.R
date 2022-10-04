# this generates the html output of diagnostics for a Hydra run given a data object and the report file
hydraplots <- function(data, report, outfile = "../outputs/junk.html", plot_comps=TRUE) {
  #data_object = "../test-data/hydraDataList.rda"
  data_object = data
  #repfile <- "../test-data/hydra_sim.rep"
  repfile <- report
  #outfile = "../outputs/junk.html"
  rmarkdown::render(input = "R/hydra-diagnostics.Rmd",
                    output_file = outfile)
}
#example call
hydraplots(data = "../test-data/hydraDataList.rda",
           report = "../test-data/hydra_sim.rep",
           outfile = "../outputs/junk.html")
