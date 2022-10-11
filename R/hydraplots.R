# this generates the html output of diagnostics for a Hydra run given a data object and the report file
hydraplots <- function(data, report, outfile = "../outputs/junk.html", plot_comps=TRUE, skill=FALSE, truth=NULL, do_osa = FALSE) {
  #data_object = "../test-data/hydraDataList.rda"
  data_object = data
  #repfile <- "../test-data/hydra_sim.rep"
  repfile <- report

  #if not doing composition plots then automagically not doing OSA either
  if (plot_comps == FALSE) do_osa <- FALSE

  if(skill){
    truth <- truth
  }
  #outfile = "../outputs/junk.html"
  rmarkdown::render(input = "R/hydra-diagnostics.Rmd",
                    output_file = outfile)
}
#example call
hydraplots(data = "../test-data/hydraDataList.rda",
           report = "../test-data/hydra_sim.rep",
           outfile = "../outputs/junk.html")
