# this generates the html output of diagnostics for a Hydra run given a data object and the report file
hydra_comparisons <- function(data, reports, outfile = "../outputs/junk2.html") {
  #data_object = "../test-data/hydraDataList.rda"
  data_object = data
  #repfile <- "../test-data/hydra_sim.rep"
  repfiles <- reports
  #outfile = "../outputs/junk.html"
  rmarkdown::render(input = "R/comparison-plots.Rmd",
                    output_file = outfile)
}
#example call
hydra_comparisons(data = "../test-data/hydraDataList.rda",
                  reports = list("../test-data/hydra_sim.rep","../../hydra_sim/results/005.rep"),
                  outfile = "../outputs/comparison-junk.html")
