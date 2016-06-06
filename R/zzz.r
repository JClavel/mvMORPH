## Echo a message when loading the package
.onAttach <- function(...) {
    # echo output to screen
    packageStartupMessage("##\n## mvMORPH package (1.0.7)")
    packageStartupMessage("## Multivariate evolutionary models")
    packageStartupMessage("##\n## See the tutorials: browseVignettes(\"mvMORPH\")")
    packageStartupMessage("##\n## Clavel et al. (2015)")
    packageStartupMessage("## (Methods in Ecology and Evolution. 6(11):1311-1319)\n##")

}
