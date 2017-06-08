## Echo a message when loading the package
.onAttach <- function(...) {
    # echo output to screen
    packageStartupMessage("##\n## mvMORPH package (1.0.9)")
    packageStartupMessage("## Multivariate evolutionary models")
    packageStartupMessage("##\n## See the tutorials: browseVignettes(\"mvMORPH\")")
    packageStartupMessage("##\n## To cite package 'mvMORPH': citation(\"mvMORPH\")\n##")

}
