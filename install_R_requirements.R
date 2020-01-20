repos="http://cran.r-project.org"
if (!suppressWarnings(require("dplyr"))) {
    install.packages("dplyr", repos=repos)
}
if (!suppressWarnings(require("doParallel"))) {
    install.packages("doParallel", repos=repos)
}
if (!suppressWarnings(require("devtools"))) {
    install.packages("devtools", repos=repos)
}
if (!suppressWarnings(require("ccube"))) {
	devtools::install_github('keyuan/ccube')
}
