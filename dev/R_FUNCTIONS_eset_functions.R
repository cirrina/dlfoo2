
# You can learn more about package authoring with RStudio at: http://r-pkgs.had.co.nz/

# devtools::use_vignette('my-vignette')

# Imports: methods, Biobase, AnnotationDbi, org.Hs.eg.db


# .onLoad <- function(libname, pkgname) { op <- options() op.devtools <- list( devtools.path = '~/R/PACKAGES/dlfoo_v1.2/', devtools.install.args = '', devtools.name = 'Your
# name goes here', devtools.desc.author = ''First Last <first.last@example.com> [aut, cre]'', devtools.desc.license = 'What license is it under?', devtools.desc.suggests =
# NULL, devtools.desc = list() ) toset <- !(names(op.devtools) %in% names(op)) if(any(toset)) options(op.devtools[toset]) packageStartupMessage('Welcome to dlfoo package -
# my functions loaded') invisible() }

