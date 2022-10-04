## =============================================================================
##                  SCRIPT TO BUILD/UPDATE CKMRcpp PACKAGE
## 
## For more info see: 
## https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-package.pdf
## =============================================================================
library(Rcpp)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##                             FIRST TIME BUILD
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## 1. Create the skeleton in the correct location
Rcpp.package.skeleton(name = "CKMRcpp", path = "source")

## 2. Paste the cpp code (e.g., 'nllCKMRcppVanilla.cpp') in the folder 'scr' and 
## remove the file 'rcpp_hello_world.cpp'. Do not touch the 'RcppExports.cpp' 
## file. Instead, run the line below which updates this file.
compileAttributes("source/CKMRcpp")

## 3. Update the DESCRIPTION file if required

## 4. Remove the manual in folder 'man' for the template function

## 5. Build the package using the lines below
devtools::install("source/CKMRcpp")

## 6. BONUS: to create a tar.gz file which is easy to move around/share, run 
## the line below
devtools::build("source/CKMRcpp")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##                            UPDATE THE PACKAGE
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## 1. Check if the package exists and can be loaded
library(CKMRcpp)
## Unload the package
detach("package:Rcpp", unload=TRUE)

## 2. Make the sure the updated/new cpp files are added to "CKMRcpp/src" 

## 3. Re-compile 
compileAttributes("source/CKMRcpp")

## 4. Install the updated package
devtools::install("source/CKMRcpp")

## 5. BONUS: to create a tar.gz file which is easy to move around/share, run 
## the line below
devtools::build("source/CKMRcpp")