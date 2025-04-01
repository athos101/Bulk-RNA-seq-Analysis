## ABOUT THE FILE
#//=======================================================================//

# Author: Athos B. Schuck

## STEP 1: LOAD LIBRARIES
#//=======================================================================//
library(Biobase)
library(SeuratObject)
suppressWarnings(library(Seurat))
suppressWarnings(library(BayesPrism))
library(immunedeconv)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(readr)
library(Matrix)


## STEP 2: LOAD DATA
#//=======================================================================//


## STEP 3: FILTERING MARKERS WITH MEANRATIO
#//=======================================================================//

