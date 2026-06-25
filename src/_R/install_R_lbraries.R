# =======================================================================
#  Authors: Aldo Hernandez-Corchado, Hamed S. Najafabadi
#
#  Copyright 2026 Aldo Hernandez-Corchado, Hamed S. Najafabadi
#
#  This file is part of MAGIX.
#
#  MAGIX is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MAGIX is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MAGIX.  If not, see <http://www.gnu.org/licenses/>.
# =======================================================================
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("GetoptLong")
install.packages("GlobalOptions")
install.packages("circlize")
BiocManager::install("ComplexHeatmap")

                 
install.packages("ggplot2")
install.packages("tictoc")
install.packages("Matrix")
install.packages("uwot")
install.packages("stringr")
install.packages("data.table")
install.packages("PRROC")
install.packages("sads")
install.packages("Rcpp")
install.packages("optparse")
