#-----------------------------------------------------------------------
# This file was retrieved from
# <https://github.com/n-kall/bayesian-workflow-wiki/blob/67fd1569f08d6d516cb18de99d5175f6950b0559/R/tikzpicture-template.R>
# on November 27, 2022. The corresponding license was:

# MIT License
# 
# Copyright (c) 2022 n-kall
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#-----------------------------------------------------------------------

require(tikzDevice)

# define plot saving method
save_tikz_plot <- function(
    plot, filename, width = NA, height = NA, asp = NA
) {
  # automatic scaling
  if (is.na(asp)) asp <- 1.618
  if (is.na(width) && is.na(height)) {
    height <- 3.71
    width <- height * asp
  }
  else if (is.na(width)) {
    width <- height * asp
  }
  else if (is.na(height)) {
    height <- width / asp
  }
  
  # make tex
  tikz(file = filename, width = width, height = height)
  print(plot)
  dev.off()
  
  # patch cropping issues
  lines <- readLines(con = filename)
  lines <- lines[-which(grepl("\\path\\[clip\\]*", lines))]
  lines <- lines[-which(grepl("\\path\\[use as bounding box*", lines))]
  writeLines(lines, con = filename)
}
