library(pdftools)
library(here)
setwd(paste0(here(),"/Abstracts"))

pdf_combine(c("Abstracts1-Hans.pdf",
              "Abstracts2-Hans.pdf",
              "Abstracts3-Hans.pdf"),
            output = "abs-Hans.pdf")
