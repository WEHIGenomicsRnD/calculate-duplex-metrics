# source("renv/activate.R")

# Set options for better CLI experience
options(
  repos = c(CRAN = "https://cran.rstudio.com/"),
  warn = 1,
  error = quote({
    if (!interactive()) {
      cat("Error occurred, exiting...\n")
      quit(status = 1)
    }
  })
)
