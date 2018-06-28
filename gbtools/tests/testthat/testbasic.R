context("Basic tests")
library(gbtools)

test_that("Input format validator catches errors", {
  expect_true(file.exists(file.path("../testdata","SampleG1_error.covstats")))
  expect_message(gbt_checkinput(covstats=file.path("../testdata","SampleG1_error.covstats"),
                                ssu=file.path("../testdata","olavius_metagenome.ssu.tab")),
                 "Errors found in input files. Check validator log for details")
  expect_message(gbt_checkinput(covstats=c(file.path("../testdata","SampleG1.covstats"),
                                           file.path("../testdata","SampleA2.covstats")),
                                ssu=file.path("../testdata","olavius_metagenome.ssu.tab"),
                                mark=c(file.path("../testdata","amphora2_results.tab"),
                                       file.path("../testdata","blobology_results.tab"))
                                ),
               "No errors found in input files")
  expect_true(file.exists(file.path(".","input_validator.log")))
})


  expect_warning(e <- gbt(covstats=c(file.path("../testdata","SampleA2.covstats"),
                                     file.path("../testdata","SampleG1.covstats")),
                          ssu=file.path("../testdata","olavius_metagenome.ssu.tab"),
                          mark=c(file.path("../testdata","amphora2_results.tab"),
                                 file.path("../testdata","blobology_results.tab"))
                          )
                 )

test_that("Data can be imported", {

  d <- gbt(covstats=c(file.path("../testdata","SampleG1.covstats"),
                    file.path("../testdata","SampleA2.covstats")),
         ssu=file.path("../testdata","olavius_metagenome.ssu.tab"),
         mark=c(file.path("../testdata","amphora2_results.tab"),
                file.path("../testdata","blobology_results.tab")),
         marksource=c("amphora2","blob")
         )
  expect_is(d,"gbt")
})

# Cleanup
file.remove(file.path(".","input_validator.log"))