## Provenance: Tests drafted with OpenAI Codex assistance and reviewed against
## the package API for this contribution.
test_that("elapsed tokens normalize by subject baseline and preserve the container", {
    se <- make_se(c("0 min", "15 min", "30 min", "baseline"),
                  c("a", "a", "a", "b"))
    assay_before <- SummarizedExperiment::assay(se)
    rowdata_before <- SummarizedExperiment::rowData(se)
    sample_names <- colnames(se)
    out <- annotateSampleTime(se, "time", subject = "subject")
    expect_s4_class(out, "SummarizedExperiment")
    expect_identical(SummarizedExperiment::assay(out), assay_before)
    expect_identical(SummarizedExperiment::rowData(out), rowdata_before)
    expect_identical(colnames(out), sample_names)
    expect_equal(SummarizedExperiment::colData(out)$chrono_elapsed_hours,
                 c(0, 0.25, 0.5, 0))
    expect_s4_class(validateSampleTime(se, "time", subject = "subject"), "DataFrame")
})

test_that("numeric elapsed times require units and report missing and duplicates", {
    se <- make_se(c(0, NA, 0, 12), c("a", "a", "a", "a"))
    report <- validateSampleTime(se, "time", subject = "subject", time_unit = "hours")
    expect_true(all(c("missing_time", "duplicate_subject_time") %in% report$code))
    expect_error(annotateSampleTime(se, "time", subject = "subject"), "unit_required")
})

test_that("common elapsed-time spellings are parsed only with justified units", {
    se <- make_se(c("baseline", "15 min", "day3", "T4"))
    report <- validateSampleTime(se, "time", time_unit = "hours")
    expect_false(any(report$severity == "error"))
    expect_true("unit_required" %in% validateSampleTime(se, "time")$code)
})

test_that("civil timestamps handle UTC, leap day, and year boundaries", {
    se <- make_se(c("2024-12-31 23:30:00", "2025-01-01T00:30:00",
                    "2024-02-29 12:00:00"))
    out <- annotateSampleTime(se, "time", timezone = "UTC")
    expect_equal(SummarizedExperiment::colData(out)$chrono_elapsed_hours[1:2], c(0, 1))
    expect_lt(SummarizedExperiment::colData(out)$chrono_elapsed_hours[[3]], 0)
    expect_equal(SummarizedExperiment::colData(out)$chrono_timestamp_local[[2]],
                 "2025-01-01 00:30:00")
})

test_that("DST gaps are rejected and repeated clock times require disambiguation", {
    spring <- make_se("2024-03-31 02:30:00")
    spring_report <- validateSampleTime(spring, "time", timezone = "Europe/Oslo")
    expect_true("nonexistent_local_time" %in% spring_report$code)
    expect_error(annotateSampleTime(spring, "time", timezone = "Europe/Oslo"),
                 "nonexistent_local_time")

    fall <- make_se(c("2024-10-27 02:30:00", "2024-10-27 02:30:00"))
    report <- validateSampleTime(fall, "time", timezone = "Europe/Oslo")
    expect_true("ambiguous_local_time" %in% report$code)
    expect_error(annotateSampleTime(fall, "time", timezone = "Europe/Oslo"),
                 "ambiguous_local_time")
    early <- annotateSampleTime(fall, "time", timezone = "Europe/Oslo",
                                ambiguous = "earliest")
    late <- annotateSampleTime(fall, "time", timezone = "Europe/Oslo",
                               ambiguous = "latest")
    early_utc <- SummarizedExperiment::colData(early)$chrono_timestamp_utc
    late_utc <- SummarizedExperiment::colData(late)$chrono_timestamp_utc
    expect_equal(as.numeric(late_utc[[1]] - early_utc[[1]]), 1)
})

test_that("New York daylight-saving rules are treated as civil-time diagnostics", {
    spring <- make_se("2024-03-10 02:30:00")
    expect_true("nonexistent_local_time" %in%
        validateSampleTime(spring, "time", timezone = "America/New_York")$code)
    fall <- make_se("2024-11-03 01:30:00")
    expect_true("ambiguous_local_time" %in%
        validateSampleTime(fall, "time", timezone = "America/New_York")$code)
})

test_that("Zeitgeber time uses an explicit light schedule", {
    se <- make_se(c("2024-06-01 07:00:00", "2024-06-01 18:00:00",
                    "2024-06-01 20:00:00"))
    out <- addZeitgeberTime(se, "time", lights_on = "07:00",
                            timezone = "Europe/Oslo", photoperiod = 12)
    cd <- SummarizedExperiment::colData(out)
    expect_equal(cd$chrono_zeitgeber_hours, c(0, 11, 13))
    expect_equal(cd$chrono_light_phase, c("light", "light", "dark"))
    expect_equal(cd$chrono_hours_until_phase_change, c(12, 1, 11))
})

test_that("window matching is deterministic and one-to-one", {
    aligned <- alignSampleTimes(c(0.1, 0.2, 6.1, 30), c(0, 6, 12),
                               tolerance = 0.5, subject = rep("mouse", 4),
                               sample = paste0("s", 1:4))
    expect_equal(sum(aligned$status == "matched"), 2)
    expect_equal(sum(aligned$status == "missing"), 1)
    expect_equal(sum(aligned$status == "unexpected"), 2)
    expect_equal(sum(!is.na(aligned$sample[aligned$status == "matched"])), 2)
})

test_that("coverage and design covariates summarize irregular sampling", {
    se <- make_se(c(0, 6, 18, 0, 12), c("a", "a", "a", "b", "b"))
    out <- annotateSampleTime(se, "time", subject = "subject", time_unit = "hours")
    cov <- temporalCoverage(out, "chrono_elapsed_hours", subject = "subject")
    expect_equal(cov$n_samples, c(3, 2))
    expect_equal(cov$median_interval, c(9, 12))
    design <- temporalDesignData(out, "chrono_elapsed_hours", "subject", cycle_length = 24)
    expect_true(all(c("cycle_sin", "cycle_cos") %in% names(design)))
    expect_equal(nrow(design), ncol(se))
})

test_that("annotations refuse to overwrite existing temporal fields", {
    se <- make_se(0)
    out <- annotateSampleTime(se, "time", time_unit = "hours")
    expect_error(annotateSampleTime(out, "time", time_unit = "hours"), "Refusing to overwrite")
    expect_s4_class(annotateSampleTime(out, "time", time_unit = "hours", overwrite = TRUE),
                    "SummarizedExperiment")
})

test_that("SingleCellExperiment subclasses retain cell-specific structures", {
    skip_if_not_installed("SingleCellExperiment")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(seq_len(8), nrow = 2)),
        colData = S4Vectors::DataFrame(time = 0:3),
        reducedDims = S4Vectors::SimpleList(PCA = matrix(seq_len(8), ncol = 2))
    )
    out <- annotateSampleTime(sce, "time", time_unit = "hours")
    expect_s4_class(out, "SingleCellExperiment")
    expect_identical(SingleCellExperiment::reducedDims(out), SingleCellExperiment::reducedDims(sce))
    expect_identical(SummarizedExperiment::assay(out), SummarizedExperiment::assay(sce))
})

test_that("the fission experiment data load and retain assay contents", {
    skip_if_not_installed("fission")
    data("fission", package = "fission", envir = environment())
    expect_s4_class(fission, "RangedSummarizedExperiment")
    expect_true(all(c("minute", "strain") %in% names(SummarizedExperiment::colData(fission))))
    fission_subject <- paste(fission$strain, fission$replicate, sep = "_")
    ## Reconstruct the packaged snapshot in the currently installed container
    ## class; the experiment package serializes an older assay wrapper.
    current_fission <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(fission)),
        rowRanges = SummarizedExperiment::rowRanges(fission),
        colData = S4Vectors::DataFrame(as.data.frame(SummarizedExperiment::colData(fission)))
    )
    before <- SummarizedExperiment::assay(current_fission)
    out <- annotateSampleTime(current_fission, "minute", subject = fission_subject,
                              time_unit = "minutes")
    expect_identical(SummarizedExperiment::assay(out), before)
    expect_equal(sort(unique(SummarizedExperiment::colData(out)$chrono_elapsed_hours)),
                 c(0, 0.25, 0.5, 1, 2, 3))
})

test_that("zero-column experiments and empty observations have defined results", {
    empty <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(numeric(), nrow = 2, ncol = 0)),
        colData = S4Vectors::DataFrame(time = numeric(), row.names = character())
    )
    expect_s4_class(annotateSampleTime(empty, "time", time_unit = "hours"),
                    "SummarizedExperiment")
    aligned <- alignSampleTimes(numeric(), c(0, 6), tolerance = 0)
    expect_equal(sum(aligned$status == "missing"), 2)
})
