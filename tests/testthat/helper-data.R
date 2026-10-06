## Provenance: Test fixtures prepared with OpenAI Codex assistance.
make_se <- function(time, subject = NULL) {
    cd <- S4Vectors::DataFrame(time = time)
    if (!is.null(subject)) cd$subject <- subject
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(seq_len(3L * length(time)), nrow = 3,
            dimnames = list(paste0("g", 1:3), paste0("s", seq_along(time))))),
        rowData = S4Vectors::DataFrame(feature = letters[1:3]),
        colData = cd,
        metadata = list(study = "synthetic test")
    )
}
