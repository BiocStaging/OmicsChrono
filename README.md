# OmicsChrono

OmicsChrono validates and standardizes temporal sample metadata for longitudinal and circadian high-throughput experiments. It works in the `colData` of `SummarizedExperiment` objects and derived classes, including `SingleCellExperiment`, and leaves assays, feature metadata, and sample order intact.

## Why temporal metadata needs explicit handling

Time-course studies may combine numeric elapsed times, civil timestamps, multiple subjects, and controlled light schedules. Treating these as interchangeable can misalign samples, obscure irregular sampling, or silently change elapsed durations around daylight-saving transitions. OmicsChrono makes the representation and assumptions explicit before downstream modeling.

## Installation

During review, install from source with `BiocManager::install("OmicsChrono")` after the package enters Bioconductor devel. The package is not yet available from Bioconductor.

## Quick example

```r
library(SummarizedExperiment)
library(OmicsChrono)

data <- SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 3)),
    colData = DataFrame(
        subject = c("A", "A", "A", "B"),
        minute = c(0, 15, 30, 0)
    )
)
validateSampleTime(data, time = "minute", subject = "subject", time_unit = "minutes")
data <- annotateSampleTime(
    data, time = "minute", subject = "subject", time_unit = "minutes"
)
temporalCoverage(data, time = "chrono_elapsed_hours", subject = "subject")
```

## Civil time and experimental time

Civil timestamps identify instants in a named time zone. Numeric experimental times describe durations from a study baseline. These representations are kept separate; converting a civil timestamp requires an explicit IANA time zone, and ambiguous or nonexistent daylight-saving local times are reported.

## Zeitgeber time

`addZeitgeberTime()` derives time since an explicitly supplied lights-on clock time. It does not infer circadian phase from expression data. Photoperiod is supplied by the user as a controlled experimental setting.

## Sampling windows and temporal coverage

`alignSampleTimes()` performs deterministic one-to-one matching to expected windows within a tolerance. `temporalCoverage()` reports observed ranges, intervals, interval variability, and sample counts by subject and condition.

## Interoperability

Temporal annotations are added to `colData()`. The package does not copy or transform assay values, `rowData`, `rowRanges`, `metadata`, or reduced dimensions.

## Scope and non-goals

OmicsChrono is an upstream metadata engineering layer. It does not test rhythmicity, infer molecular phase, or fit differential-expression models. Those analyses remain with packages such as `rain`, `DiscoRhythm`, and `timeOmics`.

## Contributing and license

Please report issues in the GitHub issue tracker. OmicsChrono is distributed
under the MIT license. The initial implementation was prepared with OpenAI
Codex assistance; the maintainer is responsible for review, reliability, and
ongoing maintenance.
