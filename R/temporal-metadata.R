# Provenance: Authored for OmicsChrono with OpenAI Codex assistance; maintained
# and validated by the package maintainer under Bioconductor's AI policy.
#' Validate temporal sample metadata
#'
#' Inspect a sample-time column for missing or malformed values, ambiguous civil
#' times, duplicate subject-time combinations, negative elapsed values, and
#' non-monotonic sample order. Civil timestamps must be ISO local times and need
#' an explicit IANA time zone. Numeric elapsed times need 
#' `time_unit`; tokens such as `"15 min"`, `"day3"`, `"T4"`, and `"baseline"`
#' are also supported.
#'
#' @param x A `SummarizedExperiment` or derived object.
#' @param time Name of a `colData(x)` column containing sample times.
#' @param subject Optional subject column name or one subject identifier per
#'   sample. If omitted, all samples are treated as one subject.
#' @param time_unit Unit for bare numeric values or `T` tokens: `"minutes"`,
#'   `"hours"`, or `"days"`.
#' @param timezone One IANA time zone or one zone per sample for civil times.
#'   Supply `NULL` for elapsed-time data or absolute `POSIXct` instants.
#' @param ambiguous How an ambiguous local civil timestamp is treated for
#'   diagnostics. The default `"error"` requires explicit disambiguation;
#'   `"earliest"` and `"latest"` select one occurrence and report a warning.
#'
#' @return An `S4Vectors::DataFrame` with one row per diagnostic and columns
#'   `sample`, `subject`, `code`, `severity`, `message`, and `value`. A report
#'   with zero rows means no issues were detected.
#'
#' @details This function reports metadata problems without changing `x`.
#'   Civil times use `%Y-%m-%d %H:%M:%S`; a `T` date/time separator is accepted.
#'   Numeric and unit-bearing elapsed-time tokens are normalized to hours only
#'   for validation. No assay data are read or modified.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = matrix(1:6, nrow = 2)),
#'     colData = S4Vectors::DataFrame(subject = c("a", "a", "b"),
#'                                    sample_time = c("0 min", "15 min", "baseline"))
#' )
#' validateSampleTime(se, "sample_time", subject = "subject")
#'
#' @export
validateSampleTime <- function(x, time, subject = NULL, time_unit = NULL,
                               timezone = NULL, ambiguous = "error") {
    cd <- .as_coldata(x)
    raw <- .column(cd, time, "time")
    n <- ncol(x)
    sample <- colnames(x)
    if (is.null(sample)) sample <- as.character(seq_len(n))
    subject_values <- .resolve_subject(cd, subject, n)
    subj <- if (is.null(subject_values)) rep(NA_character_, n) else subject_values
    grp <- .groups(subject_values, n)
    out <- .empty_issues()

    if (inherits(raw, "POSIXt")) {
        instant <- as.POSIXct(raw, tz = "UTC")
        codes <- rep(NA_character_, n)
        codes[is.na(instant)] <- "missing_time"
        out <- .append_issues(out, sample, subj, codes, value = as.character(raw))
        metric <- as.numeric(instant) / 3600
    } else if (is.numeric(raw) || is.factor(raw)) {
        parsed <- .parse_elapsed(if (is.factor(raw)) as.character(raw) else raw, time_unit)
        out <- .append_issues(out, sample, subj, parsed$issue, value = as.character(raw))
        metric <- parsed$hours
    } else {
        txt <- as.character(raw)
        civil <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}[ T]", txt)
        if (any(civil, na.rm = TRUE)) {
            parsed <- .parse_civil(txt, timezone, ambiguous)
            codes <- parsed$issue
            codes[!civil & is.na(codes)] <- "malformed_time"
            sev <- ifelse(codes == "ambiguous_time_resolved", "warning", "error")
            out <- .append_issues(out, sample, subj, codes, severity = sev, value = txt)
            metric <- as.numeric(parsed$utc) / 3600
        } else {
            parsed <- .parse_elapsed(txt, time_unit)
            out <- .append_issues(out, sample, subj, parsed$issue, value = txt)
            metric <- parsed$hours
        }
    }

    valid <- is.finite(metric)
    for (g in unique(grp)) {
        idx <- which(grp == g & valid)
        if (!length(idx)) next
        if (any(metric[idx] < metric[idx[[1L]]])) {
            codes <- rep(NA_character_, n)
            codes[idx[metric[idx] < metric[idx[[1L]]]]] <- "negative_elapsed_time"
            out <- .append_issues(out, sample, subj, codes, value = as.character(raw))
        }
        if (length(idx) > 1L && any(diff(metric[idx]) < 0)) {
            bad <- idx[c(FALSE, diff(metric[idx]) < 0)]
            codes <- rep(NA_character_, n)
            codes[bad] <- "non_monotonic_time"
            out <- .append_issues(out, sample, subj, codes, severity = "warning",
                                  value = as.character(raw))
        }
        duplicated <- duplicated(metric[idx]) | duplicated(metric[idx], fromLast = TRUE)
        if (any(duplicated)) {
            codes <- rep(NA_character_, n)
            codes[idx[duplicated]] <- "duplicate_subject_time"
            out <- .append_issues(out, sample, subj, codes, severity = "warning",
                                  value = as.character(raw))
        }
    }
    out
}

#' Annotate a SummarizedExperiment with normalized sample time
#'
#' Add normalized time fields to `colData(x)`, leaving assay values, feature
#' metadata, sample order, and the input class intact. Character civil
#' timestamps require `timezone`; DST-ambiguous times error by default and may
#' be explicitly resolved with `ambiguous = "earliest"` or `"latest"`.
#'
#' @inheritParams validateSampleTime
#' @param baseline Baseline rule. The default `"first"` sets elapsed time zero
#'   at the first valid sample for each subject, in the supplied sample order.
#' @param overwrite Replace derived columns with the same names. Defaults to
#'   `FALSE` so existing metadata is never silently overwritten.
#'
#' @return An object of the same class as `x`, with temporal fields appended to
#'   `colData(x)`.
#'
#' @details Elapsed-time inputs are converted to hours. Civil timestamps are
#'   parsed as local wall times in the supplied time zone(s), then differenced
#'   as absolute instants in UTC. Derived columns are
#'   `chrono_timestamp_utc`, `chrono_timestamp_local`, `chrono_elapsed_hours`,
#'   `chrono_study_day`, `chrono_day_index`, `chrono_clock_hour`,
#'   `chrono_day_of_week`, `chrono_is_dst`, and
#'   `chrono_sampling_interval_hours`. Non-monotonic input order is preserved
#'   and can be inspected with [validateSampleTime()].
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = matrix(1:6, nrow = 2)),
#'     colData = S4Vectors::DataFrame(subject = c("a", "a", "b"),
#'                                    minute = c(0, 15, 0))
#' )
#' out <- annotateSampleTime(se, "minute", subject = "subject", time_unit = "minutes")
#' SummarizedExperiment::colData(out)$chrono_elapsed_hours
#'
#' @export
annotateSampleTime <- function(x, time, subject = NULL, time_unit = NULL,
                               timezone = NULL, ambiguous = "error",
                               baseline = "first", overwrite = FALSE) {
    cd <- .as_coldata(x)
    raw <- .column(cd, time, "time")
    if (!identical(baseline, "first")) {
        stop("Only baseline = 'first' is currently supported.", call. = FALSE)
    }
    n <- ncol(x)
    subject_values <- .resolve_subject(cd, subject, n)
    subj <- if (is.null(subject_values)) rep(".all", n) else subject_values
    grp <- .groups(subject_values, n)
    is_civil <- !inherits(raw, "POSIXt") && !is.numeric(raw) &&
        any(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}[ T]", as.character(raw)), na.rm = TRUE)
    cols <- c("chrono_timestamp_utc", "chrono_timestamp_local",
              "chrono_elapsed_hours", "chrono_study_day", "chrono_day_index",
              "chrono_clock_hour", "chrono_day_of_week", "chrono_is_dst",
              "chrono_sampling_interval_hours")
    .check_new_columns(cd, cols, overwrite)

    if (inherits(raw, "POSIXt")) {
        utc <- as.POSIXct(raw, tz = "UTC")
        input_zone <- attr(raw, "tzone", exact = TRUE)
        if (is.null(input_zone) || !length(input_zone) || !nzchar(input_zone[[1L]])) input_zone <- "UTC"
        zones <- if (is.null(timezone)) rep(input_zone[[1L]], n) else
            if (length(timezone) == 1L) rep(timezone, n) else timezone
        issue <- ifelse(is.na(utc), "missing_time", NA_character_)
    } else if (is_civil) {
        parsed <- .parse_civil(raw, timezone, ambiguous)
        utc <- parsed$utc
        zones <- parsed$timezone
        issue <- parsed$issue
    } else {
        parsed <- .parse_elapsed(if (is.factor(raw)) as.character(raw) else raw, time_unit)
        elapsed <- parsed$hours
        issue <- parsed$issue
        utc <- as.POSIXct(rep(NA_real_, n), origin = "1970-01-01", tz = "UTC")
        zones <- rep(NA_character_, n)
    }
    errors <- unique(issue[!is.na(issue) & issue != "ambiguous_time_resolved"])
    if (length(errors)) {
        stop("Cannot annotate sample time; resolve: ", paste(errors, collapse = ", "),
             ". See validateSampleTime() for row-level diagnostics.", call. = FALSE)
    }
    if (is_civil || inherits(raw, "POSIXt")) {
        elapsed <- rep(NA_real_, n)
        for (g in unique(grp)) {
            idx <- which(grp == g & is.finite(as.numeric(utc)))
            if (length(idx)) elapsed[idx] <- as.numeric(difftime(utc[idx], utc[idx[[1L]]], units = "hours"))
        }
    }
    interval <- rep(NA_real_, n)
    for (g in unique(grp)) {
        idx <- which(grp == g & is.finite(elapsed))
        if (length(idx) > 1L) interval[idx[-1L]] <- diff(elapsed[idx])
    }
    local <- rep(NA_character_, n)
    clock_hour <- rep(NA_real_, n)
    weekday <- rep(NA_integer_, n)
    is_dst <- rep(NA, n)
    if (is_civil || inherits(raw, "POSIXt")) {
        for (zone in unique(zones[!is.na(zones)])) {
            idx <- which(zones == zone & is.finite(as.numeric(utc)))
            if (!length(idx)) next
            local[idx] <- format(utc[idx], "%Y-%m-%d %H:%M:%S", tz = zone)
            parts <- as.POSIXlt(utc[idx], tz = zone)
            clock_hour[idx] <- parts$hour + parts$min / 60 + parts$sec / 3600
            weekday[idx] <- as.integer(format(utc[idx], "%u", tz = zone))
            is_dst[idx] <- parts$isdst > 0
        }
    }
    if (is_civil || inherits(raw, "POSIXt")) {
        elapsed_out <- elapsed
        day_index <- floor(elapsed / 24)
    } else {
        elapsed_out <- elapsed
        day_index <- floor(elapsed / 24)
    }
    cd$chrono_timestamp_utc <- utc
    attr(cd$chrono_timestamp_utc, "tzone") <- "UTC"
    cd$chrono_timestamp_local <- local
    cd$chrono_elapsed_hours <- elapsed_out
    cd$chrono_study_day <- ifelse(is.na(day_index), NA_integer_, as.integer(day_index + 1))
    cd$chrono_day_index <- ifelse(is.na(day_index), NA_integer_, as.integer(day_index))
    cd$chrono_clock_hour <- clock_hour
    cd$chrono_day_of_week <- weekday
    cd$chrono_is_dst <- is_dst
    cd$chrono_sampling_interval_hours <- interval
    SummarizedExperiment::colData(x) <- cd
    x
}

#' Add Zeitgeber and photoperiod annotations
#'
#' Compute time since an explicitly specified local lights-on time and classify
#' each sample as light or dark under the supplied controlled photoperiod.
#'
#' @param x A `SummarizedExperiment` or derived object.
#' @param timestamp Name of a civil-time column in `colData(x)`.
#' @param lights_on Local clock time in `HH:MM` or `HH:MM:SS` format.
#' @param timezone One IANA time zone or one zone per sample.
#' @param photoperiod Number of hours per cycle in the light phase; must be
#'   between zero and 24.
#' @param ambiguous DST fall-back policy as in [annotateSampleTime()].
#' @param overwrite Replace existing `chrono_zeitgeber_hours`,
#'   `chrono_light_phase`, and `chrono_hours_until_phase_change` columns.
#'
#' @return An object of the same class with three Zeitgeber fields in `colData`.
#'
#' @details ZT0 is defined by the provided local lights-on clock time. The
#'   calculation uses local clock coordinates modulo 24 hours and the supplied
#'   light duration; it does not infer biological phase or query external data.
#'   Studies with non-24-hour cycles should use a separately specified
#'   experimental schedule.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = matrix(1:2, nrow = 1)),
#'     colData = S4Vectors::DataFrame(collection =
#'         c("2024-06-01 07:00:00", "2024-06-01 19:00:00"))
#' )
#' out <- addZeitgeberTime(se, "collection", lights_on = "07:00",
#'                         timezone = "Europe/Oslo", photoperiod = 12)
#' SummarizedExperiment::colData(out)$chrono_zeitgeber_hours
#'
#' @export
addZeitgeberTime <- function(x, timestamp, lights_on, timezone,
                             photoperiod = 12, ambiguous = "error",
                             overwrite = FALSE) {
    cd <- .as_coldata(x)
    raw <- .column(cd, timestamp, "timestamp")
    if (length(photoperiod) != 1L || !is.numeric(photoperiod) ||
        !is.finite(photoperiod) || photoperiod <= 0 || photoperiod >= 24) {
        stop("'photoperiod' must be a number strictly between 0 and 24 hours.", call. = FALSE)
    }
    if (!grepl("^([01][0-9]|2[0-3]):[0-5][0-9](:[0-5][0-9])?$", lights_on)) {
        stop("'lights_on' must use HH:MM or HH:MM:SS.", call. = FALSE)
    }
    cols <- c("chrono_zeitgeber_hours", "chrono_light_phase",
              "chrono_hours_until_phase_change")
    .check_new_columns(cd, cols, overwrite)
    parsed <- .parse_civil(raw, timezone, ambiguous)
    bad <- !is.na(parsed$issue) & parsed$issue != "ambiguous_time_resolved"
    if (any(bad)) {
        stop("Cannot compute Zeitgeber time; resolve: ",
             paste(unique(parsed$issue[bad]), collapse = ", "), ".", call. = FALSE)
    }
    zones <- parsed$timezone
    zt <- rep(NA_real_, length(raw))
    for (zone in unique(zones)) {
        idx <- which(zones == zone & is.finite(as.numeric(parsed$utc)))
        if (!length(idx)) next
        local <- as.POSIXlt(parsed$utc[idx], tz = zone)
        clock_hour <- local$hour + local$min / 60 + local$sec / 3600
        onset <- as.numeric(substr(lights_on, 1L, 2L)) +
            as.numeric(substr(lights_on, 4L, 5L)) / 60 +
            if (nchar(lights_on) == 8L) as.numeric(substr(lights_on, 7L, 8L)) / 3600 else 0
        zt[idx] <- (clock_hour - onset) %% 24
    }
    phase <- ifelse(is.na(zt), NA_character_, ifelse(zt < photoperiod, "light", "dark"))
    until <- ifelse(is.na(zt), NA_real_, ifelse(zt < photoperiod, photoperiod - zt, 24 - zt))
    cd$chrono_zeitgeber_hours <- zt
    cd$chrono_light_phase <- phase
    cd$chrono_hours_until_phase_change <- until
    SummarizedExperiment::colData(x) <- cd
    x
}

#' Align observed sample times to expected sampling windows
#'
#' Match samples to expected time points using a deterministic one-to-one
#' nearest-offset assignment. Samples and expected times use the same numeric
#' unit. Ties are resolved by input sample order and then expected-time order.
#'
#' @param observed Numeric vector of observed times.
#' @param expected Numeric vector of expected times.
#' @param tolerance Maximum absolute difference for a match, in the same unit
#'   as `observed` and `expected`.
#' @param subject Optional subject identifier per observation.
#' @param sample Optional sample identifier per observation.
#'
#' @return An `S4Vectors::DataFrame` with one row per expected slot and one row
#'   per unmatched observed sample. `status` is `matched`, `missing`, or
#'   `unexpected`.
#'
#' @examples
#' alignSampleTimes(c(0.1, 5.8, 12.3), c(0, 6, 12), tolerance = 0.5,
#'                  subject = c("a", "a", "a"))
#'
#' @export
alignSampleTimes <- function(observed, expected, tolerance, subject = NULL,
                             sample = NULL) {
    if (!is.numeric(observed) || !is.numeric(expected) || any(!is.finite(expected))) {
        stop("'observed' and 'expected' must be numeric; expected values must be finite.", call. = FALSE)
    }
    if (length(tolerance) != 1L || !is.numeric(tolerance) || !is.finite(tolerance) || tolerance < 0) {
        stop("'tolerance' must be one finite non-negative number.", call. = FALSE)
    }
    n <- length(observed)
    if (is.null(sample)) sample <- as.character(seq_len(n))
    if (length(sample) != n) stop("'sample' must have one value per observation.", call. = FALSE)
    subj <- if (is.null(subject)) rep(".all", n) else as.character(subject)
    grp <- if (n == 0L && is.null(subject)) ".all" else .groups(subject, n)
    if (length(subj) != n) stop("'subject' must have one value per observation.", call. = FALSE)
    rows <- list()
    for (g in unique(grp)) {
        all_o <- if (n) which(grp == g) else integer()
        oi <- all_o[is.finite(observed[all_o])]
        ei <- seq_along(expected)
        pairs <- expand.grid(observed = oi, expected = ei)
        pairs$offset <- abs(observed[pairs$observed] - expected[pairs$expected])
        pairs <- pairs[pairs$offset <= tolerance, , drop = FALSE]
        if (nrow(pairs)) {
            pairs <- pairs[order(pairs$offset, pairs$observed, pairs$expected), , drop = FALSE]
        }
        used_o <- integer()
        used_e <- integer()
        assignment <- stats::setNames(integer(length(expected)), as.character(seq_along(expected)))
        for (k in seq_len(nrow(pairs))) {
            o <- pairs$observed[[k]]
            e <- pairs$expected[[k]]
            if (o %in% used_o || e %in% used_e) next
            used_o <- c(used_o, o)
            used_e <- c(used_e, e)
            assignment[[as.character(e)]] <- o
        }
        for (e in ei) {
            if (e %in% used_e) {
                o <- assignment[[as.character(e)]]
                rows[[length(rows) + 1L]] <- data.frame(subject = subj[if (length(oi)) oi[[1L]] else 1L],
                    sample = as.character(sample[[o]]), expected_time = expected[[e]],
                    observed_time = observed[[o]], offset = observed[[o]] - expected[[e]],
                    within_tolerance = TRUE, status = "matched", stringsAsFactors = FALSE)
            } else {
                rows[[length(rows) + 1L]] <- data.frame(subject = if (length(oi)) subj[[oi[[1L]]]] else g,
                    sample = NA_character_, expected_time = expected[[e]], observed_time = NA_real_,
                    offset = NA_real_, within_tolerance = FALSE, status = "missing", stringsAsFactors = FALSE)
            }
        }
        un <- setdiff(all_o, used_o)
        for (o in un) {
            rows[[length(rows) + 1L]] <- data.frame(subject = subj[[o]], sample = as.character(sample[[o]]),
                expected_time = NA_real_, observed_time = observed[[o]], offset = NA_real_,
                within_tolerance = FALSE, status = "unexpected", stringsAsFactors = FALSE)
        }
    }
    if (!length(rows)) {
        return(S4Vectors::DataFrame(subject = character(), sample = character(),
            expected_time = numeric(), observed_time = numeric(), offset = numeric(),
            within_tolerance = logical(), status = character()))
    }
    S4Vectors::DataFrame(do.call(rbind, rows), check.names = FALSE)
}

#' Summarize temporal coverage
#'
#' Report sample counts, temporal ranges, intervals, and interval variability
#' by subject and optional condition from the sample metadata in a
#' `SummarizedExperiment`.
#'
#' @param x A `SummarizedExperiment` or derived object.
#' @param time Name of a numeric `colData(x)` time column.
#' @param subject Optional subject column name or one subject identifier per
#'   sample.
#' @param condition Optional condition column.
#'
#' @return An `S4Vectors::DataFrame` with one row per subject/condition group.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = matrix(1:8, nrow = 2)),
#'     colData = S4Vectors::DataFrame(subject = rep(c("a", "b"), each = 2),
#'                                    hour = c(0, 6, 0, 12))
#' )
#' temporalCoverage(se, "hour", subject = "subject")
#'
#' @export
temporalCoverage <- function(x, time, subject = NULL, condition = NULL) {
    cd <- .as_coldata(x)
    value <- .column(cd, time, "time")
    if (!is.numeric(value)) stop("'time' must refer to a numeric colData column.", call. = FALSE)
    subject_values <- .resolve_subject(cd, subject, ncol(x))
    subj <- if (is.null(subject_values)) rep(".all", ncol(x)) else subject_values
    cond <- if (is.null(condition)) rep(".all", ncol(x)) else as.character(.column(cd, condition, "condition"))
    group <- interaction(subj, cond, drop = TRUE, lex.order = TRUE, sep = "\r")
    groups <- levels(group)
    result <- lapply(groups, function(key) {
        idx <- which(group == key & is.finite(value))
        vals <- value[idx]
        intervals <- if (length(vals) > 1L) diff(sort(vals)) else numeric()
        cv <- if (length(intervals) > 1L && mean(intervals) != 0) stats::sd(intervals) / mean(intervals) else NA_real_
        data.frame(subject = if (length(idx)) subj[idx[[1L]]] else NA_character_,
                   condition = if (length(idx)) cond[idx[[1L]]] else NA_character_,
                   n_samples = length(idx), n_unique_times = length(unique(vals)),
                   time_start = if (length(vals)) min(vals) else NA_real_,
                   time_end = if (length(vals)) max(vals) else NA_real_,
                   median_interval = if (length(intervals)) stats::median(intervals) else NA_real_,
                   min_interval = if (length(intervals)) min(intervals) else NA_real_,
                   max_interval = if (length(intervals)) max(intervals) else NA_real_,
                   interval_cv = cv, stringsAsFactors = FALSE)
    })
    if (!length(result)) return(S4Vectors::DataFrame())
    S4Vectors::DataFrame(do.call(rbind, result), check.names = FALSE)
}

#' Construct explicit temporal design covariates
#'
#' Copy selected sample-level time, subject, and condition fields into a
#' `DataFrame`, optionally adding sine and cosine terms for a user-specified
#' cycle length. No model is fitted and phase is not inferred.
#'
#' @param x A `SummarizedExperiment` or derived object.
#' @param time Name of a numeric time column, normally
#'   `"chrono_elapsed_hours"`.
#' @param subject Optional subject column name or one subject identifier per
#'   sample.
#' @param condition Optional condition column name.
#' @param cycle_length Optional positive numeric cycle length in the same unit
#'   as `time`. When supplied, include cyclical sine/cosine covariates.
#'
#' @return An `S4Vectors::DataFrame` in original sample order.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = matrix(1:4, nrow = 1)),
#'     colData = S4Vectors::DataFrame(hour = c(0, 6, 12, 18),
#'                                    subject = rep("mouse1", 4))
#' )
#' temporalDesignData(se, "hour", subject = "subject", cycle_length = 24)
#'
#' @export
temporalDesignData <- function(x, time, subject = NULL, condition = NULL,
                               cycle_length = NULL) {
    cd <- .as_coldata(x)
    values <- .column(cd, time, "time")
    if (!is.numeric(values)) stop("'time' must refer to a numeric colData column.", call. = FALSE)
    out <- S4Vectors::DataFrame(time = values)
    if (!is.null(subject)) out$subject <- .column(cd, subject, "subject")
    if (!is.null(condition)) out$condition <- .column(cd, condition, "condition")
    if (!is.null(cycle_length)) {
        if (length(cycle_length) != 1L || !is.numeric(cycle_length) ||
            !is.finite(cycle_length) || cycle_length <= 0) {
            stop("'cycle_length' must be one positive finite number.", call. = FALSE)
        }
        angle <- 2 * pi * values / cycle_length
        out$cycle_sin <- sin(angle)
        out$cycle_cos <- cos(angle)
    }
    out
}
