# Provenance: Authored for OmicsChrono with OpenAI Codex assistance; maintained
# and validated by the package maintainer under Bioconductor's AI policy.
.as_coldata <- function(x) {
    if (!methods::is(x, "SummarizedExperiment")) {
        stop("'x' must be a SummarizedExperiment or a class derived from it.", call. = FALSE)
    }
    SummarizedExperiment::colData(x)
}

.column <- function(cd, name, what) {
    if (length(name) != 1L || !is.character(name) || !name %in% names(cd)) {
        stop("'", what, "' must name a column in colData(x).", call. = FALSE)
    }
    cd[[name]]
}

.resolve_subject <- function(cd, subject, n) {
    if (is.null(subject)) return(NULL)
    if (is.character(subject) && length(subject) == 1L && subject %in% names(cd)) {
        subject <- cd[[subject]]
    }
    if (length(subject) != n) {
        stop("'subject' must name a colData column or have one value per sample.", call. = FALSE)
    }
    as.character(subject)
}

.groups <- function(subject, n) {
    if (is.null(subject)) {
        return(rep(".all", n))
    }
    subject <- as.character(subject)
    if (length(subject) != n) {
        stop("'subject' must have one value per sample.", call. = FALSE)
    }
    missing <- is.na(subject) | !nzchar(subject)
    subject[missing] <- paste0(".missing_subject_", which(missing))
    subject
}

.to_hours <- function(x, unit) {
    factor <- switch(unit, minutes = 1 / 60, hours = 1, days = 24,
                     stop("'time_unit' must be 'minutes', 'hours', or 'days'.", call. = FALSE))
    as.numeric(x) * factor
}

.parse_elapsed <- function(x, time_unit = NULL) {
    n <- length(x)
    hours <- rep(NA_real_, n)
    issue <- rep(NA_character_, n)
    if (is.numeric(x)) {
        if (is.null(time_unit)) {
            issue[] <- "unit_required"
            return(list(hours = hours, issue = issue))
        }
        hours <- tryCatch(.to_hours(x, time_unit), error = function(e) e)
        if (inherits(hours, "error")) {
            stop(conditionMessage(hours), call. = FALSE)
        }
        issue[is.na(x)] <- "missing_time"
        issue[!is.na(x) & !is.finite(x)] <- "non_finite_time"
        hours[!is.finite(hours)] <- NA_real_
        return(list(hours = hours, issue = issue))
    }
    x <- trimws(as.character(x))
    issue[is.na(x) | !nzchar(x)] <- "missing_time"
    for (i in which(is.na(issue))) {
        token <- tolower(x[[i]])
        if (token %in% c("baseline", "base", "t0")) {
            hours[[i]] <- 0
            next
        }
        if (grepl("^[-+]?[0-9]*\\.?[0-9]+$", token)) {
            if (!is.null(time_unit)) {
                hours[[i]] <- .to_hours(as.numeric(token), time_unit)
            } else {
                issue[[i]] <- "unit_required"
            }
            next
        }
        match <- regmatches(token, regexec("^([-+]?[0-9]*\\.?[0-9]+)\\s*(min|mins|minute|minutes|m|h|hr|hrs|hour|hours|d|day|days)$", token))[[1L]]
        if (length(match)) {
            value <- as.numeric(match[[2L]])
            unit <- match[[3L]]
            multiplier <- if (unit %in% c("m", "min", "mins", "minute", "minutes")) 1 / 60 else if (unit %in% c("d", "day", "days")) 24 else 1
            hours[[i]] <- value * multiplier
            next
        }
        match <- regmatches(token, regexec("^day\\s*([-+]?[0-9]*\\.?[0-9]+)$", token))[[1L]]
        if (length(match)) {
            hours[[i]] <- as.numeric(match[[2L]]) * 24
            next
        }
        match <- regmatches(token, regexec("^t\\s*([-+]?[0-9]*\\.?[0-9]+)$", token))[[1L]]
        if (length(match) && !is.null(time_unit)) {
            hours[[i]] <- .to_hours(as.numeric(match[[2L]]), time_unit)
            next
        }
        issue[[i]] <- if (grepl("^t", token)) "unit_required" else "malformed_time"
    }
    list(hours = hours, issue = issue)
}

.parse_civil <- function(x, timezone, ambiguous = "error") {
    x <- as.character(x)
    n <- length(x)
    zones <- if (length(timezone) == 1L) rep(timezone, n) else timezone
    if (length(zones) != n) {
        stop("'timezone' must be a single IANA time zone or have one value per sample.", call. = FALSE)
    }
    utc <- as.POSIXct(rep(NA_real_, n), origin = "1970-01-01", tz = "UTC")
    issue <- rep(NA_character_, n)
    if (is.null(timezone)) {
        issue[] <- "timezone_required"
        return(list(utc = utc, issue = issue, timezone = zones))
    }
    if (!ambiguous %in% c("error", "earliest", "latest")) {
        stop("'ambiguous' must be 'error', 'earliest', or 'latest'.", call. = FALSE)
    }
    missing_zone <- is.na(zones) | !nzchar(zones)
    issue[missing_zone & is.na(issue)] <- "invalid_timezone"
    x <- gsub("T", " ", trimws(x), fixed = TRUE)
    issue[is.na(x) | !nzchar(x)] <- "missing_time"
    for (zone in unique(zones[!missing_zone])) {
        idx <- which(zones == zone & is.na(issue))
        if (!length(idx)) next
        if (is.na(zone) || !zone %in% OlsonNames()) {
            issue[idx] <- "invalid_timezone"
            next
        }
        raw <- x[idx]
        parsed_early <- tryCatch(
            suppressWarnings(clock::date_time_parse(raw, zone = zone,
                format = "%Y-%m-%d %H:%M:%S", ambiguous = "earliest", nonexistent = "NA")),
            error = function(e) NULL
        )
        parsed_late <- tryCatch(
            suppressWarnings(clock::date_time_parse(raw, zone = zone,
                format = "%Y-%m-%d %H:%M:%S", ambiguous = "latest", nonexistent = "NA")),
            error = function(e) NULL
        )
        if (is.null(parsed_early) || is.null(parsed_late)) {
            issue[idx] <- "malformed_time"
            next
        }
        early <- as.POSIXct(parsed_early, tz = "UTC")
        late <- as.POSIXct(parsed_late, tz = "UTC")
        ambiguous_idx <- is.finite(as.numeric(early)) & is.finite(as.numeric(late)) & as.numeric(early) != as.numeric(late)
        local_idx <- which(ambiguous_idx)
        if (length(local_idx)) {
            issue[idx[local_idx]] <- if (ambiguous == "error") "ambiguous_local_time" else "ambiguous_time_resolved"
        }
        missing_epoch <- !is.finite(as.numeric(early)) & is.na(issue[idx])
        if (any(missing_epoch)) {
            probe <- tryCatch(
                suppressWarnings(clock::date_time_parse(raw[missing_epoch], zone = "UTC",
                    format = "%Y-%m-%d %H:%M:%S", ambiguous = "earliest", nonexistent = "NA")),
                error = function(e) NULL
            )
            valid_civil <- if (is.null(probe)) rep(FALSE, sum(missing_epoch)) else
                is.finite(as.numeric(as.POSIXct(probe, tz = "UTC")))
            issue[idx[which(missing_epoch)[valid_civil]]] <- "nonexistent_local_time"
            issue[idx[which(missing_epoch)[!valid_civil]]] <- "malformed_time"
        }
        selected <- if (ambiguous == "latest") late else early
        keep <- is.na(issue[idx]) | issue[idx] == "ambiguous_time_resolved"
        utc[idx[keep]] <- selected[keep]
        bad_parse <- !is.finite(as.numeric(early)) & is.na(issue[idx])
        issue[idx[bad_parse]] <- "malformed_time"
    }
    list(utc = utc, issue = issue, timezone = zones)
}

.empty_issues <- function() {
    S4Vectors::DataFrame(
        sample = character(), subject = character(), code = character(),
        severity = character(), message = character(), value = character()
    )
}

.issue_rows <- function(sample, subject, code, severity, message, value) {
    S4Vectors::DataFrame(
        sample = as.character(sample), subject = as.character(subject),
        code = as.character(code), severity = as.character(severity),
        message = as.character(message), value = as.character(value)
    )
}

.append_issues <- function(out, sample, subject, codes, severity = "error", value = NA_character_) {
    idx <- which(!is.na(codes))
    if (!length(idx)) return(out)
    msgs <- c(
        missing_time = "Sample time is missing.",
        malformed_time = "Sample time is not a supported elapsed-time token or ISO civil timestamp.",
        unit_required = "A time unit is required to interpret this numeric time.",
        non_finite_time = "Numeric sample time must be finite.",
        timezone_required = "A time zone is required to interpret a local civil timestamp.",
        invalid_timezone = "The supplied time zone is not a recognized IANA time zone.",
        ambiguous_local_time = "This local time occurs twice at a daylight-saving fall-back; disambiguate it explicitly.",
        ambiguous_time_resolved = "This ambiguous local time was resolved using the requested occurrence.",
        nonexistent_local_time = "This local time does not exist because of a daylight-saving spring-forward.",
        negative_elapsed_time = "Elapsed time is negative relative to the first sample for this subject.",
        duplicate_subject_time = "More than one sample has the same subject and temporal value.",
        non_monotonic_time = "Sample times are not increasing in their supplied order within this subject."
    )
    sev <- if (length(severity) == 1L) rep(severity, length(idx)) else severity[idx]
    val <- if (length(value) == 1L) rep(value, length(idx)) else value[idx]
    rbind(out, .issue_rows(sample[idx], subject[idx], codes[idx], sev,
        unname(msgs[codes[idx]]), val))
}

.check_new_columns <- function(cd, columns, overwrite) {
    exists <- intersect(columns, names(cd))
    if (length(exists) && !overwrite) {
        stop("Refusing to overwrite existing colData column(s): ",
             paste(exists, collapse = ", "), ". Set overwrite = TRUE to replace them.", call. = FALSE)
    }
}
