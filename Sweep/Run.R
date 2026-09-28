# Exhaustive comparison searches on the supplied BioModels network matrices.
# Run from the repository folder: Rscript --vanilla Sweep/Run.R
# Display the published table:   Rscript --vanilla Sweep/Run.R --summary

# Candidate matrices and parallel search ---------------------------------------

.fqCount <- new.env(parent = emptyenv())
.fqCount$n <- 0

.pairCache <- new.env(parent = emptyenv())
candidatePairs <- function(d) {
  key <- as.character(d)
  if (is.null(.pairCache[[key]])) {
    top <- 2^d - 1
    .pairCache[[key]] <- list(I = rep(top:0, times = (top + 1):1),
                              J = unlist(lapply(top:0, function(i) i:0), use.names = FALSE),
                              bits = 2^((d - 1):0))
  }
  .pairCache[[key]]
}

## The rows of E selected by the two masks; E holds e_1..e_d then -e_1..-e_d.
preorderFromPair <- function(i, j, d, E, bits) {
  rows <- c(which(bitwAnd(i, bits) != 0), d + which(bitwAnd(j, bits) != 0))
  E[rows, , drop = FALSE]
}

findOrderingsParExact <- function(source.complexes, product.complexes,
                                  ncores = max(1, detectCores() - 1),
                                  blocksPerCore = 4) {
  d <- ncol(source.complexes)
  RC <- populateRC(source.complexes, product.complexes)   # already exact
  E <- matrix(0, nrow = 2 * d, ncol = d)
  for (s in 1:d) { E[s, s] <- 1; E[d + s, s] <- -1 }

  pairs <- candidatePairs(d)
  ncand <- length(pairs$I)
  ncores <- max(1, min(as.integer(ncores), ncand))
  nblk <- max(1, min(ncand, as.integer(ncores * blocksPerCore)))
  bounds <- unique(round(seq(0, ncand, length.out = nblk + 1)))
  blocks <- Map(function(a, b) (a + 1):b, bounds[-length(bounds)], bounds[-1])

  # Combine the distinct orderings found in each block, up to reversal.
  runBlock <- function(idx) {
    .fqCount$n <- 0
    set <- newOrderingSet()
    for (k in idx) {
      M <- preorderFromPair(pairs$I[k], pairs$J[k], d, E, pairs$bits)
      chk <- checkConditions(source.complexes, RC$R, RC$C, M)
      if (chk$result == 1)
        addOrdering(set, list(rate.inequalities = chk$rate.inequalities,
                              species.inequalities = chk$species.inequalities))
    }
    list(out = as.list(set), nfq = .fqCount$n)
  }

  t0 <- proc.time()
  res <- if (ncores == 1) lapply(blocks, runBlock)
         else mclapply(blocks, runBlock, mc.cores = ncores, mc.preschedule = FALSE)
  el <- proc.time() - t0

  bad <- vapply(res, function(z) inherits(z, "try-error") || !is.list(z) ||
                                 is.null(z$nfq), logical(1))
  if (any(bad)) stop("findOrderingsParExact(): ", sum(bad), " worker block(s) failed")

  merged <- newOrderingSet()
  for (z in res) for (k in names(z$out)) if (!exists(k, envir = merged, inherits = FALSE))
    assign(k, z$out[[k]], envir = merged)
  orderings <- orderingSetToList(merged)
  nfq <- sum(vapply(res, `[[`, numeric(1), "nfq"))
  processed <- if (length(orderings) >= 2) processCanonical(orderings)
               else list(preorders = list(), equivalences = list())

  cpu <- el[["user.self"]] + el[["sys.self"]] + el[["user.child"]] + el[["sys.child"]]
  c(processed,
    list(candidates = ncand, lp = nfq, raw.hits = length(orderings),
         wall = el[["elapsed"]], cpu = cpu,
         cores = ncores, orderable = (length(processed$preorders) > 0 ||
                                      length(processed$equivalences) > 0)))
}


# Results and command-line options --------------------------------------------

resultColumns <- c('source', 'id', 'name', 'd', 'n', 'candidates', 'lp', 'wall_s',
                   'cpu_s', 'cores', 'orderable', 'n_preorders', 'n_equivalences',
                   'raw_hits', 'status')

readResults <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (!identical(names(x), resultColumns) || anyDuplicated(x$id) ||
      anyNA(x) || any(x$status != 'ok'))
    stop('Results must have the expected columns, unique IDs and only successful computations.')
  if (any(x$candidates != (4^x$d + 2^x$d) / 2) ||
      any(x$orderable != as.integer(x$n_preorders > 0 | x$n_equivalences > 0)))
    stop('Inconsistent candidate or successful-search counts.')
  x
}

summarizeResults <- function(path) {
  x <- readResults(path)
  if (!nrow(x)) { cat('No completed searches yet.\n'); return(invisible(NULL)) }
  # Round half milliseconds upwards, as in the manuscript's table.
  tab <- do.call(rbind, lapply(split(x, x$d), function(z) data.frame(
    d = z$d[1], Networks = nrow(z), Successful = sum(z$orderable),
    Matrices = z$candidates[1],
    Median = floor(median(round(z$wall_s * 1000)) + 0.5) / 1000,
    Maximum = max(z$wall_s))))
  tab <- tab[order(tab$d), ]
  tab$Median <- formatC(tab$Median, format = 'f', digits = 3)
  tab$Maximum <- formatC(tab$Maximum, format = 'f', digits = 3)
  print(tab, row.names = FALSE)
  cat(sprintf('\nSuccessful searches: %d/%d; matrices checked: %.0f.\n',
              sum(x$orderable), nrow(x), sum(x$candidates)))
  cat('Median and maximum times are in seconds and cover candidate checks only.\n')
  invisible(tab)
}

readOptions <- function(args, sweep.dir) {
  opt <- list(workers = 7L, max.species = 10L,
              output = file.path(sweep.dir, 'Results-new.csv'), summary = NULL, help = FALSE)
  if (!length(args)) return(opt)
  if (identical(args, '--help')) { opt$help <- TRUE; return(opt) }
  if (args[1] == '--summary') {
    if (length(args) > 2L) stop('Use --summary alone or followed by a CSV path.')
    opt$summary <- if (length(args) == 2L) args[2] else file.path(sweep.dir, 'Results.csv')
    return(opt)
  }
  if (length(args) %% 2L) stop('Each option needs a value. Use --help for instructions.')
  keys <- args[seq(1L, length(args), by = 2L)]
  if (anyDuplicated(keys)) stop('Each option may be specified only once.')
  for (i in seq(1L, length(args), by = 2L)) {
    key <- args[i]; value <- args[i + 1L]
    if (!key %in% c('--workers', '--max-species', '--output'))
      stop('Unknown option: ', key, '. Use --help for instructions.')
    if (key == '--output') {
      if (!nzchar(value)) stop('--output needs a file path.')
      opt$output <- value
    } else {
      number <- suppressWarnings(as.numeric(value))
      limit <- if (key == '--max-species') 10 else .Machine$integer.max
      if (length(number) != 1L || !is.finite(number) || number != floor(number) ||
          number < 1 || number > limit)
        stop(key, ' must be an integer between 1 and ', limit, '.')
      opt[[if (key == '--workers') 'workers' else 'max.species']] <- as.integer(number)
    }
  }
  opt
}

# Run the sweep ---------------------------------------------------------------

main <- function() {
  file.arg <- grep('^--file=', commandArgs(), value = TRUE)
  if (length(file.arg) != 1L) stop('Run this file with Rscript; see README.md.')
  sweep.dir <- dirname(normalizePath(sub('^--file=', '', file.arg), mustWork = TRUE))
  opt <- readOptions(commandArgs(trailingOnly = TRUE), sweep.dir)
  if (opt$help) {
    cat(paste(c('Usage: Rscript --vanilla Sweep/Run.R [options]',
                '  --summary [file.csv]    Summarize published results, or the specified file.',
                '  --workers N            Parallel workers for d >= 6 (default: 7).',
                '  --max-species D        Include networks with at most D species (default: 10).',
                '  --output file.csv      New results (default: Sweep/Results-new.csv).',
                '  --help                 Show these instructions.',
                '', 'New runs resume completed searches in the chosen output file.',
                'Published Results.csv is never used as a sweep output.',
                'Parallel execution requires macOS or Linux; Windows uses one worker.'),
              collapse = '\n'), '\n')
    return(invisible(NULL))
  }
  if (!is.null(opt$summary)) return(invisible(summarizeResults(opt$summary)))

  if (!dir.exists(dirname(opt$output))) stop('The output directory does not exist.')
  output <- if (file.exists(opt$output)) normalizePath(opt$output, mustWork = TRUE)
            else file.path(normalizePath(dirname(opt$output), mustWork = TRUE), basename(opt$output))
  published <- normalizePath(file.path(sweep.dir, 'Results.csv'), mustWork = TRUE)
  if (output == published) stop('Choose another output file: Results.csv contains the published measurements.')
  if (dir.exists(output)) stop('The output path must be a CSV file, not a directory.')
  if (.Platform$OS.type == 'windows' && opt$workers > 1L) {
    message('Using one worker on Windows.')
    opt$workers <- 1L
  }

  data <- new.env()
  sys.source(file.path(sweep.dir, 'Networks.R'), envir = data)
  work <- Filter(function(z) z$d <= opt$max.species, data$networks)
  previous <- if (file.exists(output)) readResults(output) else NULL
  done <- character(0)
  if (!is.null(previous)) {
    ids <- vapply(work, `[[`, '', 'id')
    for (i in seq_len(nrow(previous))) {
      j <- match(previous$id[i], ids)
      if (is.na(j)) stop('Existing results contain networks outside the selected sweep; choose a new output file.')
      z <- work[[j]]
      nc <- if (z$d >= 6L) opt$workers else 1L
      if (previous$source[i] != z$source || previous$name[i] != z$name ||
          previous$d[i] != z$d || previous$n[i] != z$n || previous$cores[i] != nc)
        stop('Existing results use different network data or workers; choose a new output file.')
    }
    done <- previous$id
    cat(sprintf('Resuming %d completed searches.\n', length(done)))
  }

  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(sys.source(file.path(dirname(sweep.dir), 'Functions.R'),
                                          envir = globalenv()))
  # Count cone computations without changing their mathematical tests.
  local({
    inner <- coneRow
    assign('coneRow', function(...) { .fqCount$n <- .fqCount$n + 1; inner(...) },
           envir = globalenv())
  })
  if (is.null(previous)) writeLines(paste(resultColumns, collapse = ','), output)
  cat(sprintf('Networks: %d; maximum species: %d; workers for d >= 6: %d.\n',
              length(work), opt$max.species, opt$workers))
  cat(sprintf('%-18s %3s %4s %10s %10s %10s\n', 'BioModels ID', 'd', 'n', 'Matrices', 'Seconds', 'Success'))
  started <- proc.time()[['elapsed']]
  for (z in work) {
    if (z$id %in% done) next
    nc <- if (z$d >= 6L) opt$workers else 1L
    r <- tryCatch(findOrderingsParExact(z$S, z$P, ncores = nc),
                  error = function(e) stop(z$id, ': ', conditionMessage(e),
                    '\nCompleted searches have been saved. Resolve the error and rerun to resume.', call. = FALSE))
    row <- data.frame(source = z$source, id = z$id, name = z$name, d = z$d, n = z$n,
                      candidates = r$candidates, lp = r$lp, wall_s = round(r$wall, 3),
                      cpu_s = round(r$cpu, 3), cores = r$cores,
                      orderable = as.integer(r$orderable), n_preorders = length(r$preorders),
                      n_equivalences = length(r$equivalences), raw_hits = r$raw.hits, status = 'ok')
    write.table(row, output, sep = ',', append = TRUE, col.names = FALSE,
                row.names = FALSE, qmethod = 'double')
    cat(sprintf('%-18s %3d %4d %10.0f %10.3f %10s\n', z$id, z$d, z$n,
                r$candidates, r$wall, if (r$orderable) 'yes' else 'no'))
    flush(stdout())
  }
  cat(sprintf('\nCompleted in %.1f seconds this invocation. Results: %s\n\n',
              proc.time()[['elapsed']] - started, output))
  summarizeResults(output)
}

if (sys.nframe() == 0L) main()
