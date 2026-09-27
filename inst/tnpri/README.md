# babelmixr2 NONMEM TNPRI test kit

This kit checks babelmixr2's support for NONMEM's `$PRIOR TNPRI`
(`nonmemControl(tnpri=nonmemTnpri(...))`, issue #206) with a real NONMEM
installation. It is self-contained: it needs no internet access and no
AI tools. It runs a few NONMEM fits, and collects every control stream,
data file and NONMEM output into **one archive**, plus a text summary
that has everything needed if only text can be copied off the machine.

## What it checks

A TNPRI prior is an earlier NONMEM fit of the same model. babelmixr2
fits the model to the prior data (with `$MSFO` and `$COVARIANCE`), and
then writes a control stream with two problems. Problem 1 reads the
model specification file (`$MSFI ... ONLYREAD`) with `$PRIOR TNPRI
(PROBLEM 2)` and holds the model code. Problem 2 fits the new data
(`$DATA ... REWIND`).

The kit uses `nlmixr2data::Oral_1CPT`: subjects 1-60 are the prior
study and subjects 61-120 the new one.

| case | what it runs |
|---|---|
| `reference` | the model fit to the new data **without** a prior, to compare objective functions with |
| `tnpri` | babelmixr2 end to end: prior fit (dataset prior), then the TNPRI fit, read back into nlmixr2 |
| `tnpri-lincmt` | the same with `linCmt()` (closed form ADVAN when available), an omega block, combined error and `MODE=1` |
| `tnpri-fit` | the prior given as an nlmixr2 `focei` fit (babelmixr2 refits its data with NONMEM) |
| `tnpri-imp` | `est="imp"` with TNPRI must be refused before anything is written or run (NONMEM's help: do not use TNPRI with the NONMEM 7 methods) |
| `variant-no-plev` | the `tnpri` control stream without `PLEV=0` (NONMEM's own default) |
| `variant-no-input2` | ... without `$INPUT` in problem 2 |
| `variant-code2` | ... with the model code repeated in problem 2 |
| `variant-no-code1` | ... with the model code only in problem 2 |
| `variant-extra-theta` | ... with one more `$THETA` in problem 2 than in the prior |

The `variant-*` cases run NONMEM directly on edited copies of the
`tnpri` control stream, to learn which forms NONMEM accepts. Some of
them are *expected* to fail; that is the information they give. They
reuse the prior's model specification file, so they need the `tnpri`
case.

About 10 NONMEM runs in total; each is a small FOCEI fit with the
covariance step.

## Requirements

- R with the babelmixr2 version to test and its dependencies
  (nlmixr2est, rxode2, nonmem2rx, nlmixr2data, withr). The TNPRI support
  is on the `issue-206` branch until it is merged. On a machine with
  internet access, build the package with `R CMD build babelmixr2`, copy
  the `babelmixr2_*.tar.gz` file over, and install it with
  `R CMD INSTALL babelmixr2_*.tar.gz`.
- NONMEM and its run command (like `nmfe75`), either on the `PATH` or
  with its full path.

## Running it

Find the script:

```sh
KIT=$(Rscript -e 'cat(system.file("tnpri", "run-tnpri.R", package="babelmixr2"))')
```

Run every case:

```sh
Rscript "$KIT" --nonmem=nmfe75
# or with the full path
Rscript "$KIT" --nonmem=/opt/nm75/run/nmfe75
```

Other options:

```sh
Rscript "$KIT" --list                  # list the cases
Rscript "$KIT" --generate              # only write the files, no NONMEM
Rscript "$KIT" --nonmem=nmfe75 --cases='^tnpri$|variant'
Rscript "$KIT" --nonmem=nmfe75 --out=my-tnpri-results
```

## What to send back

The run ends with two paths:

- `babelmixr2-tnpri-<date>.tar.gz`: everything (preferred).
- `babelmixr2-tnpri-<date>/summary.md`: a text-only report. It has the
  status and objective function of every case, the key lines of every
  NONMEM listing (NM-TRAN errors, termination messages, `#OBJV`), what
  nonmem2rx reads from each output file, and the estimates babelmixr2
  read back. If files cannot leave the machine, this is the one to copy.

The output directory has one folder per case under `cases/`, with:

- `case.txt`: the result babelmixr2 read back (or its error)
- `*-nonmem/`: the control stream, data and NONMEM output of each run
  (the prior run is in `*_prior-nonmem/`)
- `*-lst-summary.txt`: the key lines of the NONMEM listing
- `*-readers.txt`: what nonmem2rx reads from the `.lst`, `.ext`, `.cov`
  and eta table files
