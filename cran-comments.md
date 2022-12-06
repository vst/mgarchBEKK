# R CMD check results

Following commands are ran on the package:

```R
devtools::check(".", incoming = TRUE)
devtools::check_win_release(".")
devtools::check_win_devel(".")
```

No errors. No warnings. Please see below.

## devtools::check

1 note for spell check which is OK:

```txt
❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Vehbi Sinan Tunalioglu <vst@vsthost.com>’

  Possibly misspelled words in DESCRIPTION:
    BEKK (2:54, 10:18)
    GJR (10:40)
    MGARCH (2:46, 9:60)
    mGJR (2:63)

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

## devtools::check_win_release

No errors. No warnings. No notes.

## devtools::check_win_devel

No errors. No warnings. No notes.

## revdepcheck results

This is a release mainly for build chores to fix CRAN warnings. So, I did not
check revdeps.
