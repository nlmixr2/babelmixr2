# babelmixr2 0.1.5

* Reduced number of items run on CRAN (reduce running time)

* Fix bug where `PopED` could error with certain `dvid` values

* Fix bug where if/else clauses in the model could cause the model to
  not predict the values correctly.

* Fix bug so that `shrinkage()` calculation works

* Fix bug so that you can mix 2 different `PopED` data bases in an
  analysis without crashing R.  While this didn't occur with every
  database clash, it more frequently occurred when you interleaved
  `PopED` code between two different `PopED` databases, like in issue
  #131.

* Added a new function `babelBpopIdx(poped.db, "par")` which will get
  the poped index for a model generated from `babelmixr2`, which is
  useful when calculating the power (as in example 11).
