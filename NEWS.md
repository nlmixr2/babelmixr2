# babelmixr2 0.1.2

* Handle algebraic `mu` expressions

* PKNCA controller now contains `rxControl` since it is used for some
  translation options

* This revision will load the pruned ui model to query the compartment
  properties (i.e. bioavailability, lag time, etc) when writing out the
  NONMEM model.  It should fix issues where the PK block does not
  define some of the variables and will have a larger calculated
  variable that can be used in the model instead.

* When `nonmem2rx` has a different `lst` file, as long as
  `nonmem2rx::nminfo(file)` works, then a successful conversion to a
  `nlmixr2` fit object will occur.

* Fix to save parameter history into `$parHistData` to accommodate
  changes in `focei`'s output (`$parHist` is now derived).

* Changed the solving options to match the new steady state options in
  `rxode2` and how NONMEM implements them.  Also changed the iwres
  model to account for the `rxerr.` instead of the `err.` which was
  updated in `rxode2` as well.


# babelmixr2 0.1.1

* Add new method `as.nlmixr2` to convert `nonmem2rx` methods to `nlmixr` fits

* Dropped `pmxTools` in favor of `nonmem2rx` to conserve some of the
  methods

# babelmixr2 0.1.0

* Babelmixr has support for "monolix", "nonmem", and "pknca" methods
  on release.

* Added a `NEWS.md` file to track changes to the package.
