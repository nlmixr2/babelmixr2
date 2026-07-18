# Control for a PopED design task

Control for a PopED design task

## Usage

``` r
popedControl(
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  maxn = NULL,
  rxControl = NULL,
  sigdig = 4,
  important = NULL,
  unimportant = NULL,
  iFIMCalculationType = c("reduced", "full", "weighted", "loc", "reducedPFIM", "fullABC",
    "largeMat", "reducedFIMABC"),
  iApproximationMethod = c("fo", "foce", "focei", "foi"),
  iFOCENumInd = 1000,
  prior_fim = matrix(0, 0, 1),
  d_switch = c("d", "ed"),
  ofv_calc_type = c("lnD", "d", "a", "Ds", "inverse"),
  strEDPenaltyFile = "",
  ofv_fun = NULL,
  iEDCalculationType = c("mc", "laplace", "bfgs-laplace"),
  ED_samp_size = 45,
  bLHS = c("hypercube", "random"),
  bUseRandomSearch = TRUE,
  bUseStochasticGradient = TRUE,
  bUseLineSearch = TRUE,
  bUseExchangeAlgorithm = FALSE,
  bUseBFGSMinimizer = FALSE,
  bUseGrouped_xt = FALSE,
  EACriteria = c("modified", "fedorov"),
  strRunFile = "",
  poped_version = NULL,
  modtit = "PopED babelmixr2 model",
  output_file = "PopED_output_summary",
  output_function_file = "PopED_output_",
  strIterationFileName = "PopED_current.R",
  user_data = NULL,
  ourzero = 1e-05,
  dSeed = NULL,
  line_opta = NULL,
  line_optx = NULL,
  bShowGraphs = FALSE,
  use_logfile = FALSE,
  m1_switch = c("central", "complex", "analytic", "ad"),
  m2_switch = c("central", "complex", "analytic", "ad"),
  hle_switch = c("central", "complex", "ad"),
  gradff_switch = c("central", "complex", "analytic", "ad"),
  gradfg_switch = c("central", "complex", "analytic", "ad"),
  grad_all_switch = c("central", "complex"),
  rsit_output = 5,
  sgit_output = 1,
  hm1 = 1e-05,
  hlf = 1e-05,
  hlg = 1e-05,
  hm2 = 1e-05,
  hgd = 1e-05,
  hle = 1e-05,
  AbsTol = 1e-06,
  RelTol = 1e-06,
  iDiffSolverMethod = NULL,
  bUseMemorySolver = FALSE,
  rsit = 300,
  sgit = 150,
  intrsit = 250,
  intsgit = 50,
  maxrsnullit = 50,
  convergence_eps = 1e-08,
  rslxt = 10,
  rsla = 10,
  cfaxt = 0.001,
  cfaa = 0.001,
  bGreedyGroupOpt = FALSE,
  EAStepSize = 0.01,
  EANumPoints = FALSE,
  EAConvergenceCriteria = 1e-20,
  bEANoReplicates = FALSE,
  BFGSProjectedGradientTol = 1e-04,
  BFGSTolerancef = 0.001,
  BFGSToleranceg = 0.9,
  BFGSTolerancex = 0.1,
  ED_diff_it = 30,
  ED_diff_percent = 10,
  line_search_it = 50,
  Doptim_iter = 1,
  iCompileOption = c("none", "full", "mcc", "mpi"),
  compileOnly = FALSE,
  iUseParallelMethod = c("mpi", "matlab"),
  MCC_Dep = NULL,
  strExecuteName = "calc_fim.exe",
  iNumProcesses = 2,
  iNumChunkDesignEvals = -2,
  Mat_Out_Pre = "parallel_output",
  strExtraRunOptions = "",
  dPollResultTime = 0.1,
  strFunctionInputName = "function_input",
  bParallelRS = FALSE,
  bParallelSG = FALSE,
  bParallelMFEA = FALSE,
  bParallelLS = FALSE,
  groupsize = NULL,
  time = "time",
  timeLow = "low",
  timeHi = "high",
  id = "id",
  m = NULL,
  x = NULL,
  ni = NULL,
  maxni = NULL,
  minni = NULL,
  maxtotni = NULL,
  mintotni = NULL,
  maxgroupsize = NULL,
  mingroupsize = NULL,
  maxtotgroupsize = NULL,
  mintotgroupsize = NULL,
  xt_space = NULL,
  a = NULL,
  maxa = NULL,
  mina = NULL,
  a_space = NULL,
  x_space = NULL,
  use_grouped_xt = FALSE,
  grouped_xt = NULL,
  use_grouped_a = FALSE,
  grouped_a = NULL,
  use_grouped_x = FALSE,
  grouped_x = NULL,
  our_zero = NULL,
  auto_pointer = "",
  user_distribution_pointer = "",
  minxt = NULL,
  maxxt = NULL,
  discrete_xt = NULL,
  discrete_a = NULL,
  fixRes = FALSE,
  script = NULL,
  overwrite = TRUE,
  literalFix = TRUE,
  opt_xt = FALSE,
  opt_a = FALSE,
  opt_x = FALSE,
  opt_samps = FALSE,
  optTime = TRUE,
  literalFixRes = FALSE,
  ...
)
```

## Arguments

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

- maxn:

  Maximum number of design points for optimization; By default this is
  declared by the maximum number of design points in the babelmixr2
  dataset (when `NULL`)

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- sigdig:

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), ODE solver tolerance (`0.5*10^(-sigdig-2)`,
  or `0.5*10^(-sigdig-1.5)` for sensitivity/steady-state with liblsoda),
  and boundary check tolerance (`5*10^(-sigdig+1)`).

- important:

  character vector of important parameters or NULL for default. This is
  used with Ds-optimality

- unimportant:

  character vector of unimportant parameters or NULL for default. This
  is used with Ds-optimality

- iFIMCalculationType:

  can be either an integer or a named value of the Fisher Information
  Matrix type:

  - 0/"full" = Full FIM

  - 1/"reduced" = Reduced FIM

  - 2/"weighted" = weighted models

  - 3/"loc" = Loc models

  - 4/"reducedPFIM" = reduced FIM with derivative of SD of sigma as in
    PFIM

  - 5/"fullABC" = FULL FIM parameterized with A,B,C matrices &
    derivative of variance

  - 6/"largeMat" = Calculate one model switch at a time, good for large
    matrices

  - 7/"reducedFIMABC" = =Reduced FIM parameterized with A,B,C matrices &
    derivative of variance

- iApproximationMethod:

  Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI

- iFOCENumInd:

  integer; number of individuals in focei solve

- prior_fim:

  matrix; prior FIM

- d_switch:

  integer or character option:

  - 0/"ed" = ED design

  - 1/"d" = D design

- ofv_calc_type:

  objective calculation type:

  - 1/"d" = D-optimality". Determinant of the FIM: det(FIM)

  - 2/"a" = "A-optimality". Inverse of the sum of the expected parameter
    variances: 1/trace_matrix(inv(FIM))

  - 4/"lnD" = "lnD-optimality". Natural logarithm of the determinant of
    the FIM: log(det(FIM))

  - 6/"Ds" = "Ds-optimality". Ratio of the Determinant of the FIM and
    the Determinant of the uninteresting rows and columns of the FIM:
    det(FIM)/det(FIM_u)

  - 7/"inverse" = Inverse of the sum of the expected parameter RSE:
    1/sum(get_rse(FIM,poped.db,use_percent=FALSE))

- strEDPenaltyFile:

  Penalty function name or path and filename, empty string means no
  penalty. User defined criterion can be defined this way.

- ofv_fun:

  User defined function used to compute the objective function. The
  function must have a poped database object as its first argument and
  have "..." in its argument list. Can be referenced as a function or as
  a file name where the function defined in the file has the same name
  as the file. e.g. "cost.txt" has a function named "cost" in it.

- iEDCalculationType:

  ED Integral Calculation type:

  - 0/"mc" = Monte-Carlo-Integration

  - 1/"laplace" = Laplace Approximation

  - 2/"bfgs-laplace" = BFGS Laplace Approximation

- ED_samp_size:

  Sample size for E-family sampling

- bLHS:

  How to sample from distributions in E-family calculations. 0=Random
  Sampling, 1=LatinHyperCube –

- bUseRandomSearch:

  - **\*\*\*\*\*\*START OF Optimization algorithm SPECIFICATION
    OPTIONS\*\*\*\*\*\*\*\*\*\***

  Use random search (1=TRUE, 0=FALSE)

- bUseStochasticGradient:

  Use Stochastic Gradient search (1=TRUE, 0=FALSE)

- bUseLineSearch:

  Use Line search (1=TRUE, 0=FALSE)

- bUseExchangeAlgorithm:

  Use Exchange algorithm (1=TRUE, 0=FALSE)

- bUseBFGSMinimizer:

  Use BFGS Minimizer (1=TRUE, 0=FALSE)

- bUseGrouped_xt:

  Use grouped time points (1=TRUE, 0=FALSE).

- EACriteria:

  Exchange Algorithm Criteria:

  - 1/"modified" = Modified

  - 2/"fedorov" = Fedorov

- strRunFile:

  Filename and path, or function name, for a run file that is used
  instead of the regular PopED call.

- poped_version:

  - **\*\*\*\*\*\*START OF Labeling and file names SPECIFICATION
    OPTIONS\*\*\*\*\*\*\*\*\*\***

  The current PopED version

- modtit:

  The model title

- output_file:

  Filename and path of the output file during search

- output_function_file:

  Filename suffix of the result function file

- strIterationFileName:

  Filename and path for storage of current optimal design

- user_data:

  - **\*\*\*\*\*\*START OF Miscellaneous SPECIFICATION
    OPTIONS\*\*\*\*\*\*\*\*\*\***

  User defined data structure that, for example could be used to send in
  data to the model

- ourzero:

  Value to interpret as zero in design

- dSeed:

  The seed number used for optimization and sampling – integer or -1
  which creates a random seed `as.integer(Sys.time())` or NULL.

- line_opta:

  Vector for line search on continuous design variables (1=TRUE,0=FALSE)

- line_optx:

  Vector for line search on discrete design variables (1=TRUE,0=FALSE)

- bShowGraphs:

  Use graph output during search

- use_logfile:

  If a log file should be used (0=FALSE, 1=TRUE)

- m1_switch:

  Method used to calculate M1:

  - 1/"central" = Central difference

  - 0/"complex" = Complex difference

  - 20/"analytic" = Analytic derivative

  - 30/"ad" = Automatic differentiation

- m2_switch:

  Method used to calculate M2:

  - 1/"central" = Central difference

  - 0/"complex" = Complex difference

  - 20/"analytic" = Analytic derivative

  - 30/"ad" = Automatic differentiation

- hle_switch:

  Method used to calculate linearization of residual error:

  - 1/"central" = Central difference

  - 0/"complex" = Complex difference

  - 30/"ad" = Automatic differentiation

- gradff_switch:

  Method used to calculate the gradient of the model:

  - 1/"central" = Central difference

  - 0/"complex" = Complex difference

  - 20/"analytic" = Analytic derivative

  - 30/"ad" = Automatic differentiation

- gradfg_switch:

  Method used to calculate the gradient of the parameter vector g:

  - 1/"central" = Central difference

  - 0/"complex" = Complex difference

  - 20/"analytic" = Analytic derivative

  - 30/"ad" = Automatic differentiation

- grad_all_switch:

  Method used to calculate all the gradients:

  - 1/"central" = Central difference

  - 0/"complex" = Complex difference

- rsit_output:

  Number of iterations in random search between screen output

- sgit_output:

  Number of iterations in stochastic gradient search between screen
  output

- hm1:

  Step length of derivative of linearized model w.r.t. typical values

- hlf:

  Step length of derivative of model w.r.t. g

- hlg:

  Step length of derivative of g w.r.t. b

- hm2:

  Step length of derivative of variance w.r.t. typical values

- hgd:

  Step length of derivative of OFV w.r.t. time

- hle:

  Step length of derivative of model w.r.t. sigma

- AbsTol:

  The absolute tolerance for the diff equation solver

- RelTol:

  The relative tolerance for the diff equation solver

- iDiffSolverMethod:

  The diff equation solver method, NULL as default.

- bUseMemorySolver:

  If the differential equation results should be stored in memory (1) or
  not (0)

- rsit:

  Number of Random search iterations

- sgit:

  Number of stochastic gradient iterations

- intrsit:

  Number of Random search iterations with discrete optimization.

- intsgit:

  Number of Stochastic Gradient search iterations with discrete
  optimization

- maxrsnullit:

  Iterations until adaptive narrowing in random search

- convergence_eps:

  Stochastic Gradient convergence value, (difference in OFV for
  D-optimal, difference in gradient for ED-optimal)

- rslxt:

  Random search locality factor for sample times

- rsla:

  Random search locality factor for covariates

- cfaxt:

  Stochastic Gradient search first step factor for sample times

- cfaa:

  Stochastic Gradient search first step factor for covariates

- bGreedyGroupOpt:

  Use greedy algorithm for group assignment optimization

- EAStepSize:

  Exchange Algorithm StepSize

- EANumPoints:

  Exchange Algorithm NumPoints

- EAConvergenceCriteria:

  Exchange Algorithm Convergence Limit/Criteria

- bEANoReplicates:

  Avoid replicate samples when using Exchange Algorithm

- BFGSProjectedGradientTol:

  BFGS Minimizer Convergence Criteria Normalized Projected Gradient
  Tolerance

- BFGSTolerancef:

  BFGS Minimizer Line Search Tolerance f

- BFGSToleranceg:

  BFGS Minimizer Line Search Tolerance g

- BFGSTolerancex:

  BFGS Minimizer Line Search Tolerance x

- ED_diff_it:

  Number of iterations in ED-optimal design to calculate convergence
  criteria

- ED_diff_percent:

  ED-optimal design convergence criteria in percent

- line_search_it:

  Number of grid points in the line search

- Doptim_iter:

  Number of iterations of full Random search and full Stochastic
  Gradient if line search is not used

- iCompileOption:

  Compile options for PopED

  - "none"/-1 = No compilation

  - "full/0 or 3 = Full compilation

  - "mcc"/1 or 4 = Only using MCC (shared lib)

  - "mpi"/2 or 5 = Only MPI,

  When using numbers, option 0,1,2 runs PopED and option 3,4,5 stops
  after compilation.

  When using characters, the option `compileOnly` determines if the
  model is only compiled (and PopED is not run).

- compileOnly:

  logical; only compile the model, do not run PopED (in conjunction with
  `iCompileOption`)

- iUseParallelMethod:

  Parallel method to use

  - 0/"matlab"= Matlab PCT

  - 1/"mpi" = MPI

- MCC_Dep:

  Additional dependencies used in MCC compilation (mat-files), if
  several space separated

- strExecuteName:

  Compilation output executable name

- iNumProcesses:

  Number of processes to use when running in parallel (e.g. 3 = 2
  workers, 1 job manager)

- iNumChunkDesignEvals:

  Number of design evaluations that should be evaluated in each process
  before getting new work from job manager

- Mat_Out_Pre:

  The prefix of the output mat file to communicate with the executable

- strExtraRunOptions:

  Extra options send to e\$g. the MPI executable or a batch script, see
  execute_parallel\$m for more information and options

- dPollResultTime:

  Polling time to check if the parallel execution is finished

- strFunctionInputName:

  The file containing the popedInput structure that should be used to
  evaluate the designs

- bParallelRS:

  If the random search is going to be executed in parallel

- bParallelSG:

  If the stochastic gradient search is going to be executed in parallel

- bParallelMFEA:

  If the modified exchange algorithm is going to be executed in parallel

- bParallelLS:

  If the line search is going to be executed in parallel

- groupsize:

  Vector defining the size of the different groups (num individuals in
  each group). If only one number then the number will be the same in
  every group.

- time:

  string that represents the time in the dataset (ie xt)

- timeLow:

  string that represents the lower design time (ie minxt)

- timeHi:

  string that represents the upper design time (ie maxmt)

- id:

  The id variable

- m:

  Number of groups in the study. Each individual in a group will have
  the same design.

- x:

  A matrix defining the initial discrete values for the model Each row
  is a group/individual.

- ni:

  Vector defining the number of samples for each group.

- maxni:

  - **\*\*\*\*\*\*START OF DESIGN SPACE OPTIONS\*\*\*\*\*\*\*\*\*\***

  Max number of samples per group/individual

- minni:

  Min number of samples per group/individual

- maxtotni:

  Number defining the maximum number of samples allowed in the
  experiment.

- mintotni:

  Number defining the minimum number of samples allowed in the
  experiment.

- maxgroupsize:

  Vector defining the max size of the different groups (max number of
  individuals in each group)

- mingroupsize:

  Vector defining the min size of the different groups (min num
  individuals in each group) –

- maxtotgroupsize:

  The total maximal groupsize over all groups

- mintotgroupsize:

  The total minimal groupsize over all groups

- xt_space:

  Cell array
  [`cell`](https://andrewhooker.github.io/PopED/reference/cell.html)
  defining the discrete variables allowed for each xt value. Can also be
  a vector of values `c(1:10)` (same values allowed for all xt), or a
  list of lists `list(1:10, 2:23, 4:6)` (one for each value in xt in row
  major order or just for one row in xt, and all other rows will be
  duplicated).

- a:

  Matrix defining the initial continuous covariate values. n_rows=number
  of groups, n_cols=number of covariates. If the number of rows is one
  and the number of groups \> 1 then all groups are assigned the same
  values.

- maxa:

  Vector defining the max value for each covariate. If a single value is
  supplied then all a values are given the same max value

- mina:

  Vector defining the min value for each covariate. If a single value is
  supplied then all a values are given the same max value

- a_space:

  Cell array
  [`cell`](https://andrewhooker.github.io/PopED/reference/cell.html)
  defining the discrete variables allowed for each a value. Can also be
  a list of values `list(1:10)` (same values allowed for all a), or a
  list of lists `list(1:10, 2:23, 4:6)` (one for each value in a).

- x_space:

  Cell array
  [`cell`](https://andrewhooker.github.io/PopED/reference/cell.html)
  defining the discrete variables for each x value.

- use_grouped_xt:

  Group sampling times between groups so that each group has the same
  values (`TRUE` or `FALSE`).

- grouped_xt:

  Matrix defining the grouping of sample points. Matching integers mean
  that the points are matched. Allows for finer control than
  `use_grouped_xt`

- use_grouped_a:

  Group continuous design variables between groups so that each group
  has the same values (`TRUE` or `FALSE`).

- grouped_a:

  Matrix defining the grouping of continuous design variables. Matching
  integers mean that the values are matched. Allows for finer control
  than `use_grouped_a`.

- use_grouped_x:

  Group discrete design variables between groups so that each group has
  the same values (`TRUE` or `FALSE`).

- grouped_x:

  Matrix defining the grouping of discrete design variables. Matching
  integers mean that the values are matched. Allows for finer control
  than `use_grouped_x`.

- our_zero:

  Value to interpret as zero in design.

- auto_pointer:

  Filename and path, or function name, for the Autocorrelation function,
  empty string means no autocorrelation

- user_distribution_pointer:

  Filename and path, or function name, for user defined distributions
  for E-family designs

- minxt:

  Matrix or single value defining the minimum value for each xt sample.
  If a single value is supplied then all xt values are given the same
  minimum value

- maxxt:

  Matrix or single value defining the maximum value for each xt sample.
  If a single value is supplied then all xt values are given the same
  maximum value.

- discrete_xt:

  Cell array
  [`cell`](https://andrewhooker.github.io/PopED/reference/cell.html)
  defining the discrete variables allowed for each xt value. Can also be
  a list of values `list(1:10)` (same values allowed for all xt), or a
  list of lists `list(1:10, 2:23, 4:6)` (one for each value in xt). See
  examples in
  [`create_design_space`](https://andrewhooker.github.io/PopED/reference/create_design_space.html).

- discrete_a:

  Cell array
  [`cell`](https://andrewhooker.github.io/PopED/reference/cell.html)
  defining the discrete variables allowed for each a value. Can also be
  a list of values `list(1:10)` (same values allowed for all a), or a
  list of lists `list(1:10, 2:23, 4:6)` (one for each value in a). See
  examples in
  [`create_design_space`](https://andrewhooker.github.io/PopED/reference/create_design_space.html).

- fixRes:

  boolean; Fix the residuals to what is specified by the model

- script:

  write a PopED/rxode2 script that can be modified for more fine
  control. The default is NULL.

  When `script` is TRUE, the script is returned as a lines that would be
  written to a file and with the class `babelmixr2popedScript`. This
  allows it to be printed as the script on screen.

  When `script` is a file name (with an R extension), the script is
  written to that file.

- overwrite:

  \[`logical(1)`\]  
  If `TRUE`, an existing file in place is allowed if it it is both
  readable and writable. Default is `FALSE`.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- opt_xt:

  boolean to indicate if this is meant for optimizing times

- opt_a:

  boolean to indicate if this is meant for optimizing covariates

- opt_x:

  boolean to indicate if the discrete design variables be optimized

- opt_samps:

  boolean to indicate if the sample optimizer is used (not implemented
  yet in `PopED`)

- optTime:

  boolean to indicate if the global time indexer inside of babelmixr2 is
  reset if the times are different. By default this is `TRUE`. If
  `FALSE` you can get slightly better run times and possibly slightly
  different results. When `optTime` is `FALSE` the global indexer is
  reset every time the PopED rxode2 is setup for a problem or when a
  poped dataset is created. You can manually reset with
  [`popedMultipleEndpointResetTimeIndex()`](https://nlmixr2.github.io/babelmixr2/reference/popedMultipleEndpointResetTimeIndex.md)

- literalFixRes:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- ...:

  other parameters for PopED control

## Value

popedControl object

## Author

Matthew L. Fidler
