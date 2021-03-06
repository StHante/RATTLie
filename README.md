# `RATTLie`: An integrator for mechanical systems with holonomic constraints on Lie group structured configuration spaces

This is a time integration method that can be used to solve differential algebraic equations on of the form

        Μπ(π‘) = dπΏ_π(π) Μπ£(π‘)
    π β Μπ£(π‘) = Μπ£α΅βπβπ£ - π(π‘,π,π£) - πα΅(π(π‘))βΞ»(π‘)
           π = Ξ¦(π(π‘))

with a Lie group πΊ of dimension π. Here π(π‘)βπΊ is the configuration of the system and π£(π‘)βββΏ the associated velocity vectors. πΏ is the left translation and dπΏ its derivative. We have the tilde operator which maps π-vectors to the Lie algebra of πΊ. The πΓπ-matrix π is a constant mass matrix, π is a function that gives the negative of all generalized forces excluding Coriolis (or inertial) forces, which are explicitly given as the term Μπ£α΅βπβπ£. Moreover, Ξ¦ is the constraint function with Ξ¦(π)ββα΅ and π is its derivative in the sense
    π(π)βπ€ = dΞ¦(π) dπΏ_π(π) Μπ€  for all  π€βββΏ.
Using the derivative operator π, we could also write in short π=πΞ¦. Furthermore, Ξ»(π‘)ββα΅ are the Lagrange multipliers.
RATTLie is a generalization of the RATTLE integration scheme to Lie group structured configuration spaces.

Alternatively, it can also solve unconstrained equations of motion of the form

        Μπ(π‘) = dπΏ_π(π) Μπ£(π‘)
    π β Μπ£(π‘) = Μπ£α΅βπβπ£ - π(π‘,π,π£)

In this case, RATTLie becomes a generalization of the StΓΆrmer-Verlet integration scheme to Lie group structured configuration spaces.

The integrator `RATTLie` was implemented in modern Fortran for Linux. (Compiling it on Windows, however should not be too hard.)

## Limitations
There are some limitations for using RATTLie, that might be improved in the future:

 * In the constrained case, no good approximation for the Lagrange multipliers at the starting time π‘β is calculated.
 * Mass matrices that actually depend on the configuration π have not been tested.
 * Problems, where the generalized forces π depend on the velocity vectors π£ are not supported. Currently, the calculatation of the velocity vectors assumes that d/dπ£ π(π‘,π,π£) = 0 and instead of a full Newton method, only one Newton step is executed (which is sufficient to solve the system, because it is linear under the condition that π is independent of π£).

## Content
`RATTLie` is a Fortran module that defines the following types:

 * The abstract type `RATTLie_problem`. When implementing a problem, this type needs to be extended to a non-abstract type by implementing all deferred procedures.
 * The type `RATTLie_options`, which contains parameters of the integration scheme.
 * The type `RATTLie_statistics`, which will contain information about the runtime after a performed integration step.

In the following, `prob` or `this` will be the variable of class `RATTLie_problem`.

### Deferred procedures to be implemented
Here is a list of all deferred procedures that must be implemented when extending `RATTLie_problem` including a short description. Note that it is helpful if a procedure that is never called or a dummy, it should contain `error stop "<this procedure> should never be called"`.

 * `RATTLie_M`: This function should return the mass matrix π in dependence of the configuration π. Note however, that mass matrices that are not constant have _not_ been tested and the results might not be correct. This function will never be called, if `prob%opts%diag_mass_matrix==1`.
  * `RATTLie_diag_M`: This function should return the diagonal of the mass matrix π as a vector. (See also `RATTLie_M`) This function will never be called, if `prob%opts%diag_mass_matrix==0`.
 * `RATTLie_f`: This function should return π(π‘,π,π£), the negative of of all applied generalized forces. Note that the result must exclude any Coriolis forces (inertial forces) as well as, of course, any constraint forces.
 * `RATTLie_inertial`: This function should return the generalized Coriolis forces Μπ£βπβπ£.
 * `RATTLie_qlpexpDqtilde`: The cryptic name means "Q Lie Product EXP of Delta Q TILDE". It should return the result of πβexp(ββtilde(Ξπ)), where β is the Lie group product, exp:πβπΊβπΊ the exponential map and tilde:ββΏβπβπΊ the tilde operator. We could also write this as πβexpt(ββΞπ).
 * `RATTLie_itlbtvtw`: The cryptic name means "Inverse Tilde of Lie Bracket of Tilde V and Tilde W". This is a very complicated way of saying that this function should return Μπ£βπ€, where Μπ£ is the application of the hat operator to the vector π£. This function will never be called if the system is unconstrained, ie. `prob%opts%constrained==0`.
 * `RATTLie_Ct`: This functions is currently a dummy function and is never called.
 * `RATTLie_Kt`: This functions is currently a dummy function and is never called.
 * `RATTLie_Kt_lambda`: This functions is currently a dummy function and is never called.
 * `RATTLie_Tg`: This functions should return the tangent operator π(ββΞπ). The tangent operator should be implemented in such a way that it works reliably, even if the norm of ββΞπ is very small. Note that for backwards compatibility reasons, this function takes two arguments: a scalar β and a vector Ξπ, although only their product ββΞπ will be used.
 * `RATTLie_Tg_inv_T`: This function should return the transpose of the inverse matrix πβ»α΅(π€) = ((π(π€))β»ΒΉ)α΅ of the tangent operator. Note that unlike `RATTLie_Tg` this function only takes one vector argument π€.
 * `RATTLie_d_Tg_inv_T`: This function should return the Jacobi matrix of πβ»α΅(π£)βπ€ with respect to π£.
 * `RATTLie_outputFunction`: Despite its name, this is not a function but rather a subroutine. It will be called after each successful integration step. Note, that all arguments are `intent(in)`, meaning that the problem object may not be altered. There are currently three possible values of `info`:
   * `info==0`: "Initialization": The `RATTLie_outputFunction` is called once after the problem was initialized (ie. after `RATTLie_init`). Note that here, we should not open files, because there is no way of storing the identifier. Open files in `RATTLie_init` or before calling the integration routine.
   * `info==1`: "Normal output": The `RATTLie_outputFunction` was called after a successful integration step.
   * `info==99`: "Termination": The `RATTLie_outputFunction` is called once at the very end of the full integration run.
 * `RATTLie_init`: This subroutine should initialize the problem. There are a few things that must happen in this subroutine:
   * Set the sizes `this%sizeq`, `this%sizev` and in the constrained case also `this%sizel`.
   * Allocate `this%q`, `this%v`, `this%p` and in the constrained case also `this%l`, `this%lm` and `this%lp`. Note that it is a good idea to check, whether the variables are already allocated and if so, deallocate them before allocating.
   * Set initial values:
     * `this%t = this%opts%t0`
     * `this%q`
     * `this%v`
   * Allocate and fill `this%opts%jour` if `this%opts%banded_iteration_matrix==1`.
 * `RATTLie_Phi`: This function should return the constraint function Ξ¦(π). It will never be called in the unconstrained case (ie. `prob%opts%constrained==0`).
 * `RATTLie_B`: This function should return the derivative of the constraint function π(π)=πΞ¦(π). In other words we have π(π)βπ€ = d/dπ Ξ¦(π) dπΏ_π(π) Μπ€ for all vectors π€. This function will never be called in the unconstrained case (ie. `prob%opts%constrained==0`).
 * `RATTLie_Z`: This function should return the curvature term π_π(π(π)βπ£)βπ£=d_π(π(π)βπ£) dπΏ_π(π) Μπ£. This function will never be called in the unconstrained case (ie. `prob%opts%constrained==0`) and only is used in the case that the preprocessor variable `USE_INDEX_1` is defined.
 * `RATTLie_matZ`: This functions is currently a dummy function and is never called.

### Public subroutines
The following subroutines are meant to be called by the problem file.

 * `RATTLie_integrate`: Call this subroutine in order to start the integration process after all important variables and integrator options have been set. This routine will call `RATTLie_init` once and `RATTLie_outputFunction` after each successful integration step. Note that `RATTLie` does not take care of saving intermediate results. This has to be done by the `RATTLie_outputFunction`.
 * `RATTLie_print_stats`: Prints the contents of `prob%RATTLie_stats` to the standard output.
 * `RATTLie_cleanup`: Resets most of the internal variables and integrator options and deallocates most internal allocatable variables. This will, of course, not reset any variables that are added to `RATTLie_problem` by its extension.

### Integrator options
Most of the integrator options are found in `prob%opts`. Here is a list of all other integrator options in `prob%opts`:

 * `constrained`: Set this to `1` if the system is constrained, for unconstrained systems set this to `0`.
 * `stab2`: Only applies in the constrained case. This has to be set to `1`. Although `RATTLie` does not use a direct discretization of the equivalent stabilized index-2 equations of motion, it always enforces the hidden constraints on velocity level.
 * `const_mass_matrix`: Set this to `1` if the mass matrix π does not depend on the configuration π. Set this to `0` if π depends on π. Note that `RATTLie` was not tested with nonconstant mass matrices, results may not be accurate.
 * `diag_mass_matrix`: Set this to `1` if the mass matrix π is a diagonal matrix. In this case `RATTLie_diag_M` will be used instead of `RATTLie_M`. Set this to `0` for full mass matrices.
 * `banded_iteration_matrix`: Set this to `1` if the iteration matrix has band structure or can be rearranged to a matrix with band structure. This is usually the case when integrating a system with several bodies that are chained but otherwise don't interact directly. Set this to `0` if the iteration matrix can be full.
 * `nr_subdiag`: Only applies to the case with banded iteration matrix. Set this to the number of subdiagonals (diagonals below the main diagonal) of the band structure not counting the main diagonal. Note that the Newton method may fail or still converge if `nr_subdiag` is smaller than the actual numer of subdiagonals. If chosen too large, the integration will be slower than it has to be.
 * `nr_superdiag`: See `nr_subdiag`, but for the superdiagonals (diagonals above the main diagonal) instead of the subdiagonals.
 * `jour`: Only applies to the case with banded iteration matrix. This should be an integer vector such that `St(jour,jour)` is a banded matrix, where `St` is the iteration matrix. The size should be:
   * Unconstrained case: `prob%sizev`
   * Constrained case: `prob%sizev+prob%sizel`
 * `use_num_Ct`: Obsolete and is never used.
 * `use_num_Kt`: Obsolete and is never used.
 * `no_Ct`: Obsolete and is never used.
 * `no_Kt`: Obsolete and is never used.
 * `atol`: Set this `real(8)` variable to the absolute tolerance to be used in the Newton method.
 * `rtol`: Set this `real(8)` variable to the relative tolerance to be used in the Newton method.
 * `imax`: Set this variable to the maximum number of iterations after which the Newton iteration is considered to not converge. If this variable is to low, integration might not suceed.
 * `t0`: Set this `real(8)` variable to the beginning of the time integration interval π‘β. Don't forget to set `this%t = this%t0` in `RATTLie_init`.
 * `te`: Set this `real(8)` variable to the end of the time integration interval π‘β. Note, however that in the constrained case, `RATTLie` will perform one more integration step to a value π‘β+β in order to calculate a valid approximation for the Lagrange multiplier at π‘β.
 * `nsteps`: Set this variable to the number of integration steps (of equal length) to be made between π‘β and π‘β. This means the step size can be calculated by β=(π‘β-π‘β)/`nsteps`.

### Compiler flags
This project has been written for `gfortran` on Linux, although porting to code to different compilers or different platforms should not be too hard.
The makefile of `RATTLie` will look for the variable `EXTRAFFLAGS` that is supposed to contain compiler flags. In order to make the variable `EXTRAFFLAGS` visible to the makefile, you should put `export EXTRAFFLAGS` in the problem makefile after defining it.
Here is a list of helpful compiler flags:

 * `-O`: Pass this flag to "optimize" the code. It will produce code that may run (a lot) faster. The optimization can be controlled in levels: `-O0` to turn it off completely, `-O1` and `-O3` exist, see `gfortran`s manual, and `-O2` is equivalent to `-O`.
 * `-Wall`: Turn on "all" warnings while compiling. Helpful for debugging.
 * `-Wextra`: Turn on even more warnings than "all". Helpful for debugging.
 * `-Dpure=''`: This defines the preprocessor variable `pure` to be the empty string. This will cause the Fortran keyword `pure` to be removed from the whole code, making all otherwise `pure` procedures non-`pure`. This is very helpful for debugging, because in a `pure` procedure, no side effects such as printing are allowed. (Exept right before `error stop`). _Note that this is black magic and probably everybody will tell you not to do such things in a program that should produce sensible output._
 * `-g`: Turn on debug mode. Extremely helpful for debugging, obviously. Slows down the program on the other hand.
 * `-fbounds-check`: Check bounds of vectors and such. May help to discover errors.
 * `-fimplicit-none`: Automatically uses `implicit none` everywhere. Use only for debugging and be sure to write `implicit none` everywhere it belongs.
 * `-fbacktrace`: Backtrace. Probably useful for debugging.
 * `-fcheck=all`: Check stuff. Probably useful for debugging.
 * `-finit-real=NaN`: Initialize every `real` variable with `NaN` (not a number). Useful for finding variables that are used before they were defined. (Usually undefined variables have a random value that happens to be in the memory before.)
 * `-finit-integer=-77`: See above. If this flag is given you find an integer with value `-77` you probably forgot to define it.
 * `-ffpe-trap=zero,invalid`: Find the most common floating point errors (fpe) such as: dividing by zero and other invalid stuff. (Usually would result in "undefined behaviour", whatever that means.)

It usually makes sense to define `EXTRAFFLAGS = -O` for producing test results and defining `EXTRAFFLAGS` to be all other mentioned flags for debugging.

### Preprocessor variables
There are a few preprocessor variables that can be defined when compiling `RATTLie`. Some of which are purely for debugging purposes, but some might be really useful. A preprocessor variable `MYVAR` can be defined by giving the compiler (`gfortran`) the compiler flag `-DMYVAR`. Sometimes it is meaningful to give the preprocessor variable a value, say `0`. This can be done by passing the compiler glag `-DMYVAR=0`.
Here is a list of some preprocessor variables:

 * `USE_INDEX_1_APPROX`: If this preprocessor variable is defined, then the Lagrange multiplier `this%l` will be calculated using an index-1 formulation of the equations of motion. The curvature term which appears will be approximated by finite differences.
 * `USE_INDEX_1`: If this preprocessor variable is defined, the Lagrange multiplier `this%l` will be calculated using an index-1 formulation of the equations of motion, but the curvature term is evaluated exactly by using `RATTLie_Z`.
 * `ONLY_LM_LP`: If this preprocessor variable is defined, the Lagrange multiplier `this%l` will be set to the `lp` of the next integration step.
 * `VARIABLE_STEPS`: If this preprocessor variable is defined, instead of equidistant time steps, a predefined set of time steps is used. In this case, the `RATTLie_options` type has an additional variable `tspan` that needs to be allocated and filled before calling the integrator. `tspan` is supposed to contain all time instances that RATTLie will step to. The first and last element should be `t0` and `te`, respectively. The elements of `tspan` don't need to be equally spaced.
 * `TWO_ITERATIONS_MINIMUM`: If this preprocessor variable is defined, in the (first) Newton iteration, at least two iteration steps are performed. This can help to increase the robustness of the integrator, but will make integration slower in places, where one iteration step would have been enough.
 * `RELAXATION`: If this preprocessor variable is defined, after the 10th Newton iteration, a relaxation factor 0.8 will be introduced. This can help to increase the robustness of the integrator, but will might make the integration slower.
 * `DEBUG_PRINT_ITERATION_MATRIX_AT`: This variable should have a valid numeric value if it is defined. In that case, integration will be stopped at the first integration step that is equal or larger than the value and the iteration matrix is printed to standard output.
 * `INC_LAMBDA`: Different way of calculating the `lp` and `lm` Lagrange multipliers. Possibly broken, don't use.
 * `STNUM`: If this preprocessor variable is defined, the iteration matrix will be determined by finite differences completely. This is usually a lot slower and less reliable. This is probably broken, don't use.


## Usage
In order to implement a problem that should be integrated with `RATTLie` the following files are probably needed to be created:

 * `problem.F90` which contains a module `problem` which used the module `RATTLie` and extends the abstract type `RATTLie_problem` to a non-abstract type.
 * `main.F90` which contains the `program main` and uses the module `problem` implemented in `problem.F90`.
 * `makefile` which contains recipies to compile and maybe run the program.

For template files for these three, refer to the documentation of the generalized-Ξ± integrator `gena`. Of course "`gena`" will have to be replaced by "`RATTLie`".

## Related projects
Integrators:

 * [The Lie group generalized-Ξ± method `gena`](https://github.com/StHante/gena)
 * [The Lie group BDF method `BLieDF`](https://github.com/StHante/BLieDF)
 * [The Lie group RATTLE method `RATTLie`](https://github.com/StHante/RATTLie)
 * [The Lie group SHAKE method `SHAKELie`](https://github.com/StHante/SHAKELie)
 * [The nonholonomic RATTLie method `RATTLie_nonhol`](https://github.com/StHante/RATTLie_nonhol)

Test problems:

 * [The heavy top example `heavy_top`](https://github.com/StHante/heavy_top)
 * [The constrained Cosserat beam model `crmS3R3`](https://github.com/StHante/crmS3R3)
 * [The rolling disk example `rolling_disk`](https://github.com/StHante/rolling_disk)

Miscellaneous:

 * [Implementation of Lie group functions `liegroup`](https://github.com/StHante/liegroup)
 * [Expand a config file with different configurations to several files `expandconfig`](https://github.com/StHante/expandconfig)
 * [Read lua files in Matlab and Octave `readLua`](https://github.com/StHante/readLua-for-Matlab-and-Octave)

Third party projects:

 * [Reading lua files in Fortran `aotus`](https://geb.sts.nt.uni-siegen.de/doxy/aotus/)
 * [GFortran](https://gcc.gnu.org/fortran/)
 * [GNU Parallel](https://www.gnu.org/software/parallel/)
