
module RATTLie
   implicit none

   ! definition of integrator options
   type  :: RATTLie_options
      ! system is constrained?
      integer :: constrained = 0
      ! use stabilized index two formulation
      integer :: stab2 = 0
      ! mass matrix is constant?
      integer :: const_mass_matrix = 0
      ! mass matrix is a diagonal matrix? If == 1 then RATTLie_diag_M
      ! is used in order to calculate the mass matrix rather than
      ! RATTLie_M
      integer :: diag_mass_matrix = 0
      ! iteration matrix is banded? (in contrained case: an ordered iteration matrix)
      integer :: banded_iteration_matrix = 0
      ! if iteration matrix is banded: number of subdiagonals
      integer :: nr_subdiag = 0
      ! if iteration matrix is banded: number of superdiagonals
      integer :: nr_superdiag = 0
      ! vector, that is used to order the iteration matrix in order to obtain a banded structure in the constrained case
      integer, dimension(:), allocatable  :: jour
      ! recalculate the iteration matrix in every Newton step
      integer :: recalc_iteration_matrix = 0
      ! use numerical approximation for Ct and Kt
      integer :: use_num_Ct = 1
      integer :: use_num_Kt = 1
      ! omit Ct and Kt resp. in the iteration matrix
      integer :: no_Ct = 0
      integer :: no_Kt = 0
      ! variables for error control of newton method (absolute and relativ tolerance)
      real(8)  :: atol = 1.0e-10_8
      real(8)  :: rtol = 1.0e-8_8
      integer :: imax   = 5
      ! integration span
      real(8) :: t0 = 0.0_8
      real(8) :: te = 1.0_8
      integer :: nsteps = 100
#ifdef VARIABLE_STEPS
      real(8), allocatable, dimension(:)   :: tspan
#endif
   end type RATTLie_options

   ! definition of integrator statistics
   type  :: RATTLie_statistics
      ! current number of newton steps TODO: private
      integer  :: newt_steps_curr = 0
      ! number of newton steps
      integer  :: newt_steps_sum = 0
      ! maximum number of newton steps
      integer  :: newt_steps_max = 0
      ! average number of newton steps
      real(8)  :: newt_steps_avg = 0.0_8
      ! number of calls
      integer  :: ngcalls = 0
      integer  :: nBcalls = 0
      ! integration time
      real(8)  :: time = 0.0_8
   end type RATTLie_statistics

   ! definition of abstract problem type
   type, abstract :: RATTLie_problem
      ! variables that define the state of the integrator
      real(8)                             :: t = 0.0_8
      integer                             :: sizeq = 0
      real(8), dimension(:), allocatable  :: q
      integer                             :: sizev = 0
      real(8), dimension(:), allocatable  :: v
      real(8), dimension(:), allocatable  :: vd
      real(8), dimension(:), allocatable  :: Dq ! this is just for output
      integer                             :: sizel = 0 ! number of constraints, irrelevant for this%opts%constrained == 0
      real(8), dimension(:), allocatable  :: l, lm, lp ! Lagrange multipliers, only needed in the constrained case
      real(8), dimension(:,:), allocatable  :: Bq
      real(8), dimension(:), allocatable    :: fq
      ! integrator options
      type(RATTLie_options)                :: opts
      ! solver statistics
      type(RATTLie_statistics)             :: RATTLie_stats
      ! constant mass matrix (if opts%const_mass_matrix == 1)
      real(8), dimension(:,:), allocatable   :: RATTLie_const_M
      ! constant diagonal mass matrix (if opts%diag_mass_matrix == 1, also)
      real(8), dimension(:), allocatable     :: RATTLie_const_diag_M
      ! internal variables:
      ! RATTLie uses the impulse instead of velocity
      real(8), dimension(:), allocatable  :: p
   ! definition of deferred procedures
   contains
         ! function $M$ in $M(q) \dot v = \hat v^T M v - f(q,v,t)$
      procedure(RATTLie_M),              deferred :: RATTLie_M
         ! function $M$ in $M(q) \dot v = \hat v^T M V - f(q,v,t)$, diagonal case
      procedure(RATTLie_diag_M),         deferred :: RATTLie_diag_M
         ! function $g$ in $M(q) \dot v = \hat v^T M V - f(q,v,t)$
      procedure(RATTLie_f),              deferred :: RATTLie_f
         ! function that gives the inertial forces
      procedure(RATTLie_inertial),       deferred :: RATTLie_inertial
         ! operation in the lie space $q_n * \exp(h\cdot\widetilde{\Delta q_n})$
      procedure(RATTLie_qlpexphDqtilde), deferred :: RATTLie_qlpexphDqtilde
         ! operation in the lie algebra $inversetilde([\tilde{v},\tilde{w}])$, may be a dummy, if system is unconstrained
      procedure(RATTLie_itlbtvtw),       deferred :: RATTLie_itlbtvtw
         ! tangent damping matrix $C_t$
      procedure(RATTLie_Ct),             deferred :: RATTLie_Ct
         ! tangent stiffness matrix $K_t$
      procedure(RATTLie_Kt),             deferred :: RATTLie_Kt
         ! tangent stiffness matrix $K_t$ in the constrained case
         ! (may depend on the Lagrange multipliers)
      procedure(RATTLie_Kt_lambda),      deferred :: RATTLie_Kt_lambda
         ! tangent operator of the exponential map
      procedure(RATTLie_Tg),             deferred :: RATTLie_Tg
         ! transposed inverse of the tangent operator of the exponential map
      procedure(RATTLie_Tg_inv_T),       deferred :: RATTLie_Tg_inv_T
         ! derivative of transposed inverse of the tangent operator of the exponential map
      procedure(RATTLie_d_Tg_inv_T),     deferred :: RATTLie_d_Tg_inv_T
         ! subroutine for output, is called after every integration step
      procedure(RATTLie_outputFunction), deferred :: RATTLie_outputFunction
         ! subroutine to initialize the problem
      procedure(RATTLie_init),           deferred :: RATTLie_init
         ! function $\Phi(q)$
      procedure(RATTLie_phi),            deferred :: RATTLie_phi
         ! function $B(q)$, where $B(q)v = d/dq \Phi(q) * (DL_q(e) v)$
      procedure(RATTLie_b),              deferred :: RATTLie_b
         ! function $F(q,v) = d/(dq) (B(q)v) * (DL_q(e) v)$
      procedure(RATTLie_Z),              deferred :: RATTLie_Z
         ! function Z, where $Zw = Z(q(Dq)) (v(Dq + BT(q)eta), T w)$
      procedure(RATTLie_matZ),           deferred :: RATTLie_matZ
         ! subroutine to calculate the difference in lambda to correct error terms
      procedure                  :: RATTLie_calculate_lambda_correction => RATTLie_calculate_lambda_correction
         ! subroutine to calculate the new lambda and vd
      procedure                  :: RATTLie_calculate_lambda_vd => RATTLie_calculate_lambda_vd
         ! subroutine to calculate correct initial values for vd and a
      procedure                  :: RATTLie_calcInitial => RATTLie_calcInitial ! TODO: private
         ! subroutine to calculate correct initial values for vd and a
      procedure                  :: RATTLie_calcInitialConstrained => RATTLie_calcInitialConstrained
         ! subroutine to integrate one time step
      procedure                  :: RATTLie_solveTimeStep => RATTLie_solveTimeStep ! TODO: private
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: RATTLie_solveConstrainedTimeStep => RATTLie_solveConstrainedTimeStep
         ! subroutine to integrate one time step in a constrained system
      procedure                  :: RATTLie_solveConstrainedTimeStep_stab2 => RATTLie_solveConstrainedTimeStep_stab2
         ! subroutine for integration
      procedure                  :: RATTLie_integrate => RATTLie_integrate
      procedure                  :: RATTLie_print_stats => RATTLie_print_stats
         ! Clean up the RATTLie_problem object (only internal variables
         ! are reset, RATTLie_alpha_m, RATTLie_alpha_f, RATTLie_beta and
         ! RATTLie_gamma are NOT reset)
      procedure                  :: RATTLie_cleanup => RATTLie_cleanup
   end type RATTLie_problem

   ! abstract interface of the problem
   abstract interface
      pure function RATTLie_M(this, q) result(rslt)
         import                                          :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)            :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev,this%sizev)       :: rslt
      end function RATTLie_M

      pure function RATTLie_diag_M(this, q) result(rslt)
         import                                          :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)            :: this
         real(8), dimension(:), intent(in)               :: q
         ! result
         real(8), dimension(this%sizev)                  :: rslt
      end function RATTLie_diag_M

      pure function RATTLie_f(this, q, v, t) result(rslt)
         import                              :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8), dimension(:), intent(in)   :: v
         real(8),               intent(in)   :: t
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function RATTLie_f

      pure function RATTLie_inertial(this, v) result(rslt)
         import                              :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function RATTLie_inertial

      pure function RATTLie_qlpexphDqtilde(this, q, h, dq) result(rslt)
         import                              :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: q
         real(8),               intent(in)   :: h
         real(8), dimension(:), intent(in)   :: dq
         ! result
         real(8), dimension(this%sizeq)      :: rslt
      end function RATTLie_qlpexphDqtilde

      pure function RATTLie_itlbtvtw(this, v, w) result(rslt)
         import                              :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         real(8), dimension(:), intent(in)   :: w
         ! result
         real(8), dimension(this%sizev)      :: rslt
      end function RATTLie_itlbtvtw

      pure function RATTLie_Ct(this, q, v, t) result(rslt)
         import                                    :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function RATTLie_Ct

      pure function RATTLie_Kt(this, q, v, vd, t) result(rslt)
         import                                    :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function RATTLie_Kt

      pure function RATTLie_Kt_lambda(this, q, v, vd, l, t) result(rslt)
         import                                    :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)         :: this
         real(8), dimension(:), intent(in)         :: q
         real(8), dimension(:), intent(in)         :: v
         real(8), dimension(:), intent(in)         :: vd
         real(8), dimension(:), intent(in)         :: l
         real(8),               intent(in)         :: t
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function RATTLie_Kt_lambda

      pure function RATTLie_Tg(this, h, dq) result(rslt)
         import                                    :: RATTLie_problem
         ! input
         class(RATTLie_problem),intent(in)         :: this
         real(8),               intent(in)         :: h
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function RATTLie_Tg

      pure function RATTLie_Tg_inv_T(this, dq) result(rslt)
         import                                    :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)      :: this
         real(8), dimension(:), intent(in)         :: dq
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function RATTLie_Tg_inv_T

      pure function RATTLie_d_Tg_inv_T(this, v, w) result(rslt)
         import                                    :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)      :: this
         real(8), dimension(:), intent(in)         :: v, w
         ! result
         real(8), dimension(this%sizev,this%sizev) :: rslt
      end function RATTLie_d_Tg_inv_T

      pure function RATTLie_norm(this, v) result(rslt)
         import                  :: RATTLie_problem
         ! input
         class(RATTLie_problem),   intent(in)   :: this
         real(8), dimension(:), intent(in)   :: v
         ! result
         real(8)                             :: rslt
      end function RATTLie_norm

      subroutine RATTLie_outputFunction(this,info)
         import                           :: RATTLie_problem
         ! input
         class(RATTLie_problem), intent(in)  :: this
         integer,             intent(in)  :: info
      end subroutine RATTLie_outputFunction

      subroutine RATTLie_init(this)
         import                              :: RATTLie_problem
         ! input/output
         class(RATTLie_problem), intent(inout)  :: this
      end subroutine RATTLie_init

      pure function RATTLie_phi(this,q) result(rslt)
         import                                       :: RATTLie_problem
         ! input
         class(RATTLie_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel)               :: rslt
      end function RATTLie_phi

      pure function RATTLie_B(this,q) result(rslt)
         import                                       :: RATTLie_problem
         ! input
         class(RATTLie_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         ! result
         real(8), dimension(this%sizel,this%sizev)    :: rslt
      end function RATTLie_B

      pure function RATTLie_Z(this,q,v) result(rslt)
         import                                       :: RATTLie_problem
         ! input
         class(RATTLie_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         ! result
         real(8), dimension(this%sizel)               :: rslt
      end function RATTLie_Z

      pure function RATTLie_matZ(this,q,v,T) result(rslt)
         import                                       :: RATTLie_problem
         ! input
         class(RATTLie_problem),              intent(in) :: this
         real(8), dimension(:),            intent(in) :: q
         real(8), dimension(:),            intent(in) :: v
         real(8), dimension(:,:),          intent(in) :: T
         ! result
         real(8), dimension(this%sizel, this%sizev)   :: rslt
      end function RATTLie_matZ

   end interface

   contains

   ! function for calculating the change in lambda to correct error terms
   function RATTLie_calculate_lambda_correction(this, q, quot) result(rslt)
      implicit none
      class(RATTLie_problem), intent(inout)  :: this
      real(8),                intent(in   )  :: q(this%sizeq)
      real(8),                intent(in   )  :: quot(this%sizel)
      ! result
      real(8)                                :: rslt(this%sizel)
      ! internal
      real(8)                                :: mat(this%sizev+this%sizel, this%sizev+this%sizel)
      real(8)                                :: rhs(this%sizev+this%sizel)
      integer                                :: ipiv(this%sizev+this%sizel), info
      integer                                :: i
      !
      associate(sv => this%sizev, &
                sl => this%sizel)
         ! Fill matrix
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               mat(1:sv,1:sv) = 0.0_8
               forall (i=1:sv)
                  mat(i,i) = this%RATTLie_const_diag_M(i)
               end forall
            else
               mat(1:sv, 1:sv) = this%RATTLie_const_M
            end if
         else
            ERROR STOP "Nonconstant mass matrix not supported yet"
         end if
         mat(sv+1:sv+sl,1:sv) = this%RATTLie_B(q)
         this%RATTLie_stats%nBcalls = this%RATTLie_stats%nBcalls + 1
         mat(1:sv,sv+1:sv+sl) = transpose(mat(sv+1:sv+sl,1:sv))
         mat(sv+1:sv+sl,sv+1:sv+sl) = 0.0_8

         ! Fill rhs
         rhs(1:sv) = 0.0_8
         rhs(sv+1:sv+sl) = quot

         ! Solve linear System
         call dgesv(        &! solve the System A*X=B and save the result in B
                     sv+sl, &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                     1,     &! number of right hand sides (=size(B,2))
                     mat,   &! matrix A
                     sv+sl, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                     ipiv,  &! integer pivot vector; it is not needed
                     rhs,   &! matrix B
                     sv+sl, &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                     info)  ! integer information flag

         if (.not. info == 0) then
            error stop "failed to solve linear system in calculate_lambda_correction"
         end if

         ! Result variable
         rslt = rhs(sv+1:sv+sl)
      end associate
   end function RATTLie_calculate_lambda_correction

   ! function for calculating the new lambda and vd
   subroutine RATTLie_calculate_lambda_vd(this, force, curv)
      implicit none
      ! input/output
      class(RATTLie_problem), intent(inout)  :: this
      ! input
      real(8),                intent(in   )  :: force(this%sizev)
      real(8),                intent(in   )  :: curv(this%sizel)
      ! internal
      real(8)                                :: mat(this%sizev+this%sizel, this%sizev+this%sizel)
      real(8)                                :: rhs(this%sizev+this%sizel)
      integer                                :: ipiv(this%sizev+this%sizel), info
      integer                                :: i
      !
      associate(sv => this%sizev, &
                sl => this%sizel)
         ! Fill matrix
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               mat(1:sv,1:sv) = 0.0_8
               forall (i=1:sv)
                  mat(i,i) = this%RATTLie_const_diag_M(i)
               end forall
            else
               mat(1:sv, 1:sv) = this%RATTLie_const_M
            end if
         else
            ERROR STOP "Nonconstant mass matrix not supported yet"
         end if
         mat(sv+1:sv+sl,1:sv) = this%RATTLie_B(this%q)
         this%RATTLie_stats%nBcalls = this%RATTLie_stats%nBcalls + 1
         mat(1:sv,sv+1:sv+sl) = transpose(mat(sv+1:sv+sl,1:sv))
         mat(sv+1:sv+sl,sv+1:sv+sl) = 0.0_8

         ! Fill rhs
         rhs(1:sv) = force
         rhs(sv+1:sv+sl) = curv

         ! Solve linear System
         call dgesv(        &! solve the System A*X=B and save the result in B
                     sv+sl, &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                     1,     &! number of right hand sides (=size(B,2))
                     mat,   &! matrix A
                     sv+sl, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                     ipiv,  &! integer pivot vector; it is not needed
                     rhs,   &! matrix B
                     sv+sl, &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                     info)  ! integer information flag

         if (.not. info == 0) then
            error stop "failed to solve linear system in calculate_lambda_correction"
         end if

         ! Set output
         this%l  = rhs(sv+1:sv+sl)
         this%vd = rhs(1:sv)
      end associate
   end subroutine RATTLie_calculate_lambda_vd


   ! subroutine for pertubing initial values
   subroutine RATTLie_calcInitial(this)
      implicit none
      class(RATTLie_problem),  intent(inout)  :: this  ! problem object
      !

      associate ( t  => this%t,   &
                  q  => this%q,   &
                  p  => this%p,   &
                  v  => this%v,   &
                  sv => this%sizev)
         ! calculate $p_0$
         if (this%opts%diag_mass_matrix == 1) then
            if (this%opts%const_mass_matrix == 1) then
               p = this%RATTLie_const_diag_M*v
            else
               p = this%RATTLie_diag_M(q)*v
            end if
         else
            if (this%opts%const_mass_matrix == 1) then
               p = matmul(this%RATTLie_const_M, v)
            else
               p = matmul(this%RATTLie_M(q), v)
            end if
         end if
      end associate
   end subroutine RATTLie_calcInitial

   subroutine RATTLie_calcInitialConstrained(this, h) !TODO TODO
      implicit none
      class(RATTLie_problem),  intent(inout)  :: this  ! problem object
      real(8),                intent(in   )  :: h     ! step size
      ! internal variables
      integer                                                              :: i
      real(8), dimension(this%sizev+this%sizel)                            :: vdl
      integer, dimension(this%sizev+this%sizel)                            :: ipiv0   ! pivot vector for dgesv for LU factorization of MBB0
      integer                                                              :: info    ! info flag for dgesv
      real(8), dimension(this%sizev+this%sizel, this%sizev+this%sizel)     :: MBB0    ! Matrix
      real(8), dimension(this%sizel, this%sizev)                           :: B0      ! $B(q_0)$
      !

      ERROR STOP 'RATTLie_calcInitialConstrained is not supported.'
   end subroutine RATTLie_calcInitialConstrained

   ! subroutine for integrating one time step
   subroutine RATTLie_solveConstrainedTimeStep(this,t1)
      implicit none
      class(RATTLie_problem), intent(inout)  :: this  ! problem object
      real(8)            , intent(in   )  :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                      :: i     ! for iteration
      integer, dimension(this%sizev+this%sizel)    :: ipiv  ! needed for dgesv
      integer                                      :: info  ! needed for dgesv

      ! intenal real variables
      real(8)                                                           :: h       ! step size
      real(8), dimension(this%sizev)                                    :: vd1     ! $\dot v_{n+1}$
      !real(8), dimension(this%sizev)                                    :: a1      ! $a_{n+1}$ not needed 2015-06-23
      real(8), dimension(this%sizev)                                    :: v1      ! $v_{n+1}$
      real(8), dimension(this%sizev)                                    :: Dq      ! $\Delta q_n$
      real(8), dimension(this%sizeq)                                    :: q1      ! $q_{n+1}}
      real(8), dimension(this%sizev+this%sizel)                         :: res     ! $res(\dots)$
      real(8), dimension(this%sizev+this%sizel,this%sizev+this%sizel)   :: St      ! $S_t$
      real(8), dimension(:,:), allocatable                              :: Stb     ! permuted $S_t$ in banded format
      !real(8)                                                           :: beta_d  ! $\beta'$ not needed 2015-06-23
      !real(8)                                                           :: gamma_d ! $\gamma'$ not needed 2015-06-23
      real(8), dimension(this%sizev+this%sizel)                         :: Dxl     ! $(\Delta x, \Delta \lambda)^T$

      ! internal logical
      logical                                                           :: converged
      ERROR STOP "Not implemented. Use stab2=1"
   end subroutine RATTLie_solveConstrainedTimeStep

   ! subroutine for integrating one time step
   subroutine RATTLie_solveConstrainedTimeStep_stab2(this,t1)
      implicit none
      class(RATTLie_problem), intent(inout)  :: this  ! problem object
      real(8)            , intent(in   )  :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                    :: i     ! for iteration
      integer, dimension(this%sizev+this%sizel)  :: ipiv  ! needed for dgesv
      integer, dimension(this%sizev)  :: ipiv_sv  ! needed for dgesv
      integer                                    :: info  ! needed for dgesv

      ! intenal real variables
      real(8)                                                          :: h       ! step size
      real(8), dimension(this%sizev)                                   :: Dq      ! $v_{k-1/2}$
      real(8), dimension(this%sizeq)                                   :: q1      ! $q_{k}}
      real(8), dimension(this%sizev+this%sizel)                        :: res     ! $res(\dots)$
      !real(8), dimension(this%sizev)                                   :: f
      !real(8), dimension(this%sizel,this%sizev)                        :: Bq
      real(8), dimension(this%sizev)                                   :: M_Dq
      real(8), dimension(this%sizev+this%sizel,this%sizev+this%sizel)  :: St      ! $S_t$
      real(8), dimension(:,:), allocatable                             :: Stb     ! permuted $S_t$ in banded format
      real(8), dimension(this%sizev+this%sizel)                        :: Dxl  ! $(\Delta x, \Delta \lambda, \Delta \eta)^T$
#ifdef INC_LAMBDA
      real(8), dimension(this%sizel)                                   :: l_inc
#endif
      real(8), dimension(this%sizev)                                   :: tmp_sv
      real(8), dimension(this%sizev,this%sizel)                        :: tmp_sv_sl
      real(8)                                                          :: crit
#ifdef RELAXATION
      real(8)                                                          :: relaxation
#endif

      ! internal logical
      logical                                                              :: converged


      ! associate construct for better readability
      associate (sv  => this%sizev,           &
                 sl  => this%sizel,           &
                 v   => this%v,               &
                 p   => this%p,               &
                 q   => this%q,               &
                 t   => this%t,               &
                 lp  => this%lp,              &
                 lm  => this%lm,              &
                 Bq  => this%Bq)
         ! calculation of step size $h$
         h = t1 - t

         ! calculate B(q_n)
         Bq = this%RATTLie_B(q)
         ! count calls
         this%RATTLie_stats%nBcalls = this%RATTLie_stats%nBcalls + 1

         ! initialization of the new values $v_{k-1/2}$ and $\lambda_{k-1}^+$
         Dq  = v
#ifdef INC_LAMBDA
         l_inc = 0.0_8
#else
         lp = lm
#endif



         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax

            ! calculate the new value $q_{k}$
            q1  = this%RATTLie_qlpexphDqtilde(q,h,Dq)

            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq = this%RATTLie_const_diag_M*Dq
               else
                  M_Dq = matmul(this%RATTLie_const_M,Dq)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq = this%RATTLie_diag_M(q)*Dq
               else
                  M_Dq = matmul(this%RATTLie_M(q),Dq)
               end if
            end if

            ! calculate the residue $res_h$
#ifdef INC_LAMBDA
            res(1:sv) = -p  + matmul(this%RATTLie_Tg_inv_T(-h*Dq),M_Dq) + h/2*this%RATTLie_f(q,v,t) + matmul(transpose(Bq),h*lm + l_inc)/2
#else
            res(1:sv) = -p  + matmul(this%RATTLie_Tg_inv_T(-h*Dq),M_Dq) + h/2*this%RATTLie_f(q,v,t) + matmul(transpose(Bq),h*lp)/2
#endif
            res(sv+1:sv+sl) = this%RATTLie_phi(q1) / h
            ! count calls
            this%RATTLie_stats%ngcalls = this%RATTLie_stats%ngcalls + 1

            ! solve the linear System $S_t \cdot \Delta xl = -res$
            Dxl = -res ! the result will be saved in the right hand side's spot


            if (mod(i,10) == 1 .or. this%opts%recalc_iteration_matrix == 1) then
#ifdef STNUM
               St = this%St_num(Dq,h,l,eta,t1)
#else
               ! calculate iteration matrix
               St(1:sv,1:sv) = this%RATTLie_Tg_inv_T(-h*Dq)
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     forall (i=1:sv)
                        St(:,i) = St(:,i)*this%RATTLie_const_diag_M(i)
                     end forall
                  else
                     St(1:sv,1:sv) = matmul(St(1:sv,1:sv), this%RATTLie_const_M)
                  end if
               else
                  ERROR STOP "Nonconstant mass matrix not supported yet"
               end if
               St(1:sv,1:sv) = St(1:sv,1:sv) - h*this%RATTLie_d_Tg_inv_T(-h*Dq, M_Dq)

               St(      1:sv      ,    sv+1:sv+sl   ) = transpose(Bq)/2
               St(sv   +1:sv+sl   ,       1:sv      ) = matmul(this%RATTLie_B(q1), this%RATTLie_Tg(h,Dq))
               St(sv   +1:sv+sl   ,    sv+1:sv+sl   ) = 0.0_8

               ! count calls
               this%RATTLie_stats%nBcalls = this%RATTLie_stats%nBcalls + 1

! endif STNUM
#endif

               ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag)
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1, sv+sl))
                     end if
                     !Stb = 0.0_8 ! DEBUG
                     forall (i=1:sv+sl)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(sv+sl,i+lo)-i+lo+1+hi),i) &
                            = St(this%opts%jour(max(1, i-hi):min(sv+sl, i + lo)), this%opts%jour(i))
                     end forall
                     Dxl = Dxl(this%opts%jour)
                     call dgbsv(           &! solve the system A*X=B and save the result in B
                                 sv+sl,    &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,       &!
                                 hi,       &!
                                 1,        &! number of right hand sides (=size(B,2))
                                 Stb,      &! matrix A
                                 2*lo+hi+1,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,     &! integer pivot vector; it is not needed
                                 Dxl,      &! matrix B
                                 sv+sl,    &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)      ! integer information flag
                     Dxl(this%opts%jour) = Dxl
                     !call print_vector(Dxl,'Dxl_b') ! DEBUG
                     !stop ! DEBUG
                  end associate
               else
                  call dgesv(          &! solve the System A*X=B and save the result in B
                              sv+sl,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,       &! number of right hand sides (=size(B,2))
                              St,      &! matrix A
                              sv+sl,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,    &! integer pivot vector; it is not needed
                              Dxl,     &! matrix B
                              sv+sl,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)     ! integer information flag
                  ! TODO: Errors abfragen
               end if

               if (info.ne.0) print *, "dgesv sagt info=", info

            else ! (St was not (re)calculated)
               ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     Dxl = Dxl(this%opts%jour) !TODO
                     call dgbtrs(               &
                                 'No transpose',&! solve the System A*X=B and save the result in B
                                 sv+sl,         &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,            &!
                                 hi,            &!
                                 1,             &! number of right hand sides (=size(B,2))
                                 Stb,           &! matrix A
                                 2*lo+hi+1,     &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv,          &! integer pivot vector; it is not needed
                                 Dxl,           &! matrix B
                                 sv+sl,         &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)           ! integer information flag
                     Dxl(this%opts%jour) = Dxl
                  end associate
               else
                  call dgetrs(                &
                              'No transpose', & ! solve the System A*X=B and save the result in B
                              sv + sl,        & ! number of linear equations (=size(A,1))  ! Vorsicht: double precision muss real(8) sein
                              1,              & ! number of right hand sides (=size(B,2))
                              St,             & ! matrix A
                              sv + sl,        & ! leading dimesnion of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv,           & ! integer pivot vector; holds the pivoting that was calculated in the dgesv call
                              Dxl,            & ! matrix B
                              sv + sl,        & ! leading dimension of B, in thie case is equal to the number of linear equations (=size(A,1))
                              info)             ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print*, 'dgetrs sagt info=', info
               ! GUBED
            end if

            ! Check for convergence
            crit = ( sum((Dxl(      1:      sv) / (this%opts%atol + this%opts%rtol*abs(Dq  )))**2 )   &
                  +sum((Dxl(sv+   1:   sv+sl) / (this%opts%atol + this%opts%rtol*abs(h*lp)))**2 ) ) &
                   / (sv+sl)
            if (crit <= 1.0_8) then
               converged = .true.
            end if

            ! update
#ifdef RELAXATION
            if (i > 10) then
               relaxation = 0.8_8
            else
               relaxation = 1.0_8
            end if
            Dq = Dq + relaxation*Dxl(1:sv)
            lp = lp + relaxation*Dxl(sv+1:sv+sl)/h
#else
            Dq   = Dq  + Dxl(   1:   sv)
#ifdef INC_LAMBDA
            l_inc = l_inc + Dxl(sv+1:sv+sl)
#else
            lp   = lp  + Dxl(sv+1:sv+sl)/h
#endif
#endif

            ! update solver stats
            this%RATTLie_stats%newt_steps_curr = i

#ifdef TWO_ITERATIONS_MINIMUM
            ! exit loop, if converged
            if ( converged .and. i > 1 ) exit
#else
            ! exit loop, if converged
            if ( converged ) exit
#endif

         end do ! end newton method

         ! DEBUG
         if (this%RATTLie_stats%newt_steps_curr >= this%opts%imax) then
            print *, "Iteration matrix:"
            !call print_matrix(St,'St')
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !! DEBUG !
            call this%RATTLie_print_stats
            print *, 'current iteration:', i
            print *, 'convergence criterium:', crit
            print *, 'norm of last residual:', norm2(res)
            ! DEBUG
            !print *, 'continuing anyways'
            error stop
            !! GUBED !
         end if

         ! The calculation of the new q is done, we now need to calulate the new p
         t  = t1
         q  = this%RATTLie_qlpexphDqtilde(q, h, Dq)
#ifdef INC_LAMBDA
         lp = lm + l_inc/h
#endif

         !! We essentially use one step of a Newton-Method in order to solve
         !! the occuring linear system

         ! initial values
         !p = p
         !lm = 0.0_8 ! TODO, DEBUG: hier war ursprünglich kein Kommentarzeichen

         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               M_Dq = this%RATTLie_const_diag_M*Dq
            else
               M_Dq = matmul(this%RATTLie_const_M,Dq)
            end if
         else
            if (this%opts%diag_mass_matrix == 1) then
               M_Dq = this%RATTLie_diag_M(q)*Dq
            else
               M_Dq = matmul(this%RATTLie_M(q),Dq)
            end if
         end if

         ! Calculate $B(q)$
         Bq = this%RATTLie_B(q)
         this%RATTLie_stats%nBcalls = this%RATTLie_stats%nBcalls + 1

         ! FIXME: Works only, if f does not actually depend on v
         ! calculate the residue $res_h$
         this%fq = this%RATTLie_f(q,v,t)
         this%RATTLie_stats%ngcalls = this%RATTLie_stats%ngcalls + 1
         res(1:sv) = -p  + matmul(this%RATTLie_Tg_inv_T(h*Dq),M_Dq) - h/2*this%fq - matmul(transpose(Bq),h*lm)/2
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               res(sv+1:sv+sl) = matmul(Bq,p/this%RATTLie_const_diag_M)
            else
               tmp_sv = p
               ! Solve linear System
               call dgesv(     &! solve the System A*X=B and save the result in B
                           sv, &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                           1,     &! number of right hand sides (=size(B,2))
                           this%RATTLie_const_M,   &! matrix A
                           sv,    &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                           ipiv_sv,  &! integer pivot vector; it is not needed
                           tmp_sv,   &! matrix B
                           sv, &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                           info)  ! integer information flag
               !
               if (.not. info == 0) then
                  error stop "failed to solve linear system"
               end if
               !
               res(sv+1:sv+sl) = matmul(Bq,tmp_sv)
            end if
         else
            ERROR STOP "Nonconstant mass matrix not supported yet"
         end if

         ! solve the linear System $S_t \cdot \Delta xl = -res$
         Dxl = -res ! the result will be saved in the right hand side's spot

         ! calculate iteration matrix
         St(1:sv,1:sv) = 0.0_8
         forall (i=1:sv) St(i,i) = -1.0_8
         St(1:sv,sv+1:sv+sl) = -transpose(Bq)/2
         St(sv+1:sv+sl,1:sv) = Bq
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               forall (i=1:sv)
                  St(sv+1:sv+sl,i) = St(sv+1:sv+sl,i)/this%RATTLie_const_diag_M(i)
               end forall
            else
               tmp_sv_sl = transpose(St(sv+1:sv+sl,1:sv))
               ! Solve linear System
               call dgesv(     &! solve the System A*X=B and save the result in B
                           sv, &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                           sl,     &! number of right hand sides (=size(B,2))
                           this%RATTLie_const_M,   &! matrix A
                           sv,    &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                           ipiv_sv,  &! integer pivot vector; it is not needed
                           tmp_sv_sl,   &! matrix B
                           sv, &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                           info)  ! integer information flag
               !
               if (.not. info == 0) then
                  error stop "failed to solve linear system"
               end if
               !
               St(sv+1:sv+sl,1:sv) = transpose(tmp_sv_sl)
            end if
         else
            ERROR STOP "Nonconstant mass matrix not supported yet"
         end if
         St(sv+1:sv+sl,sv+1:sv+sl) = 0.0_8

         ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
         if (this%opts%banded_iteration_matrix == 1) then
            associate (lo   => this%opts%nr_subdiag,   &
                       hi   => this%opts%nr_superdiag)
               if (.not. allocated(Stb)) then
                  allocate(Stb(lo*2+hi+1, sv+sl))
               end if
               !Stb = 0.0_8 ! DEBUG
               forall (i=1:sv+sl)
                  Stb((max(1, i-hi)-i+lo+1+hi):(min(sv+sl,i+lo)-i+lo+1+hi),i) &
                      = St(this%opts%jour(max(1, i-hi):min(sv+sl, i + lo)), this%opts%jour(i))
               end forall
               Dxl = Dxl(this%opts%jour)
               call dgbsv(           &! solve the system A*X=B and save the result in B
                           sv+sl,    &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                           lo,       &!
                           hi,       &!
                           1,        &! number of right hand sides (=size(B,2))
                           Stb,      &! matrix A
                           2*lo+hi+1,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                           ipiv,     &! integer pivot vector; it is not needed
                           Dxl,      &! matrix B
                           sv+sl,    &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                           info)      ! integer information flag
               Dxl(this%opts%jour) = Dxl
               !call print_vector(Dxl,'Dxl_b') ! DEBUG
               !stop ! DEBUG
            end associate
         else
            call dgesv(          &! solve the System A*X=B and save the result in B
                        sv+sl,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                        1,       &! number of right hand sides (=size(B,2))
                        St,      &! matrix A
                        sv+sl,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                        ipiv,    &! integer pivot vector; it is not needed
                        Dxl,     &! matrix B
                        sv+sl,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                        info)     ! integer information flag
            ! TODO: Errors abfragen
         end if

         if (info.ne.0) print *, "dgesv sagt info=", info

         ! update
         p   = p  + Dxl(   1:   sv)
         lm  = lm + Dxl(sv+1:sv+sl)/h

         ! Calculate velocity from impulse
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               forall (i=1:sv)
                  v(i) = p(i)/this%RATTLie_const_diag_M(i)
               end forall
            else
               tmp_sv = p
               ! Solve linear System
               call dgesv(     &! solve the System A*X=B and save the result in B
                           sv, &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                           1,     &! number of right hand sides (=size(B,2))
                           this%RATTLie_const_M,   &! matrix A
                           sv,    &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                           ipiv_sv,  &! integer pivot vector; it is not needed
                           tmp_sv,   &! matrix B
                           sv, &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                           info)  ! integer information flag
               !
               if (.not. info == 0) then
                  error stop "failed to solve linear system"
               end if
               !
               v = tmp_sv
            end if
         else
            ERROR STOP "Nonconstant mass matrix not supported yet"
         end if

#ifdef USE_INDEX_1
         !!!!! Calculate vd and lambda from index-1 formulation
         ! Set up matrix
         St = 0.0_8
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               forall (i=1:sv)
                  St(i,i) = this%RATTLie_const_diag_M(i)
               end forall
            else
               St(1:sv,1:sv) = this%RATTLie_const_M
            end if
         else
            ERROR STOP "Nonconstant mass matrix not supported yet"
         end if
         St(1:sv,sv+1:sv+sl) = transpose(Bq)
         St(sv+1:sv+sl,1:sv) = Bq

         ! Set up right hand side
         Dxl(1:sv) = this%RATTLie_inertial(v) - this%fq
         Dxl(sv+1:sv+sl) = -this%RATTLie_Z(q,v)

         ! Actually solve the linear System $S_t \cdot \Delta xl = -res$
         if (this%opts%banded_iteration_matrix == 1) then
            associate (lo   => this%opts%nr_subdiag,   &
                       hi   => this%opts%nr_superdiag)
               if (.not. allocated(Stb)) then
                  allocate(Stb(lo*2+hi+1, sv+sl))
               end if
               !Stb = 0.0_8 ! DEBUG
               forall (i=1:sv+sl)
                  Stb((max(1, i-hi)-i+lo+1+hi):(min(sv+sl,i+lo)-i+lo+1+hi),i) &
                      = St(this%opts%jour(max(1, i-hi):min(sv+sl, i + lo)), this%opts%jour(i))
               end forall
               Dxl = Dxl(this%opts%jour)
               call dgbsv(           &! solve the system A*X=B and save the result in B
                           sv+sl,    &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                           lo,       &!
                           hi,       &!
                           1,        &! number of right hand sides (=size(B,2))
                           Stb,      &! matrix A
                           2*lo+hi+1,&! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                           ipiv,     &! integer pivot vector; it is not needed
                           Dxl,      &! matrix B
                           sv+sl,    &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                           info)      ! integer information flag
               Dxl(this%opts%jour) = Dxl
               !call print_vector(Dxl,'Dxl_b') ! DEBUG
               !stop ! DEBUG
            end associate
         else
            call dgesv(          &! solve the System A*X=B and save the result in B
                        sv+sl,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                        1,       &! number of right hand sides (=size(B,2))
                        St,      &! matrix A
                        sv+sl,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                        ipiv,    &! integer pivot vector; it is not needed
                        Dxl,     &! matrix B
                        sv+sl,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                        info)     ! integer information flag
            ! TODO: Errors abfragen
         end if

         this%l = Dxl(sv+1:sv+sl)
#endif

      end associate
   end subroutine RATTLie_solveConstrainedTimeStep_stab2

   ! subroutine for integrating one time step
   subroutine RATTLie_solveTimeStep(this,t1)
      implicit none
      class(RATTLie_problem), intent(inout)      :: this  ! problem object
      real(8),               intent(in   )      :: t1    ! $t_{n+1}$

      ! internal integer variables
      integer                                   :: i     ! for iteration
      integer, dimension(this%sizev)            :: ipiv  ! needed for dgesv
      integer, dimension(this%sizev)            :: ipiv_sv  ! needed for dgesv
      integer                                   :: info  ! needed for dgesv

      ! internal real variables
      real(8)                                   :: h       ! step size
      real(8), dimension(this%sizev)            :: res     ! $res$
      real(8), dimension(this%sizev,this%sizev) :: St      ! $S_t$
      real(8), dimension(:,:), allocatable      :: Stb     ! $S_t$ in banded format
      real(8), dimension(this%sizev)            :: Dq      ! $v_{1/2}
      real(8), dimension(this%sizev)            :: Dx      ! $\Delta x$
      real(8), dimension(this%sizev)            :: M_Dq
      real(8), dimension(this%sizev)            :: tmp_sv

      ! internal logical
      logical                                   :: converged


      ! associate construct for better readability
      associate (p      => this%p,               &
                 sv     => this%sizev,           &
                 v      => this%v,               &
                 q      => this%q,               &
                 t      => this%t                )
         ! calculation of step size $h$
         h = t1 - t

         ! Initial guess
         ! DEBUG
         !Dq = 0.0_8
         Dq = v

         ! loop of the newton method
         converged = .false.
         do i=1,this%opts%imax
            ! caluclate the residue $res$
            if (this%opts%const_mass_matrix == 1) then
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq = this%RATTLie_const_diag_M*Dq
               else
                  M_Dq = matmul(this%RATTLie_const_M,Dq)
               end if
            else
               if (this%opts%diag_mass_matrix == 1) then
                  M_Dq = this%RATTLie_diag_M(q)*Dq
               else
                  M_Dq = matmul(this%RATTLie_M(q),Dq)
               end if
            end if
            ! Apply inverse tangent operator and add f
            res = -p  + matmul(this%RATTLie_Tg_inv_T(-h*Dq),M_Dq) + h/2*this%RATTLie_f(q,v,t)
            ! count calls
            this%RATTLie_stats%ngcalls = this%RATTLie_stats%ngcalls + 1


            Dx = -res ! the result will be saved in the right hand side's spot

            if (i == 1  .or. this%opts%recalc_iteration_matrix == 1) then
               ! calculate iteration matrix
               St = this%RATTLie_Tg_inv_T(-h*Dq)
               if (this%opts%const_mass_matrix == 1) then
                  if (this%opts%diag_mass_matrix == 1) then
                     forall (i=1:sv)
                        St(:,i) = St(:,i)*this%RATTLie_const_diag_M(i)
                     end forall
                  else
                     St = matmul(St, this%RATTLie_const_M)
                  end if
               else
                  ERROR STOP "Nonconstant mass matrix not supported yet"
               end if
               St = St - h*this%RATTLie_d_Tg_inv_T(-h*Dq, M_Dq)

               !! DEBUG !
               !print *, 'i=', i
               !print *, 't=', t
               !call print_vector(q,'q')
               !call print_vector(v,'v')
               !call print_vector(vd,'vd')
               !call print_vector(a,'a')
               !!call print_matrix(St,'Stn')
               !call print_vector(res, 'res')
               !!call print_vector(this%RATTLie_const_diag_M,'Mn')
               !!print *, "t=", t
               !!if (t > 2.0_8) then
               !!   call print_matrix(St,'St')
               !!   print *, "cond=", mycond(St)
               !!   exit
               !!endif
               !!print *,
               !!stop
               !! GUBED !

               ! solve the linear System $S_t \cdot \Delta x = -res$
               if (this%opts%banded_iteration_matrix == 1) then
#ifdef DEBUG_PRINT_ITERATION_MATRIX_AT
                  if (t >= DEBUG_PRINT_ITERATION_MATRIX_AT) then
                     call print_matrix(St,'St')
                     call print_vector_int(this%opts%jour,'jour')
                     stop "Successfully printed iteration matrix"
                  end if
#endif
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag, &
                             uStb => ubound(Stb,1))
                     if (.not. allocated(Stb)) then
                        allocate(Stb(lo*2+hi+1,this%sizev))
                     end if
                     forall (i=1:this%sizev)
                        Stb((max(1, i-hi)-i+lo+1+hi):(min(this%sizev,i+lo)-i+lo+1+hi),i) &
                            = St(max(1, i-hi):min(this%sizev, i + lo),i)
                     end forall
                     ! TODO: Conditional jump or move depends on unititialised value(s) in the next line?
                     call dgbsv(       &! solve the system A*X=B and save the result in B
                                 sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,   &!
                                 hi,   &!
                                 1,    &! number of right hand sides (=size(B,2))
                                 Stb,  &! matrix A
                                 2*lo+hi+1, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv, &! integer pivot vector; it is not needed
                                 Dx,   &! matrix B
                                 sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)  ! integer information flag
                  end associate
               else
                  call dgesv(       &! solve the System A*X=B and save the result in B
                              sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,    &! number of right hand sides (=size(B,2))
                              St,   &! matrix A
                              sv,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv, &! integer pivot vector; it is not needed
                              Dx,   &! matrix B
                              sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)  ! integer information flag
               end if
               ! TODO: Errors abfragen
               !!! DEBUG
               !   print *, 'St = ...'
               !   print *, '[', St(1,1), ',', St(1,2), ';'
               !   print *, ' ', St(2,1), ',', St(2,2), '];'
               !   print *, 'mres = ...'
               !   print *, '[', -res(1), ',', -res(2), '];'
               !   print *, 'Dx = ...'
               !   print *, '[', Dx(1), ',', Dx(2), '];'
               !   print *, ''
               if (info.ne.0) print *, 'dgesv sagt info=', info
               !! GEBUD
            else
               ! solve the linear System $S_t \cdot \Delta x = -res$
               if (this%opts%banded_iteration_matrix == 1) then
                  associate (lo   => this%opts%nr_subdiag,   &
                             hi   => this%opts%nr_superdiag  )
                     call dgbtrs(      &
                                 'No transpose', &! solve the System A*X=B and save the result in B
                                 sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                                 lo,   &!
                                 hi,   &!
                                 1,    &! number of right hand sides (=size(B,2))
                                 Stb,  &! matrix A
                                 2*lo+hi+1, &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                                 ipiv, &! integer pivot vector; it is not needed
                                 Dx,   &! matrix B
                                 sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                                 info)  ! integer information flag
                  end associate
               else
                  call dgetrs(      &
                              'No transpose', &! solve the System A*X=B and save the result in B
                              sv,   &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                              1,    &! number of right hand sides (=size(B,2))
                              St,   &! matrix A
                              sv,   &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                              ipiv, &!
                              Dx,   &! matrix B
                              sv,   &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                              info)  ! integer information flag
               end if
               ! DEBUG
               if (info.ne.0) print *, 'dgesv sagt info=', info
               ! GEBUD
            end if




            ! Check for convergence
            if ( sum((Dx / (this%opts%atol + this%opts%rtol*abs(Dq)))**2 )  &
                  / (sv)   <=  1.0_8 ) then
               converged = .true.
            end if

! DEBUG!!
#define REL *1
            ! update $\dot v_{n+1}$, $v_{n+1}$ and $\Delta q_n$
            Dq  = Dq  + Dx                       REL
#undef REL
! GUBED!!

            ! update solver stats
            this%RATTLie_stats%newt_steps_curr = i

#ifdef TWO_ITERATIONS_MINIMUM
            ! exit loop, if converged
            if ( converged .and. i > 1 ) exit
#else
            ! exit loop, if converged
            if ( converged ) exit
#endif

         end do ! end newton method

         !! DEBUG
         if (this%RATTLie_stats%newt_steps_curr >= this%opts%imax) then
            print *, ' ' ! DEBUG
            print *, "Newton did not converge at t=", t
            !!! DEBUG !
            !call this%RATTLie_print_stats
            !print *, "Iteration matrix:"
            !call print_matrix(St,'St')
            errorstop
            !!! GUBED !
         end if

         ! Calculate new t, q, p and v
         t = t1
         q = this%RATTLie_qlpexphDqtilde(q,h,Dq)
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               M_Dq = this%RATTLie_const_diag_M*Dq
            else
               M_Dq = matmul(this%RATTLie_const_M,Dq)
            end if
         else
            if (this%opts%diag_mass_matrix == 1) then
               M_Dq = this%RATTLie_diag_M(q)*Dq
            else
               M_Dq = matmul(this%RATTLie_M(q),Dq)
            end if
         end if
         p = matmul(this%RATTLie_Tg_inv_T(h*Dq),M_Dq) - h/2*this%RATTLie_f(q,Dq,t) ! TODO: Hier müsste ein v stehen, kein Dq, ist aber egal, wenn das Modell ungedämpft ist
         ! count calls
         this%RATTLie_stats%ngcalls = this%RATTLie_stats%ngcalls + 1

         ! Calculate velocity from impulse
         if (this%opts%const_mass_matrix == 1) then
            if (this%opts%diag_mass_matrix == 1) then
               forall (i=1:sv)
                  v(i) = p(i)/this%RATTLie_const_diag_M(i)
               end forall
            else
               tmp_sv = p
               ! Solve linear System
               call dgesv(     &! solve the System A*X=B and save the result in B
                           sv, &! number of linear equations (=size(A,1))      ! Vorsicht: double precision muss real(8) sein, sonst gibt es Probleme
                           1,     &! number of right hand sides (=size(B,2))
                           this%RATTLie_const_M,   &! matrix A
                           sv,    &! leading dimension of A, in this case is equal to the number of linear equations (=size(A,1))
                           ipiv_sv,  &! integer pivot vector; it is not needed
                           tmp_sv,   &! matrix B
                           sv, &! leading dimension of B,  in this case is equal to the number of linear equations (=size(B,1)=size(A,1))
                           info)  ! integer information flag
               !
               if (.not. info == 0) then
                  error stop "failed to solve linear system"
               end if
               !
               v = tmp_sv
            end if
         else
            ERROR STOP "Nonconstant mass matrix not supported yet"
         end if

      end associate

   end subroutine RATTLie_solveTimeStep


!   pure function cum_sum_omit_last(vec) result(rslt)
!      ! input
!      real(8), intent(in)  :: vec(:)
!      ! result
!      real(8)              :: rslt(size(vec)-1)
!      ! internal
!      integer              :: i
!      !
!      rslt(1) = vec(1)
!      do i=2,size(vec)-1
!         rslt(i) = rslt(i-1) + vec(i)
!      end do
!   end function cum_sum_omit_last


   ! subroutine to integrate
   subroutine RATTLie_integrate(this)
      implicit none
      ! input/output
      class(RATTLie_problem),  intent(inout)  :: this
      ! internal
      class(RATTLie_problem), allocatable     :: oldstep
      integer  :: n  ! needed for iteration
      real(8)  :: h  ! step size $h$
      real(8)  :: t1 ! next time $t_{n+1}$
      real(4)  :: times(2), time ! for dtime
      ! for correcting the lambda
      !real(8), dimension(:), allocatable   :: old_res, current_res, new_res
      real(8), dimension(:,:), allocatable   :: old_Bq
#ifdef VARIABLE_STEPS
      real(8)                                :: hnew, hold
#endif

      ! problem initialization
      call this%RATTLie_init()

      ! TODO: More error checking
      if (this%opts%constrained == 1             .and. &
          this%opts%banded_iteration_matrix == 1 .and. &
          .not. allocated(this%opts%jour)               ) then
         print *, '[ERROR] opts%jour not allocated, aborting'
         errorstop
      end if
      ! TODO: t0 und t vergleichen (müssen gleich sein)

      ! initialize output function
      call this%RATTLie_outputFunction(0)

#ifndef VARIABLE_STEPS
      ! Calculate step size $h$
      h = (this%opts%te - this%opts%t0)/this%opts%nsteps
#endif

      ! Set stats of solver to zero
      this%RATTLie_stats%newt_steps_curr = 0
      this%RATTLie_stats%newt_steps_sum = 0
      this%RATTLie_stats%newt_steps_max = 0
      this%RATTLie_stats%newt_steps_avg = 0
      this%RATTLie_stats%ngcalls = 0
      this%RATTLie_stats%nBcalls = 0

      ! if mass matrix is constant, calculate it
      if (this%opts%const_mass_matrix == 1) then
         if (this%opts%diag_mass_matrix == 1) then
            if (.not. allocated(this%RATTLie_const_diag_M)) then
               allocate(this%RATTLie_const_diag_M(this%sizev))
            end if
            this%RATTLie_const_diag_M = this%RATTLie_diag_M(this%q)
         else
            if (.not. allocated(this%RATTLie_const_M)) then
               allocate(this%RATTLie_const_M(this%sizev,this%sizev))
            end if
            this%RATTLie_const_M = this%RATTLie_M(this%q)
         end if
      end if

      ! allocate Bq and fq
      allocate(old_Bq(this%sizel,this%sizev))
      allocate(this%Bq(this%sizel,this%sizev))
      allocate(this%fq(this%sizev))

      ! start stopwatch
      time = dtime(times)

      if (this%opts%constrained == 0) then

         ! calculate initial values
         call this%RATTLie_calcInitial()

         ! output for the first time
         call this%RATTLie_outputFunction(1)

         ! integration loop
#ifdef VARIABLE_STEPS
         do n=1,size(this%opts%tspan)-1
            ! Get the next time $t_{n+1}$
            t1 = this%opts%tspan(n+1)
#else
         do n=1,this%opts%nsteps
            ! Calculate the next time $t_{n+1}$
            t1 = this%opts%t0 + n*h
#endif
            ! solve time step
            call this%RATTLie_solveTimeStep(t1)
            ! update solver stats
            this%RATTLie_stats%newt_steps_sum = this%RATTLie_stats%newt_steps_sum &
               + this%RATTLie_stats%newt_steps_curr
            this%RATTLie_stats%newt_steps_avg = real(this%RATTLie_stats%newt_steps_sum,8)/n
            this%RATTLie_stats%newt_steps_max = max( &
               this%RATTLie_stats%newt_steps_max,    &
               this%RATTLie_stats%newt_steps_curr  )
            ! output normally
            call this%RATTLie_outputFunction(1)
         end do


      else
         ! calculate correct initial values
         !call this%RATTLie_calcInitialConstrained(h) DEBUG
         call this%RATTLie_calcInitial()
         ! DEBUG FIXME
         this%lm = 0.0_8
         ! GUBED

         ! Remember old Bq
         this%Bq = this%RATTLie_B(this%q)
         this%RATTLie_stats%nBcalls = this%RATTLie_stats%nBcalls + 1

         !! output for the first time is done via oldstep
         !call this%RATTLie_outputFunction(1)

         !! Residual of old step (dummy)
         !allocate(old_res(this%sizel))
         !old_res = 0.0_8
         !! Residual of this step (dummy)
         !allocate(current_res(this%sizel))
         !current_res = 0.0_8
         !! Residual of next step
         !allocate(new_res(this%sizel))
         !new_res = this%RATTLie_Phi(this%q)
         !this%RATTLie_stats%ngcalls = this%RATTLie_stats%ngcalls + 1

         ! Allocate object to remember last step
         allocate(oldstep, source=this)

         ! integration loop, one extra step so correctly calculate the last lambda
#ifdef VARIABLE_STEPS
         do n=1,size(this%opts%tspan)-1
            ! Get the next time $t_{n+1}$
            t1 = this%opts%tspan(n+1)
#else
         do n=1,this%opts%nsteps+1
            ! Calculate the next time $t_{n+1}$
            t1 = this%opts%t0 + n*h
#endif
            ! Remember old values
            old_Bq = oldstep%Bq
            oldstep%t = this%t
            oldstep%q = this%q
            oldstep%v = this%v
            oldstep%l = this%l
            oldstep%lm = this%lm
            oldstep%lp = this%lp
            oldstep%Bq = this%Bq
            oldstep%fq = this%fq

            oldstep%RATTLie_stats%newt_steps_curr = this%RATTLie_stats%newt_steps_curr
            oldstep%RATTLie_stats%newt_steps_avg = this%RATTLie_stats%newt_steps_avg
            oldstep%RATTLie_stats%newt_steps_max = this%RATTLie_stats%newt_steps_max


            if (this%opts%stab2 == 1) then
               ! solve time step with the stabilized index-2 system
               call this%RATTLie_solveConstrainedTimeStep_stab2(t1)
            else
               ERROR STOP "RATTLie does not support index-3 integration"
               ! solve time step
               call this%RATTLie_solveConstrainedTimeStep(t1)
            end if

            !print *, this%opts%tspan(n), oldstep%t, this%opts%tspan(n) - oldstep%t

#ifdef ONLY_LM_LP
   !oldstep%lm = oldstep%lm
   oldstep%l = this%lp
#else
!ONLY_LM_LP is not defined
#ifndef USE_INDEX_1
#ifdef USE_INDEX_1_APPROX
            if (n==1) then
               oldstep%lm = this%lp
               oldstep%l = this%lp
            else
#ifdef VARIABLE_STEPS
               hold = this%opts%tspan(n) - this%opts%tspan(n-1)
               hnew = this%opts%tspan(n+1) - this%opts%tspan(n)
               call oldstep%RATTLie_calculate_lambda_vd(oldstep%RATTLie_inertial(oldstep%v) - oldstep%fq, matmul(hnew/(hold*(hnew+hold))*old_Bq + (1.0_8/hnew - 1.0_8/hold)*oldstep%Bq - hold/(hnew*(hnew+hold))*this%Bq, oldstep%v))
#else
!VARIABLE_STEPS is not defined
               call oldstep%RATTLie_calculate_lambda_vd(oldstep%RATTLie_inertial(oldstep%v) - oldstep%fq, -matmul(this%Bq - old_Bq, oldstep%v)/(2*h))
#endif
!for ifdef VARIABLE_STEPS
               oldstep%lm = (oldstep%lm + this%lp)/2
            end if
#else
!USE_INDEX_1_APPROX is not defined
            !! shift residuals and calculate the new
            !old_res = current_res
            !current_res = new_res
            !old_res = new_res
            !new_res = this%RATTLie_Phi(this%q)
            !this%RATTLie_stats%ngcalls = this%RATTLie_stats%ngcalls + 1

            ! Correct lambdas
            !this%lm = this%lm + this%RATTLie_calculate_lambda_correction(this%q, -2*(new_res - old_res)/h**2)
            !this%lp = this%lp + this%RATTLie_calculate_lambda_correction(oldstep%q, 2*(new_res - old_res)/h**2)

            ! calculate Lagrange multipliers of old step and output
            if (n==1) then
               ! FIXME
               oldstep%l = this%lp
            else
               !oldstep%l  = (oldstep%lm + this%lp)/2 + this%RATTLie_calculate_lambda_correction(oldstep%q, -old_res/h**2 + 2*current_res/h**2 - new_res/h**2)

               ! Calculate interpolated lambda
#ifdef VARIABLE_STEPS
               hold = this%opts%tspan(n) - this%opts%tspan(n-1)
               hnew = this%opts%tspan(n+1) - this%opts%tspan(n)
               oldstep%l = (hnew*(oldstep%lm + oldstep%lp)/2  +  hold*(this%lm + this%lp)/2)/(hnew + hold)
#else
!VARIABLE_STEPS is not defined
               oldstep%l = (oldstep%lm + oldstep%lp + this%lp + this%lm)/4
#endif
!for ifdef VARIABLE_STEPS
               oldstep%lm = (oldstep%lm + this%lp)/2

               !oldstep%l = this%lp
               !oldstep%lm = oldstep%lm
            end if
#endif
!for ifdef USE_INDEX_1_APPROX
#else
!USE_INDEX_1 is defined
            ! Lagrange multipliers oldstel%l were already calculated in  RATTLie_solveConstrainedTimeStep_stab2
            ! This is only to output the Lagrange multipliers calulated by lm and lp
            if (n==1) then
               ! FIXME
               oldstep%lm = this%lp
            else
               oldstep%lm = (oldstep%lm + this%lp)/2
            end if
#endif
!for ifndef USE_INDEX_1
#endif
!for ifdef ONLY_LM_LP

            ! update solver stats
            this%RATTLie_stats%newt_steps_sum = this%RATTLie_stats%newt_steps_sum &
               + this%RATTLie_stats%newt_steps_curr
            this%RATTLie_stats%newt_steps_avg = real(this%RATTLie_stats%newt_steps_sum,8)/n
            this%RATTLie_stats%newt_steps_max = max( &
               this%RATTLie_stats%newt_steps_max,    &
               this%RATTLie_stats%newt_steps_curr  )

            ! output old step
            call oldstep%RATTLie_outputFunction(1)
         end do

         ! Deallocate variable oldstep
         deallocate(oldstep)

         !! Deallocate variables for residuals
         !deallocate(old_res)
         !deallocate(current_res)
         !deallocate(new_res)

      end if

      ! stop stopwatch
      this%RATTLie_stats%time = dtime(times) - time

      ! deallocate Bq and fq
      deallocate(this%Bq)
      deallocate(old_Bq)
      deallocate(this%fq)

      ! output to terminate
      call this%RATTLie_outputFunction(99)

   end subroutine RATTLie_integrate

! DEBUG
   pure function RATTLie_num_B(this, q) result(rslt)
      ! input
      class(RATTLie_problem),   intent(in)         :: this
      real(8), dimension(:), intent(in)         :: q
      ! result
      real(8), dimension(this%sizel,this%sizev) :: rslt
      ! internal
      integer                                   :: i
      real(8), dimension(this%sizev)            :: w
      real(8), dimension(this%sizev)            :: Phi0
      !
      !
      ! Calculate Phi
      Phi0 = this%RATTLie_Phi(q)

      ! set w to zero
      w = 0.0_8
      ! loop over the columns of Ct
      do i=1,this%sizev
         if (.not. i == 1) w(i-1) = 0.0_8
         w(i) = max( abs(q(i))*1.0e-8_8,  1.0e-12_8)
         rslt(:,i) = (this%RATTLie_Phi(this%RATTLie_qlpexphDqtilde(q,1.0_8, w)) - Phi0) / w(i)
         !w    = 0.0_8
         !w(i) = h
         !rslt(:,i) = (this%RATTLie_g(q+this%RATTLie_tilde(w),v,t) - g0)/h
      end do
   end function RATTLie_num_B
! GUBED


   subroutine RATTLie_print_stats(this)
      ! input
      class(RATTLie_problem), intent(in)  :: this
      !
      print *, 'time:          ', this%RATTLie_stats%time
      print *, '#calls of g:   ', this%RATTLie_stats%ngcalls
      print *, '#calls of B:   ', this%RATTLie_stats%nBcalls
      print *, 'newt_steps_max:', this%RATTLie_stats%newt_steps_max
      print *, 'newt_steps_avg:', this%RATTLie_stats%newt_steps_avg
   end subroutine

   subroutine print_matrix(A,Aname)
      implicit none
      ! input
      real(8), dimension(:,:) :: A
      character(len=*)       :: Aname
      ! internal
      integer                 :: i
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      do i=1,ubound(A,1) - 1
         print *, A(i,:), ';'
      end do
      print *, A(ubound(A,1),:), '];'
   end subroutine print_matrix

   subroutine print_vector_int(A,Aname)
      implicit none
      ! input
      integer, dimension(:)   :: A
      character(len=*)        :: Aname
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      print *, A, '];'
   end subroutine print_vector_int

   subroutine print_vector(A,Aname)
      implicit none
      ! input
      real(8), dimension(:)   :: A
      character(len=*)        :: Aname
      !
      print *, '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, trim(Aname), ' = ['
      print *, A, '];'
   end subroutine print_vector

   subroutine RATTLie_cleanup(this)
      implicit none
      ! input/output
      class(RATTLie_problem), intent(inout)  :: this
      !
      this%opts%constrained = 0
      this%opts%const_mass_matrix = 0
      this%opts%diag_mass_matrix = 0
      this%opts%banded_iteration_matrix = 0
      this%opts%nr_subdiag = 0
      this%opts%nr_superdiag = 0
      this%opts%recalc_iteration_matrix = 0
      this%opts%use_num_Ct = 1
      this%opts%use_num_Kt = 1
      this%opts%atol = 1.0e-10_8
      this%opts%rtol = 1.0e-8_8
      this%opts%imax   = 5
      this%opts%t0 = 0.0_8
      this%opts%te = 1.0_8
      this%opts%nsteps = 100
      !
      this%RATTLie_stats%newt_steps_curr = 0
      this%RATTLie_stats%newt_steps_max = 0
      this%RATTLie_stats%newt_steps_avg = 0.0_8
      this%RATTLie_stats%ngcalls = 0
      this%RATTLie_stats%nBcalls = 0
      this%RATTLie_stats%time = 0.0_8
      !
      this%t = 0.0_8
      this%sizeq = 0
      this%sizev = 0
      this%sizel = 0
      !
      if (allocated(this%q))  deallocate(this%q)
      if (allocated(this%v))  deallocate(this%v)
      if (allocated(this%p))  deallocate(this%p)
      if (allocated(this%l))  deallocate(this%l)
      if (allocated(this%lm))  deallocate(this%lm)
      if (allocated(this%lp))  deallocate(this%lp)
      if (allocated(this%RATTLie_const_M)) deallocate(this%RATTLie_const_M)
      if (allocated(this%RATTLie_const_diag_M)) deallocate(this%RATTLie_const_diag_M)
      if (allocated(this%opts%jour)) deallocate(this%opts%jour)
   end subroutine RATTLie_cleanup

#if 1
   function mycond(A) result(rslt)
      implicit none
      ! input
      real(8), intent(in)  :: A(:,:)
      ! result
      real(8)              :: rslt
      ! internal
      real(8)              :: AA(size(A,1),size(A,2))
      real(8)              :: S(min(size(A,1), size(A,2)))
      real(8), allocatable :: work(:)
      integer              :: lwork
      real(8)              :: U(size(A,1),size(A,1))
      real(8)              :: VT(size(A,2),size(A,2))
      integer              :: iwork(8*min(size(A,1), size(A,2)))
      integer              :: info
      !
      AA = A
      lwork = 3*min(size(A,1),size(A,2)) + max(max(size(A,1),size(A,2)),7*min(size(A,1),size(A,2)))
      allocate(work(lwork))
      call DGESDD('A',       &!JOBZ, calculate singular vectors TODO: Will fail if they are not calculated
                  size(A,1), &!M,
                  size(A,2), &!N,
                  AA,        &!A,
                  size(A,1), &!LDA,
                  S,         &!S,
                  U,         &!U,
                  size(A,1), &!LDU,
                  VT,        &!VT,
                  size(A,2), &!LDVT,
                  work,      &!WORK,
                  -1,        &!LWORK,
                  iwork,     &!IWORK
                  info)      !INFO )
      if (info /= 0) then
         errorstop "In cond: dgesdd did not succeed"
      end if
      AA = A
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      call DGESDD('A',       &!JOBZ, calculate singular vectors TODO: Will fail if they are not calculated
                  size(A,1), &!M,
                  size(A,2), &!N,
                  AA,        &!A,
                  size(A,1), &!LDA,
                  S,         &!S,
                  U,         &!U,
                  size(A,1), &!LDU,
                  VT,        &!VT,
                  size(A,2), &!LDVT,
                  work,      &!WORK,
                  lwork,     &!LWORK,
                  iwork,     &!IWORK
                  info)      !INFO )
      if (info /= 0) then
         errorstop "In cond: dgesdd did not succeed"
      end if
      if (abs(S(min(size(A,1),size(A,2)))) < 1.0e-16_8) then
         errorstop "In cond: Matrix A is singular"
      end if
      rslt = S(1)/S(min(size(A,1),size(A,2)))
   end function mycond
#endif

end module RATTLie
