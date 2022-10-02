module ode_lsoda

  implicit none
  private

  public :: init_lsoda, DLSODA

  ! ! DLS001 is declared in subroutines DLSODA, DINTDY, DSTODA, DPRJA, and DSOLSY.
  ! double precision, private :: ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, &
  !     TN, UROUND
  ! integer, private :: INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6), &
  !     ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
  !     LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
  !     MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
  !
  ! ! DLSA01 is declared in subroutines DLSODA, DSTODA, and DPRJA.
  ! double precision, private :: TSW, ROWNS2(20), PDNORM
  ! integer, private :: INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
  ! save

contains

!----------------------------------------------------------------------
subroutine init_lsoda(pb)
! Mostly just to create the array inside the pb class
  use constants, only: FID_SCREEN
  use problem_class, only: problem_type
  use friction, only: get_Jac
  use derivs_all
  use fault_stress, only : compute_stress
  use utils, only : pack

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt !dydt
!  integer :: i !ind_stress_coupling, ind_localisation

  write (FID_SCREEN, *) "Initialising DLSODA solver..."
  !lrn = 20 + 16*pb%neqs*pb%mesh%nn
! Because there are only dv and dtheta.
  pb%lsoda%neq(1) = 2
  pb%lsoda%rtol(1) = pb%acc
!  pb%lsoda%atol(1) = pb%acc
  pb%lsoda%atol = 0d0
!  pb%lsoda%lrw = 22 + 9*pb%neqs * pb%neqs**2
  pb%lsoda%lrw = 22 + pb%lsoda%neq(1)*max(16, pb%neqs+9)
  pb%lsoda%liw = 20 + pb%lsoda%neq(1)

! Lapusta2009's dt
!  pb%lsoda%dt_min = pb%mesh%dx / pb%beta / 3

! Initialize the states to be updated in each run
  allocate(pb%lsoda%rwork(pb%lsoda%lrw))
  allocate(pb%lsoda%iwork(pb%lsoda%liw))
  allocate(pb%lsoda%t(pb%mesh%nn))
  allocate(pb%lsoda%istate(pb%mesh%nn))
  pb%lsoda%istate = 1

  yt(1::pb%neqs) = pb%theta
  yt(2::pb%neqs) = pb%v
  yt(3::pb%neqs) = pb%sigma
  pb%lsoda%t = pb%time

  call pack(yt, pb%theta, pb%v, pb%sigma, pb%theta2, pb%slip, pb)
!   call compute_stress(pb%dtau_dt, dydt(3::pb%neqs), pb%kernel, yt(2)-pb%vpl)
!   write (FID_SCREEN, *) "Stress computed"
!   do i=1, pb%mesh%nn
!     call DLSODA(derivs_lsoda, pb%lsoda%neq, yt((i-1)*pb%neqs+1:i*pb%neqs), &
!                 pb%lsoda%t(i), pb%lsoda%t(i), 1, pb%lsoda%rtol, pb%lsoda%atol, &
!                 2, pb%lsoda%istate, 0, pb%lsoda%rwork, &
!                 pb%lsoda%lrw, pb%lsoda%iwork, pb%lsoda%liw, get_Jac, 2)
!     if (pb%lsoda%istate < 0)  then
!       write(FID_SCREEN, *) "Something wrong with initialization with error:", pb%lsoda%istate
!       stop "LSODA initialization failed"
! !    else
! !      write(FID_SCREEN, *) "DLSODA initialization complete with state:", pb%lsoda%istate
!     endif
!   enddo
  write (FID_SCREEN, *) "DLSODA initialized!"
end subroutine


!----------------------------------------------------------------------
subroutine DLSODA (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
                   ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
      implicit none

      EXTERNAL F, JAC
      INTEGER :: ITOL, ITASK, ISTATE, IOPT, LRW, LIW, JT
      integer :: NEQ(*)
      integer :: IWORK(LIW)
      double precision :: RWORK(LRW)
      double precision :: TOUT, T, Y(*), RTOL(*), ATOL(*)

!-----------------------------------------------------------------------
! This is the 12 November 2003 version of
! DLSODA: Livermore Solver for Ordinary Differential Equations, with
!         Automatic method switching for stiff and nonstiff problems.
!
! This version is in double precision.
!
! DLSODA solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
!
! This a variant version of the DLSODE package.
! It switches automatically between stiff and nonstiff methods.
! This means that the user does not have to determine whether the
! problem is stiff or not, and the solver will automatically choose the
! appropriate method.  It always starts with the nonstiff method.
!
! Authors:       Alan C. Hindmarsh
!                Center for Applied Scientifi! Computing, L-561
!                Lawrence Livermore National Laboratory
!                Livermore, CA 94551
! and
!                Linda R. Petzold
!                Univ. of California at Santa Barbara
!                Dept. of Computer Science
!                Santa Barbara, CA 93106
!
! References:
! 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
!     Solvers, in Scientifi! Computing, R. S. Stepleman et al. (Eds.),
!     North-Holland, Amsterdam, 1983, pp. 55-64.
! 2.  Linda R. Petzold, Automati! Selection of Methods for Solving
!     Stiff and Nonstiff Systems of Ordinary Differential Equations,
!     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
!-----------------------------------------------------------------------
! Summary of Usage.
!
! Communication between the user and the DLSODA package, for normal
! situations, is summarized here.  This summary describes only a subset
! of the full set of options available.  See the full description for
! details, including alternative treatment of the Jacobian matrix,
! optional inputs and outputs, nonstandard options, and
! instructions for special situations.  See also the example
! problem (with program and output) following this summary.
!
! A. First provide a subroutine of the form:
!               SUBROUTINE F (NEQ, T, Y, YDOT)
!               DOUBLE PRECISION T, Y(*), YDOT(*)
! which supplies the vector function f by loading YDOT(i) with f(i).
!
! B. Write a main program which calls Subroutine DLSODA once for
! each point at which answers are desired.  This should also provide
! for possible use of logical unit 6 for output of error messages
! by DLSODA.  On the first call to DLSODA, supply arguments as follows:
! F      = name of subroutine for right-hand side vector f.
!          This name must be declared External in calling program.
! NEQ    = number of first order ODEs.
! Y      = array of initial values, of length NEQ.
! T      = the initial value of the independent variable.
! TOUT   = first point where output is desired (.ne. T).
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
! RTOL   = relative tolerance parameter (scalar).
! ATOL   = absolute tolerance parameter (scalar or array).
!          the estimated local error in y(i) will be controlled so as
!          to be less than
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution: actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of y at t = TOUT.
! ISTATE = integer flag (input and output).  Set ISTATE = 1.
! IOPT   = 0 to indicate no optional inputs used.
! RWORK  = real work array of length at least:
!             22 + NEQ * MAX(16, NEQ + 9).
!          See also Paragraph E below.
! LRW    = declared length of RWORK (in user's dimension).
! IWORK  = integer work array of length at least  20 + NEQ.
! LIW    = declared length of IWORK (in user's dimension).
! JAC     = name of subroutine for Jacobian matrix.
!          Use a dummy name.  See also Paragraph E below.
! JT     = Jacobian type indicator.  Set JT = 2.
!          See also Paragraph E below.
! Note that the main program must declare arrays Y, RWORK, IWORK,
! and possibly ATOL.
!
! C. The output from the first call (or any call) is:
!      Y = array of computed values of y(t) vector.
!      T = corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if DLSODA was successful, negative otherwise.
!          -1 means excess work done on this call (perhaps wrong JT).
!          -2 means excess accuracy requested (tolerances too small).
!          -3 means illegal input detected (see printed message).
!          -4 means repeated error test failures (check all inputs).
!          -5 means repeated convergence failures (perhaps bad Jacobian
!             supplied or wrong choice of JT or tolerances).
!          -6 means error weight became zero during problem. (Solution
!             component i vanished, and ATOL or ATOL(i) = 0.)
!          -7 means work space insufficient to finish (see messages).
!
! D. To continue the integration after a successful return, simply
! reset TOUT and call DLSODA again.  No other parameters need be reset.
!
! E. Note: If and when DLSODA regards the problem as stiff, and
! switches methods accordingly, it must make use of the NEQ by NEQ
! Jacobian matrix, J = df/dy.  For the sake of simplicity, the
! inputs to DLSODA recommended in Paragraph B above cause DLSODA to
! treat J as a full matrix, and to approximate it internally by
! difference quotients.  Alternatively, J can be treated as a band
! matrix (with great potential reduction in the size of the RWORK
! array).  Also, in either the full or banded case, the user can supply
! J in closed form, with a routine whose name is passed as the JAC
! argument.  These alternatives are described in the paragraphs on
! RWORK, JAC, and JT in the full description of the call sequence below.
!
!-----------------------------------------------------------------------
! Full description of user interface to DLSODA.
!
! The user interface to DLSODA consists of the following parts.
!
! 1.   The call sequence to Subroutine DLSODA, which is a driver
!      routine for the solver.  This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      following these descriptions is a description of
!      optional inputs available through the call sequence, and then
!      a description of optional outputs (in the work arrays).
!
! 2.   Descriptions of other routines in the DLSODA package that may be
!      (optionally) called by the user.  These provide the ability to
!      alter error message handling, save and restore the internal
!      Common, and obtain specified derivatives of the solution y(t).
!
! 3.   Descriptions of Common blocks to be declared in overlay
!      or similar environments, or to be saved when doing an interrupt
!      of the problem and continued solution later.
!
! 4.   Description of a subroutine in the DLSODA package,
!      which the user may replace with his/her own version, if desired.
!      this relates to the measurement of errors.
!
!-----------------------------------------------------------------------
! Part 1.  Call Sequence.
!
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, JT,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RWORK and IWORK are also used for conditional and
! optional inputs and optional outputs.  (The term output here refers
! to the return from Subroutine DLSODA to the user's calling program.)
!
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 on input.
!
! The descriptions of the call arguments are as follows.
!
! F      = the name of the user-supplied subroutine defining the
!          ODE system.  The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  Subroutine F is to
!          compute the function f.  It is to have the form
!               SUBROUTINE F (NEQ, T, Y, YDOT)
!               DOUBLE PRECISION T, Y(*), YDOT(*)
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output.  Y and YDOT are arrays of length NEQ.
!          Subroutine F should not alter Y(1),...,Y(NEQ).
!          F must be declared External in the calling program.
!
!          Subroutine F may access user-defined quantities in
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!          (dimensioned in F) and/or Y has length exceeding NEQ(1).
!          See the descriptions of NEQ and Y below.
!
!          If quantities computed in the F routine are needed
!          externally to DLSODA, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use DINTDY instead.
!
! NEQ    = the size of the ODE system (number of first order
!          ordinary differential equations).  Used only for input.
!          NEQ may be decreased, but not increased, during the problem.
!          If NEQ is decreased (with ISTATE = 3 on input), the
!          remaining components of Y should be left undisturbed, if
!          these are to be accessed in F and/or JAC.
!
!          Normally, NEQ is a scalar, and it is generally referred to
!          as a scalar in this user interface description.  However,
!          NEQ may be an array, with NEQ(1) set to the system size.
!          (The DLSODA package accesses only NEQ(1).)  In either case,
!          this parameter is passed as the NEQ argument in all calls
!          to F and JAC.  Hence, if it is an array, locations
!          NEQ(2),... may be used to store other integer data and pass
!          it to F and/or JAC.  Subroutines F and/or JAC must include
!          NEQ in a Dimension statement in that case.
!
! Y      = a real array for the vector of dependent variables, of
!          length NEQ or more.  Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          On the first call, Y must contain the vector of initial
!          values.  On output, Y contains the computed solution vector,
!          evaluated at T.  If desired, the Y array may be used
!          for other purposes between calls to the solver.
!
!          This array is passed as the Y argument in all calls to
!          F and JAC.  Hence its length may exceed NEQ, and locations
!          Y(NEQ+1),... may be used to store other real data and
!          pass it to F and/or JAC.  (The DLSODA package accesses only
!          Y(1),...,Y(NEQ).)
!
! T      = the independent variable.  On input, T is used only on the
!          first call, as the initial point of the integration.
!          on output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          on an error return, T is the farthest point reached.
!
! TOUT   = the next value of t at which a computed solution is desired.
!          Used only for input.
!
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should .ne. T for the next call.
!          For the initial t, an input value of TOUT .ne. T is used
!          in order to determine the direction of the integration
!          (i.e. the algebrai! sign of the step sizes) and the rough
!          scale of the problem.  Integration in either direction
!          (forward or backward in t) is permitted.
!
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT .ne. T).
!          Otherwise, TOUT is required on every call.
!
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal T interval, whose endpoints are
!          TCUR - HU and TCUR (see optional outputs, below, for
!          TCUR and HU).
!
! ITOL   = an indicator for the type of error control.  See
!          description below under ATOL.  Used only for input.
!
! RTOL   = a relative error tolerance parameter, either a scalar or
!          an array of length NEQ.  See description below under ATOL.
!          Input only.
!
! ATOL   = an absolute error tolerance parameter, either a scalar or
!          an array of length NEQ.  Input only.
!
!             The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver.  The solver will
!          control the vector E = (E(i)) of estimated local errors
!          in y, according to an inequality of the form
!                      max-norm of ( E(i)/EWT(i) )   .le.   1,
!          where EWT = (EWT(i)) is a vector of positive error weights.
!          The values of RTOL and ATOL should all be non-negative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!
!             ITOL    RTOL       ATOL          EWT(i)
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting a
!          user-supplied routine for the setting of EWT.
!          See Part 4 below.
!
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL (i.e. of EWT) should be scaled
!          down uniformly.
!
! ITASK  = an index specifying the task to be performed.
!          Input only.  ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration.  This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RWORK(1).
!
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT (exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at t = TOUT are returned first).
!
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!
!          On input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK.  Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,
!             and any optional inputs except H0, MXORDN, and MXORDS.
!             (See IWORK description for ML and MU.)
!          Note:  A preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (Such a call is sometimes useful for the
!          purpose of outputting the initial conditions.)
!          Thus the first call for which TOUT .ne. T requires
!          ISTATE = 1 on input.
!
!          On output, ISTATE has the following values and meanings.
!           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T.  (MXSTEP is an optional input
!              and is normally 500.)  To continue, the user may
!              simply reset ISTATE to a value .gt. 1 and call again
!              (the excess work step counter will be reset to 0).
!              In addition, the user may increase MXSTEP to avoid
!              this error return (see below on optional inputs).
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  This was detected before
!              completing the requested task, but the integration
!              was successful as far as T.  To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3.  The optional output TOLSF may be used for this
!              purpose.  (Note: If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  See written message for details.
!              Note:  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration.  Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!          -7  means the length of RWORK and/or IWORK was too small to
!              proceed, but the integration was successful as far as T.
!              This happens when DLSODA chooses to switch methods
!              but LRW and/or LIW is too small for the new method.
!
!          Note:  Since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other inputs, before
!          calling the solver again.
!
! IOPT   = an integer flag to specify whether or not any optional
!          inputs are being used on this call.  Input only.
!          The optional inputs are listed separately below.
!          IOPT = 0 means no optional inputs are being used.
!                   default values will be used in all cases.
!          IOPT = 1 means one or more optional inputs are being used.
!
! RWORK  = a real array (double precision) for work space, and (in the
!          first 20 words) for conditional and optional inputs and
!          optional outputs.
!          As DLSODA switches automatically between stiff and nonstiff
!          methods, the required length of RWORK can change during the
!          problem.  Thus the RWORK array passed to DLSODA can either
!          have a static (fixed) length large enough for both methods,
!          or have a dynamic (changing) length altered by the calling
!          program in response to output from DLSODA.
!
!                       --- Fixed Length Case ---
!          If the RWORK length is to be fixed, it should be at least
!               MAX (LRN, LRS),
!          where LRN and LRS are the RWORK lengths required when the
!          current method is nonstiff or stiff, respectively.
!
!          The separate RWORK length requirements LRN and LRS are
!          as follows:
!          IF NEQ is constant and the maximum method orders have
!          their default values, then
!             LRN = 20 + 16*NEQ,
!             LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2,
!             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5.
!          Under any other conditions, LRN and LRS are given by:
!             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ,
!             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT,
!          where
!             NYH    = the initial value of NEQ,
!             MXORDN = 12, unless a smaller value is given as an
!                      optional input,
!             MXORDS = 5, unless a smaller value is given as an
!                      optional input,
!             LMAT   = length of matrix work space:
!             LMAT   = NEQ**2 + 2              if JT = 1 or 2,
!             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.
!
!                       --- Dynamic Length Case ---
!          If the length of RWORK is to be dynamic, then it should
!          be at least LRN or LRS, as defined above, depending on the
!          current method.  Initially, it must be at least LRN (since
!          DLSODA starts with the nonstiff method).  On any return
!          from DLSODA, the optional output MCUR indicates the current
!          method.  If MCUR differs from the value it had on the
!          previous return, or if there has only been one call to
!          DLSODA and MCUR is now 2, then DLSODA has switched
!          methods during the last call, and the length of RWORK
!          should be reset (to LRN if MCUR = 1, or to LRS if
!          MCUR = 2).  (An increase in the RWORK length is required
!          if DLSODA returned ISTATE = -7, but not otherwise.)
!          After resetting the length, call DLSODA with ISTATE = 3
!          to signal that change.
!
! LRW    = the length of the array RWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! IWORK  = an integer array for work space.
!          As DLSODA switches automatically between stiff and nonstiff
!          methods, the required length of IWORK can change during
!          problem, between
!             LIS = 20 + NEQ   and   LIN = 20,
!          respectively.  Thus the IWORK array passed to DLSODA can
!          either have a fixed length of at least 20 + NEQ, or have a
!          dynami! length of at least LIN or LIS, depending on the
!          current method.  The comments on dynami! length under
!          RWORK above apply here.  Initially, this length need
!          only be at least LIN = 20.
!
!          The first few words of IWORK are used for conditional and
!          optional inputs and optional outputs.
!
!          The following 2 words in IWORK are conditional inputs:
!            IWORK(1) = ML     these are the lower and upper
!            IWORK(2) = MU     half-bandwidths, respectively, of the
!                       banded Jacobian, excluding the main diagonal.
!                       The band is defined by the matrix locations
!                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
!                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
!                       These are required if JT is 4 or 5, and
!                       ignored otherwise.  ML and MU may in fact be
!                       the band parameters for a matrix to which
!                       df/dy is only approximately equal.
!
! LIW    = the length of the array IWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! Note: The base addresses of the work arrays must not be
! altered between calls to DLSODA for the same problem.
! The contents of the work arrays must not be altered
! between calls, except possibly for the conditional and
! optional inputs, and except for the last 3*NEQ words of RWORK.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside DLSODA between calls, if
! desired (but not for use by F or JAC).
!
! JAC    = the name of the user-supplied routine to compute the
!          Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine
!          is optional, but if the problem is expected to be stiff much
!          of the time, you are encouraged to supply JAC, for the sake
!          of efficiency.  (Alternatively, set JT = 2 or 5 to have
!          DLSODA compute df/dy internally by difference quotients.)
!          If and when DLSODA uses df/dy, it treats this NEQ by NEQ
!          matrix either as full (JT = 1 or 2), or as banded (JT =
!          4 or 5) with half-bandwidths ML and MU (discussed under
!          IWORK above).  In either case, if JT = 1 or 4, the JAC
!          routine must compute df/dy as a function of the scalar t
!          and the vector y.  It is to have the form
!               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
!               DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!          PD is to be loaded with partial derivatives (elements of
!          the Jacobian matrix) on output.  PD must be given a first
!          dimension of NROWPD.  T and Y have the same meaning as in
!          Subroutine F.
!               In the full matrix case (JT = 1), ML and MU are
!          ignored, and the Jacobian is to be loaded into PD in
!          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!               In the band matrix case (JT = 4), the elements
!          within the band are to be loaded into PD in columnwise
!          manner, with diagonal lines of df/dy loaded into the rows
!          of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters (see IWORK).
!          The locations in PD in the two triangular areas which
!          correspond to nonexistent matrix elements can be ignored
!          or loaded arbitrarily, as they are overwritten by DLSODA.
!               JAC need not provide df/dy exactly.  A crude
!          approximation (possibly with a smaller bandwidth) will do.
!               In either case, PD is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Each call to JAC is preceded by a call to F with the same
!          arguments NEQ, T, and Y.  Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user Common block by F and not recomputed by JAC,
!          if desired.  Also, JAC may alter the Y array, if desired.
!          JAC must be declared External in the calling program.
!               Subroutine JAC may access user-defined quantities in
!          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
!          See the descriptions of NEQ and Y above.
!
! JT     = Jacobian type indicator.  Used only for input.
!          JT specifies how the Jacobian matrix df/dy will be
!          treated, if and when DLSODA requires this matrix.
!          JT has the following values and meanings:
!           1 means a user-supplied full (NEQ by NEQ) Jacobian.
!           2 means an internally generated (difference quotient) full
!             Jacobian (using NEQ extra calls to F per df/dy value).
!           4 means a user-supplied banded Jacobian.
!           5 means an internally generated banded Jacobian (using
!             ML+MU+1 extra calls to F per df/dy evaluation).
!          If JT = 1 or 4, the user must supply a Subroutine JAC
!          (the name is arbitrary) as described above under JAC.
!          If JT = 2 or 5, a dummy argument can be used.
!-----------------------------------------------------------------------
! Optional Inputs.
!
! The following is a list of the optional inputs provided for in the
! call sequence.  (See also Part 2.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of these inputs requires IOPT = 1, and in that
! case all of these inputs are examined.  A value of zero for any
! of these optional inputs will cause the default value to be used.
! Thus to use a subset of the optional inputs, simply preload
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! Name    Location      Meaning and Default Value
!
! H0      RWORK(5)  the step size to be attempted on the first step.
!                   The default value is determined by the solver.
!
! HMAX    RWORK(6)  the maximum absolute step size allowed.
!                   The default value is infinite.
!
! HMIN    RWORK(7)  the minimum absolute step size allowed.
!                   The default value is 0.  (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
!
! IXPR    IWORK(5)  flag to generate extra printing at method switches.
!                   IXPR = 0 means no extra printing (the default).
!                   IXPR = 1 means print data on each switch.
!                   T, H, and NST will be printed on the same logical
!                   unit as used for error messages.
!
! MXSTEP  IWORK(6)  maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 500.
!
! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value.  The default value is 10.
!
! MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff
!                   (Adams) method.  the default value is 12.
!                   if MXORDN exceeds the default value, it will
!                   be reduced to the default value.
!                   MXORDN is held constant during the problem.
!
! MXORDS  IWORK(9)  the maximum order to be allowed for the stiff
!                   (BDF) method.  The default value is 5.
!                   If MXORDS exceeds the default value, it will
!                   be reduced to the default value.
!                   MXORDS is held constant during the problem.
!-----------------------------------------------------------------------
! Optional Outputs.
!
! As optional additional output from DLSODA, the variables listed
! below are quantities related to the performance of DLSODA
! which are available to the user.  These are communicated by way of
! the work arrays, but also have internal mnemoni! names as shown.
! except where stated otherwise, all of these outputs are defined
! on any successful return from DLSODA, and on any return with
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, outputs relevant to the error will be defined,
! as noted below.
!
! Name    Location      Meaning
!
! HU      RWORK(11) the step size in t last used (successfully).
!
! HCUR    RWORK(12) the step size to be attempted on the next step.
!
! TCUR    RWORK(13) the current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  On output, TCUR
!                   will always be at least as far as the argument
!                   T, but may be farther (if interpolation was done).
!
! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise).  If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
!
! TSW     RWORK(15) the value of t at the time of the last method
!                   switch, if any.
!
! NST     IWORK(11) the number of steps taken for the problem so far.
!
! NFE     IWORK(12) the number of f evaluations for the problem so far.
!
! NJE     IWORK(13) the number of Jacobian evaluations (and of matrix
!                   LU decompositions) for the problem so far.
!
! NQU     IWORK(14) the method order last used (successfully).
!
! NQCUR   IWORK(15) the order to be attempted on the next step.
!
! IMXER   IWORK(16) the index of the component of largest magnitude in
!                   the weighted local error vector ( E(i)/EWT(i) ),
!                   on an error return with ISTATE = -4 or -5.
!
! LENRW   IWORK(17) the length of RWORK actually required, assuming
!                   that the length of RWORK is to be fixed for the
!                   rest of the problem, and that switching may occur.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENIW   IWORK(18) the length of IWORK actually required, assuming
!                   that the length of IWORK is to be fixed for the
!                   rest of the problem, and that switching may occur.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! MUSED   IWORK(19) the method indicator for the last successful step:
!                   1 means Adams (nonstiff), 2 means BDF (stiff).
!
! MCUR    IWORK(20) the current method indicator:
!                   1 means Adams (nonstiff), 2 means BDF (stiff).
!                   This is the method to be attempted
!                   on the next step.  Thus it differs from MUSED
!                   only if a method switch has just been made.
!
! The following two arrays are segments of the RWORK array which
! may also be of interest to the user as optional outputs.
! For each array, the table below gives its internal name,
! its base address in RWORK, and its description.
!
! Name    Base Address      Description
!
! YH      21             the Nordsieck history array, of size NYH by
!                        (NQCUR + 1), where NYH is the initial value
!                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
!                        of YH contains HCUR**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the solution,
!                        evaluated at T = TCUR.
!
! ACOR     LACOR         array of size NEQ used for the accumulated
!         (from Common   corrections on each step, scaled on output
!           as noted)    to represent the estimated local error in y
!                        on the last step.  This is the vector E in
!                        the description of the error control.  It is
!                        defined only on a successful return from
!                        DLSODA.  The base address LACOR is obtained by
!                        including in the user's program the
!                        following 2 lines:
!                           COMMON /DLS001/ RLS(218), ILS(37)
!                           LACOR = ILS(22)
!
!-----------------------------------------------------------------------
! Part 2.  Other Routines Callable.
!
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with DLSODA.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATE! error handling package.)
!
!     Form of Call                  Function
!   CALL XSETUN(LUN)          set the logical unit number, LUN, for
!                             output of messages from DLSODA, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!
!   CALL XSETF(MFLAG)         set a flag to control the printing of
!                             messages by DLSODA.
!                             MFLAG = 0 means do not print. (Danger:
!                             This risks losing valuable information.)
!                             MFLAG = 1 means print (the default).
!
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!
!   CALL DSRCMA(RSAV,ISAV,JOB) saves and restores the contents of
!                             the internal Common blocks used by
!                             DLSODA (see Part 3 below).
!                             RSAV must be a real array of length 240
!                             or more, and ISAV must be an integer
!                             array of length 46 or more.
!                             JOB=1 means save Common into RSAV/ISAV.
!                             JOB=2 means restore Common from RSAV/ISAV.
!                                DSRCMA is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with DLSODA.
!
!   CALL DINTDY(,,,,,)        provide derivatives of y, of various
!        (see below)          orders, at a specified point t, if
!                             desired.  It may be called only after
!                             a successful return from DLSODA.
!
! The detailed instructions for using DINTDY are as follows.
! The form of the call is:
!
!   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
!
! The input parameters are:
!
! T         = value of independent variable where answers are desired
!             (normally the same as the T last returned by DLSODA).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional outputs for TCUR and HU.)
! K         = integer order of the derivative desired.  K must satisfy
!             0 .le. K .le. NQCUR, where NQCUR is the current order
!             (see optional outputs).  The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by DLSODA directly.  Since NQCUR .ge. 1, the first
!             derivative dy/dt is always available with DINTDY.
! RWORK(21) = the base address of the history array YH.
! NYH       = column length of YH, equal to the initial value of NEQ.
!
! The output parameters are:
!
! DKY       = a real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
!-----------------------------------------------------------------------
! Part 3.  Common Blocks.
!
! If DLSODA is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in:
!   (1) the call sequence to DLSODA, and
!   (2) the two internal Common blocks
!         /DLS001/  of length  255  (218 double precision words
!                      followed by 37 integer words),
!         /DLSA01/  of length  31    (22 double precision words
!                      followed by  9 integer words).
!
! If DLSODA is used on a system in which the contents of internal
! Common blocks are not preserved between calls, the user should
! declare the above Common blocks in the calling program to insure
! that their contents are preserved.
!
! If the solution of a given problem by DLSODA is to be interrupted
! and then later continued, such as when restarting an interrupted run
! or alternating between two or more problems, the user should save,
! following the return from the last DLSODA call prior to the
! interruption, the contents of the call sequence variables and the
! internal Common blocks, and later restore these values before the
! next DLSODA call for that problem.  To save and restore the Common
! blocks, use Subroutine DSRCMA (see Part 2 above).
!
!-----------------------------------------------------------------------
! Part 4.  Optionally Replaceable Solver Routines.
!
! Below is a description of a routine in the DLSODA package which
! relates to the measurement of errors, and can be
! replaced by a user-supplied version, if desired.  However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note: The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) DEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above:
!     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODA call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by DEWSET.
!
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
! in y(i) to.  The EWT array returned by DEWSET is passed to the
! DMNORM routine, and also used by DLSODA in the computation
! of the optional output IMXER, and the increments for difference
! quotient Jacobians.
!
! In the user-supplied version of DEWSET, it may be desirable to use
! the current values of derivatives of y.  Derivatives up to order NQ
! are available from the history array YH, described above under
! optional outputs.  In DEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of H**j/factorial(j).  On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
! can be obtained by including in DEWSET the statements:
!     DOUBLE PRECISION RLS
!     COMMON /DLS001/ RLS(218),ILS(37)
!     NQ = ILS(33)
!     NST = ILS(34)
!     H = RLS(212)
! Thus, for example, the current value of dy/dt can be obtained as
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
! unnecessary when NST = 0).
!-----------------------------------------------------------------------
!     These subroutines and functions are now within the modules
!      EXTERNAL DPRJA, DSOLSY
!      DOUBLE PRECISION :: DUMACH, DMNORM
      INTEGER :: INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS, &
          ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER :: INSUFR, INSUFI, IXPR, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
      INTEGER :: I, I1, I2, IFLAG, IMXER, KGO, LF0, &
          LENIW, LENRW, LENWM, ML, MORD(2), MU, MXHNL0, MXSTP0
      INTEGER :: LEN1, LEN1C, LEN1N, LEN1S, LEN2, LENIWC, LENRWC
      DOUBLE PRECISION :: ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION :: TSW, ROWNS2, PDNORM
      DOUBLE PRECISION :: ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI, &
            TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      LOGICAL :: IHIT
      CHARACTER*60 MSG
      SAVE MORD, MXSTP0, MXHNL0
!-----------------------------------------------------------------------
! The following two internal Common blocks contain
! (a) variables which are local to any subroutine but whose values must
!     be preserved between calls to the routine ("own" variables), and
! (b) variables which are communicated between subroutines.
! The block DLS001 is declared in subroutines DLSODA, DINTDY, DSTODA,
! DPRJA, and DSOLSY.
! The block DLSA01 is declared in subroutines DLSODA, DSTODA, and DPRJA.
! Groups of variables are replaced by dummy arrays in the Common
! declarations in routines where those variables are not used.
!-----------------------------------------------------------------------
      COMMON /DLS001/ ROWNS(209), &
          CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
          INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6), &
          ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU

      COMMON /DLSA01/ TSW, ROWNS2(20), PDNORM, &
          INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS

      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .gt. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
!-----------------------------------------------------------------------

!     Initialization to avoid complains
      LENWM = 0
      RH = 0d0
      IHIT = .false.

      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all inputs and various initializations.
!
! First check legality of the non-optional inputs NEQ, ITOL, IOPT,
! JT, ML, and MU.
!-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      IF (JT .EQ. 3 .OR. JT .LT. 1 .OR. JT .GT. 5) GO TO 608
      JTYP = JT
      IF (JT .LE. 2) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
! Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      IXPR = 0
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      HMXI = 0.0D0
      HMIN = 0.0D0
      IF (ISTATE .NE. 1) GO TO 60
      H0 = 0.0D0
      MXORDN = MORD(1)
      MXORDS = MORD(2)
      GO TO 60
 40   IXPR = IWORK(5)
      IF (IXPR .LT. 0 .OR. IXPR .GT. 1) GO TO 611
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      MXORDN = IWORK(8)
      IF (MXORDN .LT. 0) GO TO 628
      IF (MXORDN .EQ. 0) MXORDN = 100
      MXORDN = MIN(MXORDN,MORD(1))
      MXORDS = IWORK(9)
      IF (MXORDS .LT. 0) GO TO 629
      IF (MXORDS .EQ. 0) MXORDS = 100
      MXORDS = MIN(MXORDS,MORD(2))
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0D0) GO TO 615
      HMXI = 0.0D0
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0D0) GO TO 616
!-----------------------------------------------------------------------
! Set work array pointers and check lengths LRW and LIW.
! If ISTATE = 1, METH is initialized to 1 here to facilitate the
! checking of work space lengths.
! Pointers to segments of RWORK and IWORK are named by prefixing L to
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
! If the lengths provided are insufficient for the current method,
! an error return occurs.  This is treated as illegal input on the
! first call, but as a problem interruption with ISTATE = -7 on a
! continuation call.  If the lengths are sufficient for the current
! method but not for both methods, a warning message is sent.
!-----------------------------------------------------------------------
 60   IF (ISTATE .EQ. 1) METH = 1
      IF (ISTATE .EQ. 1) NYH = N
      LYH = 21
      LEN1N = 20 + (MXORDN + 1)*NYH
      LEN1S = 20 + (MXORDS + 1)*NYH
      LWM = LEN1S + 1
      IF (JT .LE. 2) LENWM = N*N + 2
      IF (JT .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEN1S = LEN1S + LENWM
      LEN1C = LEN1N
      IF (METH .EQ. 2) LEN1C = LEN1S
      LEN1 = MAX(LEN1N,LEN1S)
      LEN2 = 3*N
      LENRW = LEN1 + LEN2
      LENRWC = LEN1C + LEN2
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      LENIWC = 20
      IF (METH .EQ. 2) LENIWC = LENIW
      IWORK(18) = LENIW
      IF (ISTATE .EQ. 1 .AND. LRW .LT. LENRWC) GO TO 617
      IF (ISTATE .EQ. 1 .AND. LIW .LT. LENIWC) GO TO 618
      IF (ISTATE .EQ. 3 .AND. LRW .LT. LENRWC) GO TO 550
      IF (ISTATE .EQ. 3 .AND. LIW .LT. LENIWC) GO TO 555
      LEWT = LEN1 + 1
      INSUFR = 0
      IF (LRW .GE. LENRW) GO TO 65
      INSUFR = 2
      LEWT = LEN1C + 1
      MSG='DLSODA-  Warning.. RWORK length is sufficient for now, but  '
      CALL XERRWD (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      may not be later.  Integration will proceed anyway.   '
      CALL XERRWD (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      Length needed is LENRW = I1, while LRW = I2.'
      CALL XERRWD (MSG, 50, 103, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
 65   LSAVF = LEWT + N
      LACOR = LSAVF + N
      INSUFI = 0
      IF (LIW .GE. LENIW) GO TO 70
      INSUFI = 2
      MSG='DLSODA-  Warning.. IWORK length is sufficient for now, but  '
      CALL XERRWD (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      may not be later.  Integration will proceed anyway.   '
      CALL XERRWD (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      Length needed is LENIW = I1, while LIW = I2.'
      CALL XERRWD (MSG, 50, 104, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
 70   CONTINUE
! Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 75 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0D0) GO TO 619
        IF (ATOLI .LT. 0.0D0) GO TO 620
 75     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
! If ISTATE = 3, set flag to signal parameter changes to DSTODA. -------
      JSTART = -1
      IF (N .EQ. NYH) GO TO 200
! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1, I2
 95     RWORK(I) = 0.0D0
      GO TO 200
!-----------------------------------------------------------------------
! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! and the calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.
!-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      TSW = T
      MAXORD = MXORDN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0) &
       H0 = TCRIT - T
 110  JSTART = 0
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0D0
      NQU = 0
      MUSED = 0
      MITER = 0
      CCMAX = 0.3D0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
! Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0D0
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
      IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
!-----------------------------------------------------------------------
! The coding below computes the step size, H0, to be attempted on the
! first step, unless the user has supplied a value for this.
! First check that TOUT - T differs significantly from zero.
! A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
! if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
! so as to be between 100*UROUND and 1.0E-3.
! Then the computed value H0 is given by:
!
!   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2
!
! where   w0     = MAX ( ABS(T), ABS(TOUT) ),
!         F      = the initial value of the vector f(t,y), and
!         norm() = the weighted vector norm used throughout, given by
!                  the DMNORM function routine, and weighted by the
!                  tolerances initially loaded into the EWT array.
! The sign of H0 is inferred from the initial values of TOUT and T.
! ABS(H0) is made .le. ABS(TOUT-T) in any case.
!-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0D0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
      IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
      AYI = ABS(Y(I))
      IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0D0*UROUND)
      TOL = MIN(TOL,0.001D0)
      SUM = DMNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0D0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0D0) H0 = H0/RH
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      GO TO 270
!-----------------------------------------------------------------------
! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      T = TN
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) T = TCRIT
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (ISTATE .EQ. 2 .AND. JSTART .GE. 0) JSTART = -2
!-----------------------------------------------------------------------
! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator DSTODA.
!
! This is a looping point for the integration steps.
!
! First check for too many steps being taken, update EWT (if not at
! start of problem), check for too much accuracy being requested, and
! check for H below the roundoff level in T.
!-----------------------------------------------------------------------
 250  CONTINUE
      IF (METH .EQ. MUSED) GO TO 255
      IF (INSUFR .EQ. 1) GO TO 550
      IF (INSUFI .EQ. 1) GO TO 555
 255  IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DMNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0D0) GO TO 280
      TOLSF = TOLSF*2.0D0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DLSODA-  Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DLSODA-  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
 290  CONTINUE
!-----------------------------------------------------------------------
!   CALL DSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPRJA,DSOLSY)
!-----------------------------------------------------------------------
      CALL DSTODA (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT), &
             RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM), &
             F, JAC, DPRJA, DSOLSY)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540), KGO
!-----------------------------------------------------------------------
! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0).
! If a method switch was just made, record TSW, reset MAXORD,
! set JSTART to -1 to signal DSTODA to complete the switch,
! and do extra printing of data if IXPR = 1.
! Then, in any case, check for stop conditions.
!-----------------------------------------------------------------------
 300  INIT = 1
      IF (METH .EQ. MUSED) GO TO 310
      TSW = TN
      MAXORD = MXORDN
      IF (METH .EQ. 2) MAXORD = MXORDS
      IF (METH .EQ. 2) RWORK(LWM) = SQRT(UROUND)
      INSUFR = MIN(INSUFR,1)
      INSUFI = MIN(INSUFI,1)
      JSTART = -1
      IF (IXPR .EQ. 0) GO TO 310
      IF (METH .EQ. 2) THEN
        MSG='DLSODA- A switch to the BDF (stiff) method has occurred     '
        CALL XERRWD (MSG, 60, 105, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      ENDIF
      IF (METH .EQ. 1) THEN
        MSG='DLSODA- A switch to the Adams (nonstiff) method has occurred'
        CALL XERRWD (MSG, 60, 106, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      ENDIF
      MSG='     at T = R1,  tentative step size H = R2,  step NST = I1 '
      CALL XERRWD (MSG, 60, 107, 0, 1, NST, 0, 2, TN, H)
 310  GO TO (320, 400, 330, 340, 350), ITASK
! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 320  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400
      GO TO 250
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (JSTART .GE. 0) JSTART = -2
      GO TO 250
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
!-----------------------------------------------------------------------
! Block G.
! The following block handles all successful returns from DLSODA.
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional outputs are loaded into the
! work arrays before returning.
!-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      RWORK(15) = TSW
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = MUSED
      IWORK(20) = METH
      RETURN
!-----------------------------------------------------------------------
! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input.  First the error message routine is called.
! If there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH and T is set to TN.
! The optional outputs are loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DLSODA-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)
      ISTATE = -1
      GO TO 580
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODA-  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
! Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DLSODA-  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DLSODA-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DLSODA-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
! RWORK length too small to proceed. -----------------------------------
 550  MSG = 'DLSODA-  At current T(=R1), RWORK length too small'
      CALL XERRWD (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      to proceed.  The integration was otherwise successful.'
      CALL XERRWD (MSG, 60, 206, 0, 0, 0, 0, 1, TN, 0.0D0)
      ISTATE = -7
      GO TO 580
! IWORK length too small to proceed. -----------------------------------
 555  MSG = 'DLSODA-  At current T(=R1), IWORK length too small'
      CALL XERRWD (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      to proceed.  The integration was otherwise successful.'
      CALL XERRWD (MSG, 60, 207, 0, 0, 0, 0, 1, TN, 0.0D0)
      ISTATE = -7
      GO TO 580
! Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0D0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
! Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      RWORK(15) = TSW
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = MUSED
      IWORK(20) = METH
      RETURN

!-----------------------------------------------------------------------
! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! First the error message routine is called.  If the illegal input
! is a negative ISTATE, the run is aborted (apparent infinite loop).
!-----------------------------------------------------------------------

 601  MSG = 'DLSODA-  ISTATE (=I1) illegal.'
      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DLSODA-  ITASK (=I1) illegal. '
      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 603  MSG = 'DLSODA-  ISTATE .gt. 1 but DLSODA not initialized.'
      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 604  MSG = 'DLSODA-  NEQ (=I1) .lt. 1     '
      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 605  MSG = 'DLSODA-  ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 606  MSG = 'DLSODA-  ITOL (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 607  MSG = 'DLSODA-  IOPT (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 608  MSG = 'DLSODA-  JT (=I1) illegal.    '
      CALL XERRWD (MSG, 30, 8, 0, 1, JT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 609  MSG = 'DLSODA-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '
      CALL XERRWD (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 610  MSG = 'DLSODA-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '
      CALL XERRWD (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 611  MSG = 'DLSODA-  IXPR (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 11, 0, 1, IXPR, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 612  MSG = 'DLSODA-  MXSTEP (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 613  MSG = 'DLSODA-  MXHNIL (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 614  MSG = 'DLSODA-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)
      GO TO 700
 615  MSG = 'DLSODA-  HMAX (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)
      GO TO 700
 616  MSG = 'DLSODA-  HMIN (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)
      GO TO 700
 617  MSG='DLSODA-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 618  MSG='DLSODA-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
      GO TO 700
 619  MSG = 'DLSODA-  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)
      GO TO 700
 620  MSG = 'DLSODA-  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODA-  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)
      GO TO 700
 622  MSG='DLSODA-  TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='DLSODA-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DLSODA-  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'DLSODA-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)
      GO TO 700
 628  MSG = 'DLSODA-  MXORDN (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 28, 0, 1, MXORDN, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 629  MSG = 'DLSODA-  MXORDS (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 29, 0, 1, MXORDS, 0, 0, 0.0D0, 0.0D0)

 700  ISTATE = -3
      RETURN

 800  MSG = 'DLSODA-  Run aborted.. apparent infinite loop.    '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
      RETURN
      END subroutine DLSODA


!----------------------- End of Subroutine DLSODA ----------------------

SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG)
!***BEGIN PROLOGUE  DINTDY
!***SUBSIDIARY
!***PURPOSE  Interpolate solution derivatives.
!***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DINTDY computes interpolated values of the K-th derivative of the
!  dependent variable vector y, and stores it in DKY.  This routine
!  is called within the package with K = 0 and T = TOUT, but may
!  also be called by the user for any K up to the current order.
!  (See detailed instructions in the usage documentation.)
!
!  The computed values in DKY are gotten by interpolation using the
!  Nordsieck history array YH.  This array corresponds uniquely to a
!  vector-valued polynomial of degree NQCUR or less, and DKY is set
!  to the K-th derivative of this polynomial at T.
!  The formula for DKY is:
!               q
!   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
!              j=K
!  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
!  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
!  communicated by COMMON.  The above sum is done in reverse order.
!  IFLAG is returned negative if either K or T is out of bounds.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  XERRWD
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDO! format.  (FNF)
!   890503  Minor cosmeti! changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!   050427  Corrected roundoff decrement in TP. (ACH)
!***END PROLOGUE  DINTDY
!**End
      INTEGER :: K, IFLAG, NYH
      DOUBLE PRECISION :: T, YH, DKY
      DIMENSION YH(NYH,*), DKY(*)
      INTEGER :: IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION :: ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
          IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER :: I, IC, J, JB, JB2, JJ, JJ1, JP1
      DOUBLE PRECISION :: C, R, S, TP
      CHARACTER*80 MSG

!***FIRST EXECUTABLE STATEMENT  DINTDY
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TP = TN - HU -  100.0D0*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
      IF ((T-TP)*(T-TN) .GT. 0.0D0) GO TO 90

      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO JJ = JJ1,NQ
        IC = IC*JJ
      ENDDO
 15   C = IC

      DO I = 1,N
        DKY(I) = C*YH(I,L)
      ENDDO
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
 30       IC = IC*JJ
 35     C = IC
        DO 40 I = 1,N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
 60     DKY(I) = R*DKY(I)
      RETURN

 80   MSG = 'DINTDY-  K (=I1) illegal      '
      CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0)
      IFLAG = -1
      RETURN
 90   MSG = 'DINTDY-  T (=R1) illegal      '
      CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
      END subroutine DINTDY
!----------------------- END OF SUBROUTINE DINTDY ----------------------

      SUBROUTINE DSOLSY (WM, IWM, X, TEM)
!***BEGIN PROLOGUE  DSOLSY
!***SUBSIDIARY
!***PURPOSE  ODEPACK linear system solver.
!***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This routine manages the solution of the linear system arising from
!  a chord iteration.  It is called if MITER .ne. 0.
!  If MITER is 1 or 2, it calls DGESL to accomplish this.
!  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
!  matrix, and then computes the solution.
!  If MITER is 4 or 5, it calls DGBSL.
!  Communication with DSOLSY uses the following variables:
!  WM    = real work space containing the inverse diagonal matrix if
!          MITER = 3 and the LU decomposition of the matrix otherwise.
!          Storage of matrix elements starts at WM(3).
!          WM also contains the following matrix-related data:
!          WM(1) = SQRT(UROUND) (not used here),
!          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
!  IWM   = integer work space containing pivot information, starting at
!          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
!          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
!  X     = the right-hand side vector on input, and the solution vector
!          on output, of length N.
!  TEM   = vector of work space of length N, not used in this version.
!  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
!          IERSL = 1 if a singular matrix arose with MITER = 3.
!  This routine also uses the COMMON variables EL0, H, MITER, and N.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  DGBSL, DGESL
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDO! format.  (FNF)
!   890503  Minor cosmeti! changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!***END PROLOGUE  DSOLSY
!**End
      INTEGER :: IWM(*)
      DOUBLE PRECISION :: WM(*), X(*), TEM(*)
      INTEGER :: IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION :: ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
          IOWND(6), IOWNS(6), ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER :: I, MEBAND, ML, MU
      DOUBLE PRECISION :: DI, HL0, PHL0, R

!***FIRST EXECUTABLE STATEMENT  DSOLSY
      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)
      RETURN

 300  PHL0 = WM(2)
      HL0 = H*EL0
      WM(2) = HL0
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
        IF (ABS(DI) .EQ. 0.0D0) GO TO 390
 320    WM(I+2) = 1.0D0/DI
 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN

 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)
      RETURN
      end subroutine DSOLSY
!----------------------- END OF SUBROUTINE DSOLSY ----------------------


      SUBROUTINE DSTODA (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR, &
          WM, IWM, F, JAC, PJAC, SLVS)
      implicit none

      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER :: NEQ(*), IWM(*), NYH
      DOUBLE PRECISION :: Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*), ACOR(*), WM(*)
!      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*), &
!       ACOR(*), WM(*), IWM(*)
      INTEGER :: IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, &
          ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER :: IOWND2, ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
      DOUBLE PRECISION :: CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO, &
          CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION :: ROWND2, CM1, CM2, PDEST, PDLAST, RATIO, PDNORM
      COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12), &
          HOLD, RMAX, TESCO(3,12), &
          CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
          IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP, &
          ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSA01/ ROWND2, CM1(12), CM2(5), PDEST, PDLAST, RATIO, &
          PDNORM, &
          IOWND2(3), ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
      INTEGER :: I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
      INTEGER :: LM1, LM1P1, LM2, LM2P1, NQM1, NQM2
      DOUBLE PRECISION :: DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP, &
          R, RH, RHDN, RHSM, RHUP, TOLD  !DMNORM
      DOUBLE PRECISION :: ALPHA, DM1,DM2, EXM1,EXM2, &
          PDH, PNORM, RATE, RH1, RH1IT, RH2, RM, SM1(12)
      SAVE SM1
      DATA SM1/0.5D0, 0.575D0, 0.55D0, 0.45D0, 0.35D0, 0.25D0, &
        0.20D0, 0.15D0, 0.10D0, 0.075D0, 0.050D0, 0.025D0/
!-----------------------------------------------------------------------
! DSTODA performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! Note: DSTODA is independent of the value of the iteration method
! indicator MITER, when this is .ne. 0, and hence is independent
! of the type of chord method used, or the Jacobian structure.
! Communication with DSTODA is done with the following variables:
!
! Y      = an array of length .ge. N used as the Y argument in
!          all calls to F and JAC.
! NEQ    = integer array containing problem size in NEQ(1), and
!          passed as the NEQ argument in all calls to F and JAC.
! YH     = an NYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! NYH    = a constant integer .ge. N, the first dimension of YH.
! YH1    = a one-dimensional array occupying the same space as YH.
! EWT    = an array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = an array of working storage, of length N.
! ACOR   = a work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = real and integer work arrays associated with matrix
!          operations in chord iteration (MITER .ne. 0).
! PJAC   = name of routine to evaluate and preprocess Jacobian matrix
!          and P = I - H*EL0*Jac, if a chord method is being used.
!          It also returns an estimate of norm(Jac) in PDNORM.
! SLVS   = name of routine to solve linear system in chord iteration.
! CCMAX  = maximum relative change in H*EL0 before PJAC is called.
! H      = the step size to be attempted on the next step.
!          H is altered by the error control algorithm during the
!          problem.  H can be either positive or negative, but its
!          sign must remain constant throughout the problem.
! HMIN   = the minimum absolute value of the step size H to be used.
! HMXI   = inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
!          HMIN and HMXI may be changed at any time, but will not
!          take effect until the next change of H is considered.
! TN     = the independent variable. TN is updated on each step taken.
! JSTART = an integer used for input only, with the following
!          values and meanings:
!               0  perform the first step.
!           .gt.0  take a new step continuing from the last.
!              -1  take the next step with a new value of H,
!                    N, METH, MITER, and/or matrix parameters.
!              -2  take the next step with a new value of H,
!                    but with other inputs unchanged.
!          On return, JSTART is set to 1 to facilitate continuation.
! KFLAG  = a completion code with the following meanings:
!               0  the step was succesful.
!              -1  the requested error could not be achieved.
!              -2  corrector convergence could not be achieved.
!              -3  fatal error in PJAC or SLVS.
!          A return with KFLAG = -1 or -2 means either
!          ABS(H) = HMIN or 10 consecutive failures occurred.
!          On a return with KFLAG negative, the values of TN and
!          the YH array are as of the beginning of the last
!          step, and H is the last step size attempted.
! MAXORD = the maximum order of integration method to be allowed.
! MAXCOR = the maximum number of corrector iterations allowed.
! MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
! MXNCF  = maximum number of convergence failures allowed.
! METH   = current method.
!          METH = 1 means Adams method (nonstiff)
!          METH = 2 means BDF method (stiff)
!          METH may be reset by DSTODA.
! MITER  = corrector iteration method.
!          MITER = 0 means functional iteration.
!          MITER = JT .gt. 0 means a chord iteration corresponding
!          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
!          communicated here as JTYP, but is not used in DSTODA
!          except to load MITER following a method switch.)
!          MITER may be reset by DSTODA.
! N      = the number of first-order differential equations.
!-----------------------------------------------------------------------
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      DSM = 0d0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  RMAX is the maximum ratio by which H can be increased
! in a single step.  It is initially 1.E4 to compensate for the small
! initial H, but then is normally equal to 10.  If a failure
! occurs (in corrector convergence or error test), RMAX is set at 2
! for the next increase.
! DCFODE is called to get the needed coefficients for both methods.
!-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      NSLP = 0
      IPUP = MITER
      IRET = 3
      IREDO = 0
! Initialize switching parameters.  METH = 1 is assumed initially. -----
      ICOUNT = 20
      IRFLAG = 0
      PDEST = 0.0D0
      PDLAST = 0.0D0
      RATIO = 5.0D0
      CALL DCFODE (2, ELCO, TESCO)
      DO 10 I = 1,5
 10     CM2(I) = TESCO(2,I)*ELCO(I+1,I)
      CALL DCFODE (1, ELCO, TESCO)
      DO 20 I = 1,12
 20     CM1(I) = TESCO(2,I)*ELCO(I+1,I)
      GO TO 150
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! IPUP is set to MITER to force a matrix update.
! If an order increase is about to be considered (IALTH = 1),
! IALTH is reset to 2 to postpone consideration one more step.
! If the caller has changed METH, DCFODE is called to reset
! the coefficients of the method.
! If H is to be changed, YH must be rescaled.
! If H or METH is being changed, IALTH is reset to L = NQ + 1
! to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MUSED) GO TO 160
      CALL DCFODE (METH, ELCO, TESCO)
      IALTH = L
      IRET = 1
!-----------------------------------------------------------------------
! The el vector and related constants are reset
! whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      GO TO (160, 170, 200), IRET
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
! L = NQ + 1 to prevent a change of H for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
!-----------------------------------------------------------------------
! If METH = 1, also restrict the new step size by the stability region.
! If this reduces H, set IRFLAG to 1 so that if there are roundoff
! problems later, we can assume that is the cause of the trouble.
!-----------------------------------------------------------------------
      IF (METH .EQ. 2) GO TO 178
      IRFLAG = 0
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (RH*PDH*1.00001D0 .LT. SM1(NQ)) GO TO 178
      RH = SM1(NQ)/PDH
      IRFLAG = 1
 178  CONTINUE
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! R! is the ratio of new to old values of the coefficient  H*EL(1).
! When R! differs from 1 by more than CCMAX, IPUP is set to MITER
! to force PJAC to be called, if a Jacobian is involved.
! In any case, PJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
! DIR$ IVDEP
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
      PNORM = DMNORM (N, YH1, EWT)
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the RMS-norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 220  M = 0
      RATE = 0.0D0
      DEL = 0.0D0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
 270  IF (MITER .NE. 0) GO TO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DMNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.
!-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL SLVS (WM, IWM, Y, SAVF)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DMNORM (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!
! We first check for a change of iterates that is the size of
! roundoff error.  If this occurs, the iteration has converged, and a
! new rate estimate is not formed.
! In all other cases, force at least two iterations to estimate a
! local Lipschitz constant estimate for Adams methods.
! On convergence, form PDEST = local maximum Lipschitz constant
! estimate.  PDLAST is the most recent nonzero estimate.
!-----------------------------------------------------------------------
 400  CONTINUE
      IF (DEL .LE. 100.0D0*PNORM*UROUND) GO TO 450
      IF (M .EQ. 0 .AND. METH .EQ. 1) GO TO 405
      IF (M .EQ. 0) GO TO 402
      RM = 1024.0D0
      IF (DEL .LE. 1024.0D0*DELP) RM = DEL/DELP
      RATE = MAX(RATE,RM)
      CRATE = MAX(0.2D0*CRATE,RM)
 402  DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .GT. 1.0D0) GO TO 405
      PDEST = MAX(PDEST,RATE/ABS(H*EL(1)))
      IF (PDEST .NE. 0.0D0) PDLAST = PDEST
      GO TO 450
 405  CONTINUE
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
!-----------------------------------------------------------------------
! The corrector iteration failed to converge.
! If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
! the next try.  Otherwise the YH array is retracted to its values
! before prediction, and H is reduced, if possible.  If H cannot be
! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
!-----------------------------------------------------------------------
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
! DIR$ IVDEP
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.25D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0
! to signal that the Jacobian involved may need updating later.
! The local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DMNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH array.
! Decrease ICOUNT by 1, and if it is -1, consider switching methods.
! If a method switch is made, reset various parameters,
! rescale the YH array, and exit.  If there is no switch,
! consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
! use in a possible order increase on the next step.
! If a change in H is considered, an increase or decrease in order
! by one is considered also.  A change in H is made only if it is by a
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      MUSED = METH
      DO 460 J = 1,L
        DO 460 I = 1,N
 460      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      ICOUNT = ICOUNT - 1
      IF (ICOUNT .GE. 0) GO TO 488
      IF (METH .EQ. 2) GO TO 480
!-----------------------------------------------------------------------
! We are currently using an Adams method.  Consider switching to BDF.
! If the current order is greater than 5, assume the problem is
! not stiff, and skip this section.
! If the Lipschitz constant and error estimate are not polluted
! by roundoff, go to 470 and perform the usual test.
! Otherwise, switch to the BDF methods if the last step was
! restricted to insure stability (irflag = 1), and stay with Adams
! method if not.  When switching to BDF with polluted error estimates,
! in the absence of other information, double the step size.
!
! When the estimates are OK, we make the usual test by computing
! the step size we could have (ideally) used on this step,
! with the current (Adams) method, and also that for the BDF.
! If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
! Compare the two step sizes to decide whether to switch.
! The step size advantage must be at least RATIO = 5 to switch.
!-----------------------------------------------------------------------
      IF (NQ .GT. 5) GO TO 488
      IF (DSM .GT. 100.0D0*PNORM*UROUND .AND. PDEST .NE. 0.0D0) &
        GO TO 470
      IF (IRFLAG .EQ. 0) GO TO 488
      RH2 = 2.0D0
      NQM2 = MIN(NQ,MXORDS)
      GO TO 478
 470  CONTINUE
      EXSM = 1.0D0/L
      RH1 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RH1IT = 2.0D0*RH1
      PDH = PDLAST*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQ)/PDH
      RH1 = MIN(RH1,RH1IT)
      IF (NQ .LE. MXORDS) GO TO 474
         NQM2 = MXORDS
         LM2 = MXORDS + 1
         EXM2 = 1.0D0/LM2
         LM2P1 = LM2 + 1
         DM2 = DMNORM (N, YH(1,LM2P1), EWT)/CM2(MXORDS)
         RH2 = 1.0D0/(1.2D0*DM2**EXM2 + 0.0000012D0)
         GO TO 476
 474  DM2 = DSM*(CM1(NQ)/CM2(NQ))
      RH2 = 1.0D0/(1.2D0*DM2**EXSM + 0.0000012D0)
      NQM2 = NQ
 476  CONTINUE
      IF (RH2 .LT. RATIO*RH1) GO TO 488
! THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
 478  RH = RH2
      ICOUNT = 20
      METH = 2
      MITER = JTYP
      PDLAST = 0.0D0
      NQ = NQM2
      L = NQ + 1
      GO TO 170
!-----------------------------------------------------------------------
! We are currently using a BDF method.  Consider switching to Adams.
! Compute the step size we could have (ideally) used on this step,
! with the current (BDF) method, and also that for the Adams.
! If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
! Compare the two step sizes to decide whether to switch.
! The step size advantage must be at least 5/RATIO = 1 to switch.
! If the step size for Adams would be so small as to cause
! roundoff pollution, we stay with BDF.
!-----------------------------------------------------------------------
 480  CONTINUE
      EXSM = 1.0D0/L
      IF (MXORDN .GE. NQ) GO TO 484
        NQM1 = MXORDN
        LM1 = MXORDN + 1
        EXM1 = 1.0D0/LM1
        LM1P1 = LM1 + 1
        DM1 = DMNORM (N, YH(1,LM1P1), EWT)/CM1(MXORDN)
        RH1 = 1.0D0/(1.2D0*DM1**EXM1 + 0.0000012D0)
        GO TO 486
 484  DM1 = DSM*(CM2(NQ)/CM1(NQ))
      RH1 = 1.0D0/(1.2D0*DM1**EXSM + 0.0000012D0)
      NQM1 = NQ
      EXM1 = EXSM
 486  RH1IT = 2.0D0*RH1
      PDH = PDNORM*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQM1)/PDH
      RH1 = MIN(RH1,RH1IT)
      RH2 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      IF (RH1*RATIO .LT. 5.0D0*RH2) GO TO 488
      ALPHA = MAX(0.001D0,RH1)
      DM1 = (ALPHA**EXM1)*DM1
      IF (DM1 .LE. 1000.0D0*UROUND*PNORM) GO TO 488
! The switch test passed.  Reset relevant quantities for Adams. --------
      RH = RH1
      ICOUNT = 20
      METH = 1
      MITER = 0
      PDLAST = 0.0D0
      NQ = NQM1
      L = NQ + 1
      GO TO 170

! No method switch is being made.  Do the usual step/order selection. --
 488  CONTINUE
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for this or
! one lower order.  After 2 or more failures, H is forced to decrease
! by a factor of 0.2 or less.
!-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
!DIR$ IVDEP
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
! at order NQ - 1, order NQ, or order NQ + 1, respectively.
! In the case of failure, RHUP = 0.0 to avoid an order increase.
! The largest of these is determined and the new order chosen
! accordingly.  If the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DMNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 550
      DDN = DMNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
! If METH = 1, limit RH according to the stability region also. --------
 550  IF (METH .EQ. 2) GO TO 560
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (L .LT. LMAX) RHUP = MIN(RHUP,SM1(L)/PDH)
      RHSM = MIN(RHSM,SM1(NQ)/PDH)
      IF (NQ .GT. 1) RHDN = MIN(RHDN,SM1(NQ-1)/PDH)
      PDEST = 0.0D0
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
! If METH = 1 and H is restricted by stability, bypass 10 percent test.
 620  IF (METH .EQ. 2) GO TO 622
      IF (RH*PDH*1.00001D0 .GE. SM1(NEWQ)) GO TO 625
 622  IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 610
 625  IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, L, and the coefficients.
! In any case H is reset according to RH and the YH array is rescaled.
! Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more failures have occured.
! If 10 failures have occurred, exit with KFLAG = -1.
! It is assumed that the derivatives that have accumulated in the
! YH array have errors of the wrong order.  Hence the first
! derivative is recomputed, and the order is set to 1.  Then
! H is reduced by a factor of 10, and the step is retried,
! until it succeeds or H reaches HMIN.
!-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD
! to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
!----------------------- End of Subroutine DSTODA ----------------------
      END subroutine DSTODA


      SUBROUTINE DPRJA (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, F, JAC)
      EXTERNAL F, JAC
      INTEGER :: NEQ, IWM, NYH
      DOUBLE PRECISION :: Y, YH, EWT, FTEM, SAVF, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*), WM(*), IWM(*)
      INTEGER :: IOWND, IOWNS, ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER :: IOWND2, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
      DOUBLE PRECISION :: ROWNS, CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION :: ROWND2, ROWNS2, PDNORM
      COMMON /DLS001/ ROWNS(209), &
          CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND, &
          IOWND(6), IOWNS(6), &
          ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, &
          LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER, &
          MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSA01/ ROWND2, ROWNS2(20), PDNORM, &
          IOWND2(3), IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
      INTEGER :: I, I1, I2, IER, II, J, J1, JJ, LENP, &
          MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
      DOUBLE PRECISION :: CON, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ
!         Comment out the function names
!          DMNORM, DFNORM, DBNORM
!-----------------------------------------------------------------------
! DPRJA is called by DSTODA to compute and process the matrix
! P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
! Here J is computed by the user-supplied routine JAC if
! MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
! J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
! matrix norm consistent with the weighted max-norm on vectors given
! by DMNORM) is computed, and J is overwritten by P.  P is then
! subjected to LU decomposition in preparation for later solution
! of linear systems with P as coefficient matrix.  This is done
! by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
!
! In addition to variables described previously, communication
! with DPRJA uses the following:
! Y     = array containing predicted values on entry.
! FTEM  = work array of length N (ACOR in DSTODA).
! SAVF  = array containing f evaluated at predicted y.
! WM    = real work space for matrices.  On output it contains the
!         LU decomposition of P.
!         Storage of matrix elements starts at WM(3).
!         WM also contains the following matrix-related data:
!         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
! IWM   = integer work space containing pivot information, starting at
!         IWM(21).   IWM also contains the band parameters
!         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
! EL0   = EL(1) (input).
! PDNORM= norm of Jacobian matrix. (Output).
! IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
!         P matrix found to be singular.
! JCUR  = output flag = 1 to indicate that the Jacobian matrix
!         (or approximation) is now current.
! This routine also uses the Common variables EL0, H, TN, UROUND,
! MITER, N, NFE, and NJE.
!-----------------------------------------------------------------------
      NJE = NJE + 1
      IERPJ = 0
      JCUR = 1
      HL0 = H*EL0
      GO TO (100, 200, 300, 400, 500), MITER
! If MITER = 1, call JAC and multiply by scalar. -----------------------
 100  LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 240
! If MITER = 2, make N calls to F to approximate J. --------------------
 200  FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL F (NEQ, TN, Y, FTEM)
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
 240  CONTINUE
! Compute norm of Jacobian. --------------------------------------------
      PDNORM = DFNORM (N, WM(3), EWT)/ABS(HL0)
! Add identity matrix. -------------------------------------------------
      J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + 1.0D0
 250    J = J + NP1
! Do LU decomposition on P. --------------------------------------------
      CALL DGEFA(WM(3), N, N, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
! Dummy block only, since MITER is never 3 in this routine. ------------
 300  RETURN
! If MITER = 4, call JAC and multiply by scalar. -----------------------
 400  ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
 410    WM(I+2) = 0.0D0
      CALL JAC(NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
      CON = -HL0
      DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
      GO TO 570
! If MITER = 5, make MBAND calls to F to approximate J. ----------------
 500  ML = IWM(1)
      MU = IWM(2)
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (NEQ, TN, Y, FTEM)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
 570  CONTINUE
! Compute norm of Jacobian. --------------------------------------------
      PDNORM = DBNORM (N, WM(ML+3), MEBAND, ML, MU, EWT)/ABS(HL0)
! Add identity matrix. -------------------------------------------------
      II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + 1.0D0
 580    II = II + MEBAND
! Do LU decomposition of P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
      end subroutine DPRJA
!----------------------- End of Subroutine DPRJA -----------------------

SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
!***BEGIN PROLOGUE  XERRWD
!***SUBSIDIARY
!***PURPOSE  Write error message with values.
!***CATEGORY  R3C
!***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
!  as given here, constitute a simplified version of the SLATEC error
!  handling package.
!
!  All arguments are input arguments.
!
!  MSG    = The message (character array).
!  NMES   = The length of MSG (number of characters).
!  NERR   = The error number (not used).
!  LEVEL  = The error level..
!           0 or 1 means recoverable (control returns to caller).
!           2 means fatal (run is aborted--see note below).
!  NI     = Number of integers (0, 1, or 2) to be printed with message.
!  I1,I2  = Integers to be printed, depending on NI.
!  NR     = Number of reals (0, 1, or 2) to be printed with message.
!  R1,R2  = Reals to be printed, depending on NR.
!
!  Note..  this routine is machine-dependent and specialized for use
!  in limited context, in the following ways..
!  1. The argument MSG is assumed to be of type CHARACTER, and
!     the message is printed with a format of (1X,A).
!  2. The message is assumed to take only one line.
!     Multi-line messages are generated by repeated calls.
!  3. If LEVEL = 2, control passes to the statement   STOP
!     to abort the run.  This statement may be machine-dependent.
!  4. R1 and R2 are assumed to be in double precision and are printed
!     in D21.13 format.
!
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   920831  DATE WRITTEN
!   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
!   930329  Modified prologue to SLATEC format. (FNF)
!   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
!   930922  Minor cosmetic change. (FNF)
!***END PROLOGUE  XERRWD
!
!*Internal Notes:
!
! For a different default logical unit number, IXSAV (or a subsidiary
! routine that it calls) will need to be modified.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERRWD.. None
! Function routine called by XERRWD.. IXSAV
!-----------------------------------------------------------------------
!**End
!
!  Declare arguments.

      DOUBLE PRECISION R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
      CHARACTER*(*) MSG

!  Declare local variables.

      INTEGER LUNIT, MESFLG !IXSAV

!  Get logical unit number and message print flag.

!***FIRST EXECUTABLE STATEMENT  XERRWD
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100

!  Write the message.

      WRITE (LUNIT,10)  MSG
 10   FORMAT(1X,A)
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',D21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)

!  Abort the run if LEVEL = 2.

 100  IF (LEVEL .NE. 2) RETURN
      STOP

      END subroutine XERRWD
!----------------------- End of Subroutine XERRWD ----------------------


      SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
!***BEGIN PROLOGUE  DGBFA
!***PURPOSE  Factor a band matrix using Gaussian elimination.
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGBFA factors a double precision band matrix by elimination.
!
!     DGBFA is usually called by DGBCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. 2*ML + MU + 1 .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 .LE. ML .LT.  N .
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 .LE. MU .LT.  N .
!                More efficient if  ML .LE. MU .
!     On Return
!
!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGBSL will divide by zero if
!                     called.  Use  RCOND  in DGBCO for a reliable
!                     indication of singularity.
!
!     Band Storage
!
!           If  A  is a band matrix, the following program segment
!           will set up the input.
!
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGBFA
      INTEGER :: LDA,N,ML,MU,IPVT(*),INFO
      DOUBLE PRECISION :: ABD(LDA, *)
!
      DOUBLE PRECISION :: T
      INTEGER :: I,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1 !IDAMAX
!
!***FIRST EXECUTABLE STATEMENT  DGBFA
      M = ML + MU + 1
      INFO = 0
!
!     ZERO INITIAL FILL-IN COLUMNS
!
      J0 = MU + 2
      J1 = MIN(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
        I0 = M + 1 - JZ
        DO 10 I = I0, ML
          ABD(I,JZ) = 0.0D0
   10   CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
        KP1 = K + 1
!
!        ZERO NEXT FILL-IN COLUMN
!
        JZ = JZ + 1
        IF (JZ .GT. N) GO TO 50
        IF (ML .LT. 1) GO TO 50
          DO 40 I = 1, ML
            ABD(I,JZ) = 0.0D0
   40     CONTINUE
   50   CONTINUE
!
!        FIND L = PIVOT INDEX
!
        LM = MIN(ML,N-K)
        L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
        IPVT(K) = L + K - M
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
        IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
!
!           INTERCHANGE IF NECESSARY
!
          IF (L .EQ. M) GO TO 60
            T = ABD(L,K)
            ABD(L,K) = ABD(M,K)
            ABD(M,K) = T
   60     CONTINUE
!
!           COMPUTE MULTIPLIERS
!
          T = -1.0D0/ABD(M,K)
          CALL DSCAL(LM,T,ABD(M+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
          JU = MIN(MAX(JU,MU+IPVT(K)),N)
          MM = M
          IF (JU .LT. KP1) GO TO 90
          DO 80 J = KP1, JU
            L = L - 1
            MM = MM - 1
            T = ABD(L,J)
            IF (L .EQ. MM) GO TO 70
              ABD(L,J) = ABD(MM,J)
              ABD(MM,J) = T
   70       CONTINUE
            CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80     CONTINUE
   90     CONTINUE
        GO TO 110
  100   CONTINUE
          INFO = K
  110   CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END


      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
!   BEGIN PROLOGUE  DGEFA
!   PURPOSE  Factor a matrix using Gaussian elimination.
!   CATEGORY  D2A1
!   TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
!   KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!   AUTHOR  Moler, C. B., (U. of New Mexico)
!   DESCRIPTION
!
!     DGEFA factors a double precision matrix by Gaussian elimination.
!
!     DGEFA is usually called by DGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGESL or DGEDI will divide by zero
!                     if called.  Use  RCOND  in DGECO for a reliable
!                     indication of singularity.
!
!   REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!   ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
!   REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)

      DOUBLE PRECISION T
      INTEGER J,K,KP1,L,NM1 !IDAMAX

    !     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

    !   FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
        KP1 = K + 1

!            FIND L = PIVOT INDEX

        L = IDAMAX(N-K+1,A(K,K),1) + K - 1
        IPVT(K) = L

    !        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED

        IF (A(L,K) .EQ. 0.0D0) GO TO 40

    !           INTERCHANGE IF NECESSARY

          IF (L .EQ. K) GO TO 10
             T = A(L,K)
             A(L,K) = A(K,K)
             A(K,K) = T
 10       CONTINUE

    !           COMPUTE MULTIPLIERS

          T = -1.0D0/A(K,K)
          CALL DSCAL(N-K,T,A(K+1,K),1)

    !           ROW ELIMINATION WITH COLUMN INDEXING

          DO 30 J = KP1, N
            T = A(L,J)
            IF (L .EQ. K) GO TO 20
              A(L,J) = A(K,J)
              A(K,J) = T
   20       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30     CONTINUE
        GO TO 50
   40   CONTINUE
          INFO = K
   50   CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END subroutine DGEFA


      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
!***BEGIN PROLOGUE  DGESL
!***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
!            factors computed by DGECO or DGEFA.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGESL solves the double precision system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DGECO or DGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGECO has set RCOND .GT. 0.0
!        or DGEFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGESL
      INTEGER :: LDA, N, IPVT(*), JOB
      DOUBLE PRECISION :: A(LDA, *), B(*)

      DOUBLE PRECISION :: DDOT,T
      INTEGER :: K,KB,L,NM1
!***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50

!            JOB = 0 , SOLVE  A * X = B
!            FIRST SOLVE  L*Y = B

        IF (NM1 .LT. 1) GO TO 30
        DO 20 K = 1, NM1
          L = IPVT(K)
          T = B(L)
          IF (L .EQ. K) GO TO 10
            B(L) = B(K)
            B(K) = T
   10     CONTINUE
          CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20   CONTINUE
   30   CONTINUE

!            NOW SOLVE  U*X = Y

        DO 40 KB = 1, N
          K = N + 1 - KB
          B(K) = B(K)/A(K,K)
          T = -B(K)
          CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40   CONTINUE
      GO TO 100
   50 CONTINUE

!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B

        DO 60 K = 1, N
          T = DDOT(K-1,A(1,K),1,B(1),1)
          B(K) = (B(K) - T)/A(K,K)
   60   CONTINUE

    !        NOW SOLVE TRANS(L)*X = Y

        IF (NM1 .LT. 1) GO TO 90
        DO 80 KB = 1, NM1
          K = N - KB
          B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
          L = IPVT(K)
          IF (L .EQ. K) GO TO 70
             T = B(L)
             B(L) = B(K)
             B(K) = T
   70     CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END subroutine DGESL


      SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
!***BEGIN PROLOGUE  DGBSL
!***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
!            the factors computed by DGBCO or DGBFA.
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGBSL solves the double precision band system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGBCO or DGBFA.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                the output from DGBCO or DGBFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGBCO or DGBFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B , where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGBCO has set RCOND .GT. 0.0
!        or DGBFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGBSL
      INTEGER :: LDA, N, ML, MU, JOB
      integer :: IPVT(*)
      DOUBLE PRECISION :: ABD(LDA, *), B(*)

      DOUBLE PRECISION :: DDOT, T
      INTEGER :: K, KB, L, LA, LB, LM, M, NM1
!***FIRST EXECUTABLE STATEMENT  DGBSL
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50

    !        JOB = 0 , SOLVE  A * X = B
    !        FIRST SOLVE L*Y = B

        IF (ML .EQ. 0) GO TO 30
        IF (NM1 .LT. 1) GO TO 30
          DO 20 K = 1, NM1
            LM = MIN(ML,N-K)
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
              B(L) = B(K)
              B(K) = T
   10       CONTINUE
            CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20     CONTINUE
   30   CONTINUE

    !        NOW SOLVE  U*X = Y

        DO 40 KB = 1, N
          K = N + 1 - KB
          B(K) = B(K)/ABD(M,K)
          LM = MIN(K,M) - 1
          LA = M - LM
          LB = K - LM
          T = -B(K)
          CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40   CONTINUE
      GO TO 100
   50 CONTINUE

    !        JOB = NONZERO, SOLVE  TRANS(A) * X = B
    !        FIRST SOLVE  TRANS(U)*Y = B

        DO 60 K = 1, N
          LM = MIN(K,M) - 1
          LA = M - LM
          LB = K - LM
          T = DDOT(LM,ABD(LA,K),1,B(LB),1)
          B(K) = (B(K) - T)/ABD(M,K)
   60   CONTINUE

    !        NOW SOLVE TRANS(L)*X = Y

        IF (ML .EQ. 0) GO TO 90
        IF (NM1 .LT. 1) GO TO 90
          DO 80 KB = 1, NM1
            K = N - KB
            LM = MIN(ML,N-K)
            B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
              T = B(L)
              B(L) = B(K)
              B(K) = T
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
      RETURN
      END subroutine DGBSL


      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
!***BEGIN PROLOGUE  DAXPY
!***PURPOSE  Compute a constant times a vector plus a vector.
!***CATEGORY  D1A7
!***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scalar multiplier
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  double precision result (unchanged if N .LE. 0)
!
!     Overwrite double precision DY with double precision DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DAXPY
      integer, intent(in) :: N, INCY, INCX
!     The variables originally not declared:
      integer :: M, MP1, NS, NX, I, IY, IX
      DOUBLE PRECISION :: DX(N), DY(N)
      double precision :: DA
!***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60

!     Code for unequal or nonpositive increments.

    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

!     Code for both increments equal to 1.

!     Clean-up loop so remaining vector length is a multiple of 4.

   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN

!     Code for equal, positive, non-unit increments.

   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END subroutine DAXPY



      SUBROUTINE DSCAL (N, DA, DX, INCX)
!***BEGIN PROLOGUE  DSCAL
!***PURPOSE  Multiply a vector by a constant.
!***CATEGORY  D1A6
!***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scale factor
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!       DX  double precision result (unchanged if N.LE.0)
!
!     Replace double precision DX by double precision DA*DX.
!     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSCAL
      DOUBLE PRECISION :: DA, DX(*)
      INTEGER :: I, INCX, IX, M, MP1, N
!***FIRST EXECUTABLE STATEMENT  DSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20

!     Code for increment not equal to 1.

      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN

!     Code for increment equal to 1.

!     Clean-up loop so remaining vector length is a multiple of 5.

   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END subroutine DSCAL


      SUBROUTINE DCFODE (METH, ELCO, TESCO)
!***BEGIN PROLOGUE  DCFODE
!***SUBSIDIARY
!***PURPOSE  Set ODE integrator coefficients.
!***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DCFODE is called by the integrator routine to set coefficients
!  needed there.  The coefficients for the current method, as
!  given by the value of METH, are set for all orders and saved.
!  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
!  (A smaller value of the maximum order is also allowed.)
!  DCFODE is called once at the beginning of the problem,
!  and is not called again unless and until METH is changed.
!
!  The ELCO array contains the basic method coefficients.
!  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
!  order nq are stored in ELCO(i,nq).  They are given by a genetrating
!  polynomial, i.e.,
!      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
!  For the implicit Adams methods, l(x) is given by
!      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
!  For the BDF methods, l(x) is given by
!      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
!  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!
!  The TESCO array contains test constants used for the
!  local error test and the selection of step size and/or order.
!  At order nq, TESCO(k,nq) is used for the selection of step
!  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
!  nq + 1 if k = 3.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DCFODE
!**End
      INTEGER :: METH
      INTEGER :: I, IB, NQ, NQM1, NQP1
      DOUBLE PRECISION :: ELCO(13, 12), TESCO(3, 12), PC(12)
      DOUBLE PRECISION :: AGAMQ, FNQ, FNQM1, PINT, RAGQ, RQFAC, RQ1FAC, TSIGN, XPIN

!***FIRST EXECUTABLE STATEMENT  DCFODE
      GO TO (100, 200), METH

 100  ELCO(1,1) = 1.0D0
      ELCO(2,1) = 1.0D0
      TESCO(1,1) = 0.0D0
      TESCO(2,1) = 2.0D0
      TESCO(1,2) = 1.0D0
      TESCO(3,12) = 0.0D0
      PC(1) = 1.0D0
      RQFAC = 1.0D0
      DO 140 NQ = 2,12
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq-1).
! Initially, p(x) = 1.
!-----------------------------------------------------------------------
        RQ1FAC = RQFAC
        RQFAC = RQFAC/NQ
        NQM1 = NQ - 1
        FNQM1 = NQM1
        NQP1 = NQ + 1
! Form coefficients of p(x)*(x+nq-1). ----------------------------------
        PC(NQ) = 0.0D0
        DO 110 IB = 1,NQM1
          I = NQP1 - IB
 110      PC(I) = PC(I-1) + FNQM1*PC(I)
        PC(1) = FNQM1*PC(1)
! Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        PINT = PC(1)
        XPIN = PC(1)/2.0D0
        TSIGN = 1.0D0
        DO 120 I = 2,NQ
          TSIGN = -TSIGN
          PINT = PINT + TSIGN*PC(I)/I
 120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
! Store coefficients in ELCO and TESCO. --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC
        ELCO(2,NQ) = 1.0D0
        DO 130 I = 2,NQ
 130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
        AGAMQ = RQFAC*XPIN
        RAGQ = 1.0D0/AGAMQ
        TESCO(2,NQ) = RAGQ
        IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
        TESCO(3,NQM1) = RAGQ
 140    CONTINUE
      RETURN

 200  PC(1) = 1.0D0
      RQ1FAC = 1.0D0
      DO 230 NQ = 1,5
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq).
! Initially, p(x) = 1.
!-----------------------------------------------------------------------
        FNQ = NQ
        NQP1 = NQ + 1
! Form coefficients of p(x)*(x+nq). ------------------------------------
        PC(NQP1) = 0.0D0
        DO 210 IB = 1,NQ
          I = NQ + 2 - IB
 210      PC(I) = PC(I-1) + FNQ*PC(I)
        PC(1) = FNQ*PC(1)
! Store coefficients in ELCO and TESCO. --------------------------------
        DO 220 I = 1,NQP1
 220      ELCO(I,NQ) = PC(I)/PC(2)
        ELCO(2,NQ) = 1.0D0
        TESCO(1,NQ) = RQ1FAC
        TESCO(2,NQ) = NQP1/ELCO(1,NQ)
        TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
        RQ1FAC = RQ1FAC/FNQ
 230    CONTINUE
      RETURN
      END subroutine DCFODE
!----------------------- END OF SUBROUTINE DCFODE ----------------------


      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
!   BEGIN PROLOGUE  DEWSET
!   SUBSIDIARY
!   PURPOSE  Set error weight vector.
!   TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
!   AUTHOR  Hindmarsh, Alan C., (LLNL)
!   DESCRIPTION
!
!  This subroutine sets the error weight vector EWT according to
!      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
!  depending on the value of ITOL.
!
!   SEE ALSO  DLSODE
!   ROUTINES CALLED  (NONE)
!   REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDO! format.  (FNF)
!   890503  Minor cosmeti! changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   END PROLOGUE  DEWSET
!  End
      INTEGER :: N, ITOL
      INTEGER :: I
      DOUBLE PRECISION :: RTOL(*), ATOL(*), YCUR(N), EWT(N)

!   FIRST EXECUTABLE STATEMENT  DEWSET
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1,N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      RETURN
      END subroutine DEWSET

!----------------------- END OF SUBROUTINE DEWSET ----------------------


      DOUBLE PRECISION FUNCTION DUMACH()
!***BEGIN PROLOGUE  DUMACH
!***PURPOSE  Compute the unit roundoff of the machine.
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        DOUBLE PRECISION  A, DUMACH
!        A = DUMACH()
!
! *Function Return Values:
!     A : the unit roundoff of the machine.
!
! *Description:
!     The unit roundoff is defined as the smallest positive machine
!     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
!     in a machine-independent manner.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DUMSUM
!***REVISION HISTORY  (YYYYMMDD)
!   19930216  DATE WRITTEN
!   19930818  Added SLATEC-format prologue.  (FNF)
!   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
!***END PROLOGUE  DUMACH
!
      DOUBLE PRECISION :: U, COMP
!***FIRST EXECUTABLE STATEMENT  DUMACH
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      DUMACH = U*2.0D0
      RETURN
      END
!----------------------- End of Function DUMACH ------------------------

      SUBROUTINE DUMSUM(A,B,C)
!     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION :: A, B, C
      C = A + B
      RETURN
      END subroutine DUMSUM


      DOUBLE PRECISION FUNCTION DMNORM (N, V, W)
!-----------------------------------------------------------------------
! This function routine computes the weighted max-norm
! of the vector of length N contained in the array V, with weights
! contained in the array w of length N:
!   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
!-----------------------------------------------------------------------
      implicit none
      INTEGER, intent(in) :: N
      DOUBLE PRECISION, intent(in) :: V(N), W(N)
      integer :: I
      double precision :: VM
      VM = 0.0D0
      DO 10 I = 1,N
 10     VM = MAX(VM,ABS(V(I))*W(I))
      DMNORM = VM
      RETURN
      END
!----------------------- End of Function DMNORM ------------------------

      DOUBLE PRECISION FUNCTION DFNORM (N, A, W)
!-----------------------------------------------------------------------
! This function computes the norm of a full N by N matrix,
! stored in the array A, that is consistent with the weighted max-norm
! on vectors, with weights stored in the array W:
!   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
!-----------------------------------------------------------------------
      INTEGER :: N, I, J
      DOUBLE PRECISION :: A(N, N), W(N), AN, SUM
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        DO 10 J = 1,N
 10       SUM = SUM + ABS(A(I,J))/W(J)
          AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DFNORM = AN
      RETURN
!----------------------- End of Function DFNORM ------------------------
      END


      DOUBLE PRECISION FUNCTION DBNORM (N, A, NRA, ML, MU, W)
!-----------------------------------------------------------------------
! This function computes the norm of a banded N by N matrix,
! stored in the array A, that is consistent with the weighted max-norm
! on vectors, with weights stored in the array W.
! ML and MU are the lower and upper half-bandwidths of the matrix.
! NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
! In terms of the matrix elements a(i,j), the norm is given by:
!   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
!-----------------------------------------------------------------------
      INTEGER :: N, NRA, ML, MU
      INTEGER :: I, I1, JLO, JHI, J
      DOUBLE PRECISION :: A(NRA, N), W(N)
      DOUBLE PRECISION :: AN, SUM

      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        I1 = I + MU + 1
        JLO = MAX(I-ML,1)
        JHI = MIN(I+MU,N)
        DO 10 J = JLO,JHI
 10       SUM = SUM + ABS(A(I1-J,J))/W(J)
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DBNORM = AN
      RETURN
!----------------------- End of Function DBNORM ------------------------
      END


      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
!***BEGIN PROLOGUE  IXSAV
!***SUBSIDIARY
!***PURPOSE  Save and recall error message control parameters.
!***CATEGORY  R3C
!***TYPE      ALL (IXSAV-A)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  IXSAV saves and recalls one of two error message parameters:
!     LUNIT, the logical unit number to which messages are printed, and
!     MESFLG, the message print flag.
!  This is a modification of the SLATEC library routine J4SAVE.
!
!  Saved local variables..
!     LUNIT  = Logical unit number for messages.  The default is obtained
!             by a call to IUMACH (may be machine-dependent).
!     MESFLG = Print control flag..
!             1 means print all messages (the default).
!             0 means no printing.
!
!  On input..
!     IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
!     IVALUE = The value to be set for the parameter, if ISET = .TRUE.
!     ISET   = Logical flag to indicate whether to read or write.
!              If ISET = .TRUE., the parameter will be given
!              the value IVALUE.  If ISET = .FALSE., the parameter
!              will be unchanged, and IVALUE is a dummy argument.
!
!  On return..
!     IXSAV = The (old) value of the parameter.
!
!***SEE ALSO  XERRWD, XERRWV
!***ROUTINES CALLED  IUMACH
!***REVISION HISTORY  (YYMMDD)
!     921118  DATE WRITTEN
!     930329  Modified prologue to SLATEC format. (FNF)
!     930915  Added IUMACH call to get default output unit.  (ACH)
!     930922  Minor cosmetic changes. (FNF)
!     010425  Type declaration for IUMACH added. (ACH)
!***END PROLOGUE  IXSAV
!
! Subroutines called by IXSAV.. None
! Function routine called by IXSAV.. IUMACH
!-----------------------------------------------------------------------
!**End
      LOGICAL :: ISET
      INTEGER :: IPAR, IVALUE
!-----------------------------------------------------------------------
      INTEGER :: LUNIT, MESFLG !IUMACH
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/-1/, MESFLG/1/
!
!***FIRST EXECUTABLE STATEMENT  IXSAV
      IF (IPAR .EQ. 1) THEN
        IF (LUNIT .EQ. -1) LUNIT = IUMACH()
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
      ENDIF
!
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
      ENDIF
!
      RETURN
      end
!----------------------- End of Function IXSAV -------------------------


      INTEGER FUNCTION IUMACH()
!***BEGIN PROLOGUE  IUMACH
!***PURPOSE  Provide standard output unit number.
!***CATEGORY  R1
!***TYPE      INTEGER (IUMACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        INTEGER  LOUT, IUMACH
!        LOUT = IUMACH()
!
! *Function Return Values:
!     LOUT : the standard logical unit for Fortran output.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   930915  DATE WRITTEN
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  IUMACH
!
!*Internal Notes:
!  The built-in value of 6 is standard on a wide range of Fortran
!  systems.  This may be machine-dependent.
!**End
!***FIRST EXECUTABLE STATEMENT  IUMACH
      IUMACH = 6
!
      RETURN
!----------------------- End of Function IUMACH ------------------------
      END

!      *DECK IDAMAX
      INTEGER FUNCTION IDAMAX (N, DX, INCX)
!***BEGIN PROLOGUE  IDAMAX
!***PURPOSE  Find the smallest index of that component of a vector
!            having the maximum magnitude.
!***CATEGORY  D1A2
!***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!   IDAMAX  smallest index (zero if N .LE. 0)
!
!     Find smallest index of maximum magnitude of double precision DX.
!     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IDAMAX
      DOUBLE PRECISION :: DX(*), DMAX, XMAG
      INTEGER :: I, INCX, IX, N
!***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF (N .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN

      IF (INCX .EQ. 1) GOTO 20
!
!     Code for increments not equal to 1.
!
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = ABS(DX(IX))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
!
!     Code for increments equal to 1.
!
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
        XMAG = ABS(DX(I))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END

end module ode_lsoda
