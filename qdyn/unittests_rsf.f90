module unittests_rsf

  use constants, only : SOLVER_TYPE
  use problem_class, only : problem_type
  use mesh, only : mesh_get_size
  use solver, only : init_rk45
  use friction
  use fault_stress, only : init_kernel
  use unittests_aux
  use logger, only : log_msg

  implicit none
  private

  public :: test_rsf_friction

contains

!===============================================================================
! Instantiate variables for RSF testing
!===============================================================================
subroutine initiate_RSF(pb)
  type(problem_type) :: pb

  pb%a = 0.001
  pb%inv_a = 1 / pb%a
  pb%b = 0.0015
  pb%dc = 1e-5
  pb%v_star = 1e-6
  pb%mu_star = 0.6
  pb%v1 = 1e-3
  pb%v2 = 1e-1
  pb%inv_visc = 0d0

  pb%i_rns_law = 0
  pb%itheta_law = 1
  pb%slip = 0.d0
  pb%sigma = 2e7
  pb%tau = pb%mu_star * pb%sigma

  call set_theta_star(pb)
  pb%v = compute_velocity_RSF(pb%tau, pb%sigma, pb%theta, pb)
  pb%theta = pb%dc / pb%v

  ! Initiate solvers
  call initiate_solver(pb)

  call log_msg(" * RSF model set-up")

end subroutine initiate_RSF

!===============================================================================
! Subroutine to test code units using the rate-and-state friction framework
!===============================================================================
subroutine test_rsf_friction(pb)

  type(problem_type) :: pb
  double precision, dimension(pb%mesh%nn) :: dtheta
  double precision, dimension(pb%mesh%nn) :: tau0, tau1, v0, v1, v1_true
  double precision, dimension(pb%mesh%nn) :: dV_dtau, dV_dtheta, dV_dP
  double precision, dimension(pb%mesh%nn) :: dV_dtau2, dV_dtheta2, dV_dP2
  double precision, dimension(pb%mesh%nn) :: zero
  double precision :: atol, rtol, randno
  integer :: num_tests, num_passed, i
  logical :: pass, subpass1, subpass2, subpass3, subpass4
  character(255) :: msg

  num_tests = 0
  num_passed = 0
  atol = pb%acc
  rtol = pb%acc
  zero = 0.d0

  call log_msg("")
  call log_msg("Testing rate-and-state friction...")
  call log_msg("")

  call initiate_RSF(pb)

  ! - Initiate solver, etc. (initial values)
  ! - Solve ODE with a = 0 or b = 0
  ! - Compare with analytical solutions
  ! - Solve ODE with a, b > 0
  ! - Compare with interpolated data
  ! - Request output at fixed t, compare with data
  ! - Test both solvers
  ! - V-steps: test steady-state mu/theta

  ! Steady-state tests: ensure that dtheta/dt = 0 for theta = Dc/V
  ! Subtest ageing law
  pb%itheta_law = 1
  call dtheta_dt(pb%v, pb%theta, dtheta, pb)
  subpass1 = abs_assert_close(dtheta, zero, atol)

  ! Subtest slip law
  pb%itheta_law = 2
  call dtheta_dt(pb%v, pb%theta, dtheta, pb)
  subpass2 = abs_assert_close(dtheta, zero, atol)

  ! Collect results of subtests and print to screen
  pass = subpass1 .and. subpass2
  pb%test%test_passed = pb%test%test_passed .and. pass
  call print_result("Steady-state evolution laws", pass, num_passed, num_tests)
  call print_subresult("Steady-state ageing law", subpass1)
  call print_subresult("Steady-state slip law", subpass2)

  ! Steady-state friction tests: ensure that V1 = V0 * exp(dmu / (a-b))
  ! Initial velocity
  pb%tau = pb%mu_star * pb%sigma
  pb%theta = pb%dc / pb%v_star
  tau0 = pb%tau
  v0 = compute_velocity_RSF(tau0, pb%sigma, pb%theta, pb)
  ! Stress perturbation
  tau1 = 1.01 * tau0
  ! Anticipated final velocity
  v1_true = v0 * exp((tau1 - tau0) / (pb%sigma * (pb%a - pb%b)))
  ! Corresponding steady-state
  pb%theta = pb%dc / v1_true
  v1 = compute_velocity_RSF(tau1, pb%sigma, pb%theta, pb)
  pass = abs_assert_close(v1, v1_true, atol)
  pb%test%test_passed = pb%test%test_passed .and. pass
  call print_result("Steady-state friction (classical RSF)", pass, num_passed, num_tests)

  ! ! Steady-state friction tests for cut-off velocity model
  ! ! For V1 = V2: mu = mu* - (a-b)*log(V1/V + 1)
  ! pb%tau = pb%mu_star * pb%sigma
  ! tau0 = pb%tau
  ! pb%theta = pb%dc / pb%v_star
  ! v0 = compute_velocity_RSF(tau0, pb%sigma, pb%theta, pb)
  !
  ! pb%v1 = pb%v2
  ! pb%i_rns_law = 1
  ! call set_theta_star(pb)
  ! mu = friction_mu(pb%v, pb%theta, pb)
  ! mu_truth = pb%mu_star - (pb%a - pb%b)*log(pb%v1/pb%v + 1d0)
  ! subpass1 = abs_assert_close(mu, mu_truth, atol)

  ! ! For V1 = V2 << V: mu -> mu*
  ! pb%v = 1e5*pb%v1
  ! pb%theta = pb%dc / pb%v
  ! mu = friction_mu(pb%v, pb%theta, pb)
  ! mu_truth = pb%mu_star
  ! subpass2 = abs_assert_close(mu, mu_truth, atol)
  !
  ! ! For V = V2 >> V1: mu -> mu* + b*log(2)
  ! pb%v1 = 1e-9
  ! pb%v2 = 1e0
  ! pb%v = pb%v2
  ! pb%theta = pb%dc / pb%v
  ! call set_theta_star(pb)
  ! mu = friction_mu(pb%v, pb%theta, pb)
  ! mu_truth = pb%mu_star + pb%b*log(2d0)
  ! subpass3 = abs_assert_close(mu, mu_truth, atol)
  !
  ! ! Collect results of subtests and print to screen
  ! pass = subpass1 .and. subpass2 .and. subpass3
  ! pb%test%test_passed = pb%test%test_passed .and. pass
  ! call print_result("Steady-state friction (cut-off velocities)", pass, num_passed, num_tests)
  ! call print_subresult("V1 = V2", subpass1)
  ! call print_subresult("V1 = V2 << V", subpass2)
  ! call print_subresult("V1 << V2 = V", subpass3)

  ! Ensure that the classical rate-and-state friction law (i_rns_law = 0)
  ! is identical to the regularised RSF law (i_rns_law = 2), for various
  ! V and theta (random) -- NOTE that the random seed is set to default,
  ! so that randno is deterministic between consecutive test runs
  pass = .true.
  subpass1 = .true.
  subpass2 = .true.
  subpass3 = .true.
  subpass4 = .true.

  pb%tau = 0.99 * pb%mu_star * pb%sigma
  pb%theta = pb%dc / pb%v_star

  do i = 1, 1000
    call random_number(randno)
    randno = randno + 0.5

    pb%i_rns_law = 0
    call set_theta_star(pb)
    v0 = compute_velocity_RSF(pb%tau, pb%sigma, pb%theta*randno, pb)
    call RSF_derivs(dV_dtau, dV_dtheta, dV_dP, v0, pb%theta*randno, pb%tau, &
                    pb%sigma, pb)
    pb%i_rns_law = 2
    call set_theta_star(pb)
    v1 = compute_velocity_RSF(pb%tau, pb%sigma, pb%theta*randno, pb)
    call RSF_derivs(dV_dtau2, dV_dtheta2, dV_dP2, v1, pb%theta*randno, pb%tau, &
                    pb%sigma, pb)
    ! Test if mu, dmu/dv, and dmu/dtheta are ok for all random values of theta
    ! If one of the assertions fails, subpass will become .false.
    subpass1 = subpass1 .and. abs_assert_close(v0, v1, atol)
    subpass2 = subpass2 .and. rel_assert_close(dV_dtau, dV_dtau2, rtol)
    subpass3 = subpass3 .and. rel_assert_close(dV_dtheta, dV_dtheta2, rtol)
    subpass4 = subpass4 .and. rel_assert_close(dV_dP, dV_dP2, rtol)

  enddo

  pass = subpass1 .and. subpass2 .and. subpass3 .and. subpass4
  pb%test%test_passed = pb%test%test_passed .and. pass
  call print_result("Compare classical / regularised RSF", pass, num_passed, num_tests)
  call print_subresult("Slip rate", subpass1)
  call print_subresult("dV/dtau", subpass2)
  call print_subresult("dV/dtheta", subpass3)
  call print_subresult("dV/dP", subpass4)


  write(msg, "(A, I0, A, I0, A)") " Rate-and-state friction: ", num_passed, " / ", num_tests, " passed"
  call log_msg(msg)

end subroutine test_rsf_friction

end module unittests_rsf
