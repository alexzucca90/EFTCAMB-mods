!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 09p3_Designer_GBD.f90
!! This file contains the relevant code for designer GBD models.


!----------------------------------------------------------------------------------------
!> This module contains the relevant code for designer GBD models.

!> @author Bin Hu, Marco Raveri, Simone Peirone

!> @author Alex Zucca

module EFTCAMB_designer_GBD

    use precision
    use IniFile
    use AMLutils
    use equispaced_linear_interpolation_1D
    use EFT_def
    use EFTCAMB_rootfind
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_neutral_parametrization_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_CPL_parametrizations_1D
    use EFTCAMB_JBP_parametrizations_1D
    use EFTCAMB_turning_point_parametrizations_1D
    use EFTCAMB_taylor_parametrizations_1D
    use EFTCAMB_abstract_model_designer

!> adding the interpolated function
use EFTCAMB_interpolated_function_1D

    implicit none

    private

    public EFTCAMB_GBD_designer

    !----------------------------------------------------------------------------------------
    !> This is the designer GBD model. Inherits from the abstract designer model and has the
    !! freedom of defining the expansion history.
    type, extends ( EFTCAMB_designer_model ) :: EFTCAMB_GBD_designer

        ! theory parameters:
        real(dl) :: phi_ini, dphi_ini                                     !< The initial values of the scalar dof and its derivative
        real(dl) :: xi                                                    !< The value of the coupling constant xi


        ! the pure EFT functions model selection flags:
        integer  :: EFTwDE                                                !< Model selection flag for designer GBD w DE.

        !> GBD coupling type
        integer :: coupling_type                                          !< which coupling is being used for the reconstruction

        ! the pure EFT functions:
        class( parametrized_function_1D ), allocatable ::  DesGBDwDE     !< The pure EFT function w_DE.

        ! the interpolated EFT functions that come out of the background sover:
        type(equispaced_linear_interpolate_function_1D) :: EFTOmega       !< The interpolated function Omega (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda      !< The interpolated function Lambda (and derivatives).
        type(equispaced_linear_interpolate_function_1D) :: EFTc           !< The interpolated funtcion c (and derivatives).

        ! some designer parameters:
        integer  :: designer_num_points = 1000                            !< Number of points sampled by the designer code.
        real(dl) :: x_initial           = log(10._dl**(-8._dl))           !< log(a start)
        real(dl) :: x_final             = 0.0_dl                          !< log(a final)

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => EFTCAMBDesignerGBDReadModelSelectionFromFile   !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => EFTCAMBDesignerGBDAllocateModelSelection       !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => EFTCAMBDesignerGBDInitModelParameters          !< subroutine taht initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => EFTCAMBDesignerGBDInitModelParametersFromFile  !< subroutine that reads the parameters of the model from file.

        ! background solver:
        procedure :: initialize_background           => EFTCAMBDesignerGBDInitBackground               !< subroutine that initializes the background of designer GBD.
        procedure :: solve_designer_equations        => EFTCAMBDesignerGBDSolveDesignerEquations       !< subroutine that solves the designer GBD background equations.
        !procedure :: find_initial_conditions         => EFTCAMBDesignerGBDFindInitialConditions        !< subroutine that solves the background equations several time to determine the values of the initial conditions.
        procedure :: omega_phi                       => EFTCAMBDesignerGBDCoupling                      !< function that computes omega(phi) and its derivatives

        ! utility functions:
        procedure :: compute_param_number  => EFTCAMBDesignerFRComputeParametersNumber     !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => EFTCAMBDesignerFRFeedback                    !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => EFTCAMBDesignerFRParameterNames              !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => EFTCAMBDesignerFRParameterNamesLatex         !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => EFTCAMBDesignerFRParameterValues             !< subroutine that returns the i-th parameter value.

        ! CAMB related procedures:
        procedure :: compute_background_EFT_functions  => EFTCAMBDesignerFRBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure :: compute_secondorder_EFT_functions => EFTCAMBDesignerFRSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure :: compute_dtauda                    => EFTCAMBDesignerFRComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure :: compute_adotoa                    => EFTCAMBDesignerFRComputeAdotoa            !< subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        procedure :: compute_H_derivs                  => EFTCAMBDesignerFRComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

        ! stability procedures:
        procedure :: additional_model_stability        => EFTCAMBDesignerFRAdditionalModelStability !< function that computes model specific stability requirements.

    end type EFTCAMB_GBD_designer

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBDesignerGBDReadModelSelectionFromFile( self, Ini )

        implicit none

        class(EFTCAMB_GBD_designer)  :: self   !< the base class
        type(TIniFile)               :: Ini     !< Input ini file

        ! read model selection flags:
        self%EFTwDE             = Ini_Read_Int_File( Ini, 'EFTwDE', 0 )

        !> read coupling type
        self%coupling_type      = Ini_Read_Int_File( Ini, 'GBD_coupling_type', 3)

    end subroutine EFTCAMBDesignerFRReadModelSelectionFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine EFTCAMBDesignerGBDAllocateModelSelection( self )

        implicit none

        class(EFTCAMB_GBD_designer)                       :: self              !< the base class
        character, allocatable, dimension(:)              :: param_names       !< an array of strings containing the names of the function parameters
        character, allocatable, dimension(:)              :: param_names_latex !< an array of strings containing the latex names of the function parameters

        ! allocate wDE:
        if ( allocated(self%DesGBDwDE) ) deallocate(self%DesGBDwDE)
        select case ( self%EFTwDE )
            case(0)
                allocate( wDE_LCDM_parametrization_1D::self%DesGBDwDE )
            case(1)
                allocate( constant_parametrization_1D::self%DesGBDwDE )
            case(2)
                allocate( CPL_parametrization_1D::self%DesGBDwDE )
                call self%DesGBDwDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
            case(3)
                allocate( JBP_parametrization_1D::self%DesGBDwDE )
                call self%DesGBDwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
            case(4)
                allocate( turning_point_parametrization_1D::self%DesGBDwDE )
                call self%DesGBDwDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
            case(5)
                allocate( taylor_parametrization_1D::self%DesGBDwDE )
                call self%DesGBDwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
            case(6)
                allocate( interpolated_function_1D::self%DesGBDwDE )
                call self%DesGBDwDE%set_param_names(['wDE_filename  '])
            case default
                write(*,'(a,I3)') 'No model corresponding to EFTwDE =', self%EFTwDE
                write(*,'(a)')    'Please select an appropriate model.'
        end select

        ! initialize the names:
        call self%DesGBDwDE%set_name( 'EFTw', 'w' )

    end subroutine EFTCAMBDesignerGBDAllocateModelSelection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the model parameters based on the values found in an input array.
    subroutine EFTCAMBDesignerGBDInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_GBD_designer)                            :: self   !< the base class
        real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.
        real(dl), dimension(self%parameter_number -1)          :: temp
        integer                                                :: i

        self%B0 = array(1)

        do i = 1, self%parameter_number -1
            temp(i) = array(i+1)
        end do
        call self%DesGBDwDE%init_parameters(temp)

    end subroutine EFTCAMBDesignerGBDInitModelParameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine EFTCAMBDesignerGBDInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_GBD_designer)  :: self   !< the base class
        type(TIniFile)              :: Ini    !< Input ini file

        !> read phi_ini
        self%phi_ini = Ini_Read_Double_File( Ini, 'GBD_phi_ini', 0._dl)

        !> read dphi_ini
        self%phi_ini = Ini_Read_Double_File( Ini, 'GBD_phi_ini', 0._dl)

        !> read xi
        self%xi      = Ini_Read_Double_File( Ini, 'GBD_xi',      1._dl)


        !> read w_DE parameters:
        call self%DesGBDwDE%init_from_file( Ini )

    end subroutine EFTCAMBDesignerGBDInitModelParametersFromFile

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes the coupling function and its derivatives given the value of the field
    function EFTCAMBDesignerGBDCoupling( self, phi, deriv) result(Omega)

        implicit none

        class(EFTCAMB_GBD_designer) :: self     !< the base class

        real(dl) :: phi         !< the value of the field
        integer  :: deriv       !< the derivative of the coupling function that needs to be computed

        real(dl) :: Omega       !< the result, Omega(phi) or its derivatives w.r.t. phi

        if(self%coupling_type == 1) then

            !> Linear coupling

            if (deriv == 0) then
                Omega = 1._dl + self%xi*phi

            else if (deriv == 1) then
                Omega = self%xi

            else if (deriv == 2) then
                Omega = 0._dl

            else if (deriv == 3) then
                Omega = 0._dl

            else
                write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong derivative"
                stop
            end if


        else if (self%coupling_type == 2) then

            !> Quadratic coupling

            if (deriv == 0) then
                Omega = 1._dl + self%xi*phi**2

            else if (deriv == 1) then
                Omega = 2._dl*self%xi*phi

            else if (deriv == 2) then
                Omega = 2.+dl*self%xi

            else if (deriv == 3) then
                Omega = 0._dl

            else
                write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong derivative"
                stop
            end if


        else if (self%coupling_type == 3) then

            !> Exponential coupling

            if (deriv == 0) then
                Omega = exp(self%xi*phi)

            else if (deriv == 1) then
                Omega = self%xi*exp(self%xi*phi)

            else if (deriv == 2) then
                Omega = self%xi**2 * exp(self%xi*phi)

            else if (deriv == 3) then
                Omega = self%xi**3 * exp(self%xi*phi)

            else
                write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong derivative"
                stop
            end if

        else if (self%coupling_type == 4) then

            !> Negative exponential coupling
            if (deriv == 0) then
                Omega = exp(-self%xi*phi)

            else if (deriv == 1) then
                Omega = - self%xi * exp(-self%xi*phi)

            else if (deriv == 2) then
                Omega = self%xi**2 * exp(-self%xi*phi)

            else if (deriv == 3) then
                Omega = -self%xi**3 * exp(-self%xi*phi)

            else
                write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong derivative"
                stop
            end if

        else
            write(*,*) "ERROR in EFTCAMBDesignerGBDCoupling: wrong coupling type"
            stop
        end if

    end function EFTCAMBDesignerGBDCoupling

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of designer GBD
    subroutine EFTCAMBDesignerGBDInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(EFTCAMB_GBD_designer)                     :: self             !< the base class
        type(EFTCAMB_parameter_cache),  intent(in)      :: params_cache     !< a EFTCAMB parameter cache containing cosmological parameters
        integer,                        intent(in)      :: feedback_level   !< a level of feedback from the background code. 0=none; 1=some; 2=chatty
        logical,                        intent(out)     :: success          !< whether the background initialization succeded or not

        real(dl) :: phi_i, dphi_i


        !> some feedback
        if (feedback_level>0) then
            write(*,'(a)') "***************************************************************"
            write(*,'(a)') "EFTCAMB designer GBD background solver"
            write(*,'(a)')
        end if

        !> initialize interpolating functions
        call self%EFTOmega%initialize   ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTLambda%initialize  ( self%designer_num_points, self%x_initial, self%x_final )
        call self%EFTc%initialize       ( self%designer_num_points, self%x_initial, self%x_final )


        !> solve the background equations and store the solution:
        call self%solve_designer_equations( params_cache, success=success )

    end subroutine EFTCAMBDesignerGBDInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the designer GBD background equations.
    subroutine EFTCAMBDesignerGBDSolveDesignerEquations( self, params_cache, success )

        implicit none

        class(EFTCAMB_GBD_designer)                 :: self          !< the base class.
        type(EFTCAMB_parameter_cache), intent(in)   :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
        logical , intent(out)                       :: success       !< whether the calculation ended correctly or not

        integer, parameter :: num_eq = 2   !<  Number of equations

        real(dl) :: Omegam_EFT, Omegavac_EFT, OmegaMassiveNu_EFT, OmegaGamma_EFT, OmegaNu_EFT
        real(dl) :: Omegarad_EFT, EquivalenceScale_fR, Ratio_fR, Initial_B_fR, Initial_C_fR
        real(dl) :: PPlus, yPlus, CoeffA_Part, yStar, x

        real(dl) :: y(num_eq+1), ydot(num_eq+1)

        integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
        real(dl) :: rtol, atol, t1, t2, B
        real(dl), allocatable :: rwork(:)
        integer,  allocatable :: iwork(:)

        ! digedt the input:
        !if ( .not. present(only_B0) ) only_B0 = .False.

        ! 1) Cosmological parameters:
        Omegam_EFT         = params_cache%omegab + params_cache%omegac
        Omegavac_EFT       = params_cache%omegav
        OmegaMassiveNu_EFT = params_cache%omegan
        OmegaGamma_EFT     = params_cache%omegag
        OmegaNu_EFT        = params_cache%omegar

        Omegarad_EFT       = OmegaGamma_EFT + OmegaNu_EFT

        EquivalenceScale_GBD = (Omegarad_EFT+ OmegaMassiveNu_EFT)/Omegam_EFT
        Ratio_GBD = EquivalenceScale_GBD/Exp( self%x_initial )

        ! 2) Set initial conditions:
        !> need to modify this
        x    = self%x_initial
        y(1) = self%phi_ini
        y(2) = self%dphi_ini
        ydot = 0._dl

        !> fill the interpolated EFT functions here

        ! 3) Solve the equation of motion
        !> defining the time-step


        !> Loop to fill the interpolation arrays
        do  i = 1, self%EFTOmega%num_points-1

            !> calling the solver, in this case gl10
            call gl10(num_eq+1,derivs, y, dN)

            !> check if the solution is acceptable
            if (.not.(y .ge. 0._dl .or. y .le. 0._dl)) then
                success = .False.
                return
            end if

            !> filling the interpolation arrays
            call output(num_eq+1,y,i)

        end do

        !> end of the reconstruction
        return

    contains

        subroutine derivs(num,y,yprime)

            implicit none

            integer :: num                          !< number of equations
            real(dl), dimension(num) :: y, yprime   !< input status of the system and output derivative
            real(dl) :: N,a,a2                      !< efolds number, scale factor and squared scale factor

            real(dl) :: om, omp, ompp               !< coupling functoin and its derivatives
            real(dl) :: H, adotdot                  !< expansion history parameters

            real(dl) :: rhonu, presnu               !< massive neutrinos background variables
            real(dl) :: rhonu_tot, presnu_tot
            real(dl) :: grhormass_t
            integer  :: nu_i

            real(dl) :: Em,   Er,   Enu,   X        !< normalized energy densities
            real(dl) :: Em_p, Er_p, Enu_p, X_p      !< derivatives of normalized energy densities

            real(dl) :: EFT_E_gfun                  !< effective dark energy variables
            real(dl) :: EFT_E_gfunp
            real(dl) :: EFT_E_gfunpp
            real(dl) :: EFT_E_gfunppp


            associate (N => y(1),  phi => y(2), pi=>y(3), dN => yprime(1), dphi => yprime(2) , ddphi => yprime(3))

                !> convert N in a
                a = exp(N)

                !> compute coupling function and its derivatives
                om      = self%omega_phi(phi, 0)
                omp     = self%omega_phi(phi, 1)
                ompp    = self%omega_phi(phi, 2)

                !> normalized energy densities
                Em = Omegam_EFT   * exp(-3._dl*N)
                Er = OmegaRad_EFT * exp(-4._dl*N)

                !> compute the function g(x) and its derivatives:
                EFT_E_gfun    = -(Log( self%DesGBDwDE%integral(a) ) -2._dl*x)/3._dl
                EFT_E_gfunp   = 1._dl +self%DesGBDwDE%value(a)
                EFT_E_gfunpp  = Exp(x)*self%DesGBDwDE%first_derivative(a)
                EFT_E_gfunppp = Exp(x)*self%DesGBDwDE%first_derivative(a) +Exp(2._dl*x)*self%DesGBDwDE%second_derivative(a)

                !> compute the normalized dark energy density
                X   = Omegavac_EFT*exp(-3._dl*EFT_E_gfun)
                X_p = -3._dl * Omegavac_EFT* exp(-3._dl*EFT_E_gfun)* EFT_E_gfunp

                !> compute massive neutrinos contribution:
                rhonu_tot   = 0._dl
                presnu_tot  = 0._dl
                Enu         = 0._dl
                Enu_p       = 0._dl

                if ( params_cache%Num_Nu_Massive /= 0) then
                    do nu_i = 1, params_cache%Nu_mass_eigenstates

                        rhonu  = 0._dl
                        presnu = 0._dl
                        grhormass_t= params_cache%grhormass(nu_i)/a**2

                        call params_cache%Nu_background(a*params_cache%nu_masses(nu_i),rhonu,presnu)

                        rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                        presnu_tot = presnu_tot + grhormass_t*presnu

                        Enu    = Enu   + params_cache%grhormass(nu_i)/3._dl/a**4/params_cache%h0_Mpc**2*rhonu
                        Enu_p  = Enu_p - params_cache%grhormass(nu_i)/params_cache%h0_Mpc**2/a**4*(rhonu +presnu)

                    end do
                end if

                !> now the equations of motion

                !> derivative of N=ln(a)
                dN = 1.d0

                !> phi prime
                dphi = pi

                !> phi prime prime
                ddphi = ((om - 1._dl)*(3._dl * Em + 4._dl * Er - Enu_p) - X_p)/omp/(Em+Er+Enu+X) &
                        - (1._dl + ompp)*dphi**2/omp + (1._dl+0.5_dl*(3._dl*Em+4._dl*Er-Enu_p-X_p)/(Em+Er+Enu+X))*dphi

            end associate

        end subroutine derivs

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that takes the solution of the background GBD equations and stores the values of
        !> the EFT functions.
        subroutine output( num, y, ind)

            implicit none

            integer , intent(in)                    :: num  !< number of equations in the ODE system.
            integer , intent(in)                    :: ind  !< index of the EFT functions interpolation tables to fill.
            real(dl), intent(in) , dimension(num)   :: y    !< input status of the system.

            real(dl) :: ydot(num)                           !< array of derivatives of the system
            real(dl) :: EFT_E_gfun, EFT_E_gfunp, EFT_E_gfunpp, EFT_E_gfunppp !< effective dark energy variables

            real(dl) :: om, omp ompp, omppp         !< coupling function and its derivatives

            !> some massive neutrinos variables
            real(dl) :: rhonu_tot, presnu_tot, presnudot_tot, presnudotdot_tot
            real(dl) :: rhonu, presnu, grhormass_t

            real(dl) :: Em,    Er,    Enu,    X     !< normalized energy densities
            real(dl) :: Em_p,  Er_p,  Enu_p,  X_p   !< derivatives of normalized energy densities
            real(dl) :: Em_pp, Er_pp, Enu_pp, X_pp  !< second derivatives of normalized energy densities
            real(dl) :: Enu_ppp                     !< third derivative of normalized energy density

            integer  :: nu_i
            logical  :: is_open

            !> other parameters
            real(dl) :: a, N, phi, dphi, ddphi
            real(dl) :: adotoa, Hdot, adotdotoa
            real(dl) :: calF, V, Vprime, Vdot
            real(dl) :: Etot

            !> convert N in a
            N = y(1)
            a = exp(N)

            !> extract values of the field
            phi     = y(2)
            dphi    = y(3)

            !> compute coupling function and its derivatives
            om      = self%omega_phi(phi, 0)
            omp     = self%omega_phi(phi, 1)
            ompp    = self%omega_phi(phi, 2)
            omppp   = self%omega_phi(phi, 3)

            !> normalized energy densities
            Em = Omegam_EFT   * exp(-3._dl*N)
            Er = OmegaRad_EFT * exp(-4._dl*N)

            !> compute the function g(x) and its derivatives:
            EFT_E_gfun    = -(Log( self%DesfRwDE%integral(a) ) -2._dl*x)/3._dl
            EFT_E_gfunp   = 1._dl +self%DesfRwDE%value(a)
            EFT_E_gfunpp  = Exp(x)*self%DesfRwDE%first_derivative(a)
            EFT_E_gfunppp = Exp(x)*self%DesfRwDE%first_derivative(a) +Exp(2._dl*x)*self%DesfRwDE%second_derivative(a)

            !> compute the normalized dark energy density
            X       = Omegavac_EFT*exp(-3._dl*EFT_E_gfun)
            X_p     = -3._dl * Omegavac_EFT * exp(-3._dl*EFT_E_gfun)* EFT_E_gfunp
            X_pp    =  3._dl * Omegavac_EFT * exp(-3._dl*EFT_E_gfun)*(3._dl* EFT_E_gfunp**2 -  EFT_E_gfunpp)

            !> compute massive neutrinos contribution:
            rhonu_tot   = 0._dl
            presnu_tot  = 0._dl
            Enu         = 0._dl
            Enu_p       = 0._dl

            if ( params_cache%Num_Nu_Massive /= 0) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates

                    rhonu  = 0._dl
                    presnu = 0._dl
                    grhormass_t= params_cache%grhormass(nu_i)/a**2
                    call params_cache%Nu_background(a*params_cache%nu_masses(nu_i),rhonu,presnu)
                    rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                    presnu_tot = presnu_tot + grhormass_t*presnu

                    Enu     = Enu   + params_cache%grhormass(nu_i)/3._dl/a**4/params_cache%h0_Mpc**2*rhonu
                    Enu_p   = Enu_p - params_cache%grhormass(nu_i)/params_cache%h0_Mpc**2/a**4*(rhonu +presnu)

                end do
            end if

            !> calculate H, Hdot and adotdotoa
            adotoa      = a * params_cache%h0_Mpc* sqrt(Em+Er+Enu+X)
            Hdot        = a**2 * params_cache%h0_Mpc**2 * ( (Em+Er+Enu+X) + &
                          0.5_dl * (-3._dl * Em - 4._dl*Er + Enu_p + X_p) )
            adotdotoa   = Hdot + adotoa**2

            !> start filling the EFT interpolated functions
            call derivs( num, y, ydot )

            !> extract ddphi
            ddphi = ydot(3)

            !> calculating the potential - needed for the function Lambda
            !> begin with {\cal F} in the notes
            calF = om - dphi**2 / 6._dl + omp * dphi

            !> inverted Friedmann equation
            V = 3._dl*(Em+Er+Enu+X)*calF - 3.d0*(Em+Er+Enu)
            V = V * params_cache%h0_Mpc**2 ! this is in Mpc-2

            !> get Vprime from the equation of motion for the scalar field
            Vprime = (3._dl * omp * adotdotoa - adotoa**2 *ddphi - (adotoa**2 + adotdotoa)*dphi)/a**2 ! in Mpc-2

            !> compute Vdot - needed for Lambda dot -
            Vdot = Vprime * dphi * adotoa

            !------> calculating dddphi: this requires many steps
            ! 1) Compute everything of massive nu again to get the time derivatives:
            rhonu_tot        = 0._dl
            presnu_tot       = 0._dl
            presnudot_tot    = 0._dl
            presnudotdot_tot = 0._dl

            Enu     = 0._dl
            Enu_p   = 0._dl
            Enu_pp  = 0._dl
            Enu_ppp = 0._dl

            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates

                    rhonu        = 0._dl
                    presnu       = 0._dl
                    presnudot    = 0._dl
                    presnudotdot = 0._dl

                    grhormass_t  = params_cache%grhormass(nu_i)/a**2

                    call params_cache%Nu_background(a*params_cache%nu_masses(nu_i),rhonu,presnu)
                    presnudot = params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),adotoa,presnu)
                    presnudotdot = params_cache%Nu_pidotdot(a*params_cache%nu_masses(nu_i),adotoa,Hdot,presnu,presnudot)

                    rhonu_tot  = rhonu_tot + grhormass_t*rhonu
                    presnu_tot = presnu_tot + grhormass_t*presnu
                    presnudot_tot  = presnudot_tot + grhormass_t*(presnudot -4._dl*adotoa*presnu)
                    presnudotdot_tot = presnudotdot_tot + grhormass_t*(presnudotdot &
                                    & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))

                    Enu     = Enu   + params_cache%grhormass(nu_i)/3._dl/a**4/params_cache%h0_Mpc**2*rhonu
                    Enu_p   = Enu_p - params_cache%grhormass(nu_i)/params_cache%h0_Mpc**2/a**4*(rhonu +presnu)
                    Enu_pp  = Enu_pp + 3._dl/params_cache%h0_Mpc**2*params_cache%grhormass(nu_i)/a**4*(rhonu +presnu)&
                                & -grhormass_t*(presnudot -4._dl*adotoa*presnu)/params_cache%h0_Mpc**3/sqrt(EFunction)/a**3
                    Enu_ppp = Enu_ppp -9._dl/params_cache%h0_Mpc**2*params_cache%grhormass(nu_i)/a**4*(rhonu +presnu)&
                                & +(3._dl/adotoa/params_cache%h0_Mpc**2/a**2+Hdot/adotoa**3/params_cache%h0_Mpc**2/a**2)&
                                &*grhormass_t*(presnudot -4._dl*adotoa*presnu)&
                                & -grhormass_t*(presnudotdot &
                                & -8._dl*adotoa*presnudot +4._dl*presnu*(+4._dl*adotoa**2-Hdot))/adotoa**2/params_cache%h0_Mpc**2/a**2
                end do
            end if

            !> E tot
            Etot = Em + Er + Enu + X

            !> now compute dddphi
            dddphi = (3._dl*Em + 4._dl * Er - Enu_p)/Etot*dphi                                              &
                    -(om - 1._dl)*(3._dl*Em + 4._dl*Er - Enu_p)/Etot*ompp/omp**2 *dphi                      &
                    +(om-1._dl)/Etot/omp * (- 9._dl*Em -16._dl*Er - Enu_pp)                                 &
                    -(om-1._dl)/omp/Etot**2*(3._dl*Em+4._dl*Er - Enu_p)*(-3._dl*Em-4._dl*Er+Enu_p + X_p)    &
                    + X_p * ompp * dphi / omp**2 / Etot                                                     &
                    - X_pp / omp / Etot                                                                     &
                    + X_p * (-3._dl*Em-4._dl*Er + Enu_p + X_p)/omp/Etot**2                                  &
                    +(1._dl+ompp)/omp**2 * dphi**3 * ompp                                                   &
                    - omppp*dphi**3/omp                                                                     &
                    -2._dl * (1._dl+ompp)*dphi*ddphi/omp                                                    &
                    +(0.5_dl * (-9._dl*Em - 16._dl*Er - Enu_pp - X_pp)/Etot                                 &
                    -0.5_dl*(3._dl*Em+4._dl*Er-Enu_p-X_p)*(-3._dl*Em-4._dl*Er+Enu_p+X_p)/Etot**2)*dphi      &
                    +0.5_dl* (5._dl*Em+6._dl*Er +2._dl*Enu+2._dl*X - Enu_p-X_p)/Etot*ddphi


            !> Filling the EFT functions:
            !> Omega
            self%EFTOmega%y(ind)    = om
            self%EFTOmega%yp(ind)   = omp * dphi / a
            self%EFTOmega%ypp(ind)  = (ompp * dphi**2 + omp*ddphi - omp * dphi)/a**2
            self%EFTOmega%yppp(ind) = (omppp * dphi**3 + 3._dl * ompp * dphi * ddphi - 3._dl* ompp*dphi**2 &
                                        - 3._dl*omp * ddphi + omp * dddphi + 2._dl*omp *dphi)/a**3

            !> c
            self%EFTc%y(ind)    = 0.5_dl*adotoa**2 * dphi**2
            self%EFTc%yp(ind)   = (- dphi**2 + dphi * ddphi)*adotoa**3 + dphi**2 * adotoa * Hdot

            !> Lambda
            self%EFTLambda%y(ind)   = self%EFTc%y(ind)  - V    * a**2
            self%EFTLambda%yp(ind)  = self%EFTc%yp(ind) - Vdot * a**2


            !> this has to be added later
            !if ( DebugEFTCAMB ) then
            !    inquire( unit=33, opened=is_open )
            !    if ( is_open ) then
            !        write (33,'(20E15.5)') x, Exp(x), Ricci, y(1), &
            !        & self%EFTOmega%y(ind), self%EFTOmega%yp(ind), self%EFTOmega%ypp(ind), self%EFTOmega%yppp(ind), self%EFTLambda%y(ind), self%EFTLambda%yp(ind),&
            !        & B
            !    end if
            !end if


        end subroutine output

    end subroutine EFTCAMBDesignerGBDSolveDesignerEquations


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that solves the background equations several time to determine the values of
    !! the initial conditions. Maps A at initial times to B0 today.
    subroutine EFTCAMBDesignerFRFindInitialConditions( self, params_cache, feedback_level, A_ini, success )

        implicit none

        class(EFTCAMB_fR_designer)                   :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        real(dl)                     , intent(out)   :: A_ini          !< value of A_initial that gives self%B0 today
        logical                      , intent(out)   :: success        !< whether the calculation ended correctly or not

        real(dl) :: ATemp1, ATemp2, BTemp1, BTemp2, HorizAsyntB
        real(dl) :: VertAsyntA, ATemp3, ATemp4, BTemp3, BTemp4, realAp
        integer  :: ind

        !  Find initial conditions for designer F(R).
        !  Initial conditions are found solving B0(A) = B0wanted.
        !  The next part of code is complicated by the fact that B0(A) is not continuous and we have to adopt ad-hoc strategies.

        ! 1) Find the horizontal asymptote of B0(A)
        ATemp1 = 0._dl
        ATemp2 = 0._dl
        do ind=10, 100, 1
            ATemp1 = ATemp2
            ATemp2 = 10._dl**REAL(ind)
            BTemp1 = DesFR_BfuncA(ATemp1)
            BTemp2 = DesFR_BfuncA(ATemp2)
            if (Abs(BTemp1-BTemp2)/Abs(ATemp1-ATemp2).lt.1.d-100) exit
        end do
        HorizAsyntB = BTemp2

        if ( feedback_level>1 ) write(*,'(a,E13.4)') '   horizontal asymptote = ', HorizAsyntB

        ! 2) Check that the value of B0 given is not the forbidden one: the one corresponding to the horizontal asynt.
        if ( ABS(HorizAsyntB-self%B0)<1.d-15 ) then
            success=.false.
            return
        end if

        ! 3) Bracket the vertical asyntote
        ATemp1 = -10._dl
        ATemp2 = 10._dl
        call zbrac(DesFR_BfuncA,ATemp1,ATemp2,success,HorizAsyntB)
        if (.not.success) then
            if ( feedback_level>1 ) write(*,'(a)') '   FAILURE of vertical asymptote bracketing'
            return
        end if

        ! 4) Find the vertical asyntote by tricking the root finding algorithm.
        VertAsyntA = zbrent(DesFR_BfuncA,ATemp1,ATemp2,1.d-100,HorizAsyntB,success)
        if (.not.success) then
            if ( feedback_level>1 ) write(*,'(a)') '   FAILURE of vertical asyntote finding'
            return
        end if

        if ( feedback_level>1 ) write(*,'(a,E13.4)') '   vertical asymptote   = ', VertAsyntA

        ! 5) Find values for A that are on the left and on the right of the asyntote.
        do ind=-10, -1, 1
            ATemp1 = VertAsyntA+10._dl**ind
            ATemp2 = VertAsyntA-10._dl**ind
            BTemp1 = DesFR_BfuncA(ATemp1)
            BTemp2 = DesFR_BfuncA(ATemp2)
            if (BTemp1*BTemp2<0._dl) exit
        end do

        ! 6) Extablish on which side of the asyntote to look for the solution.
        ATemp3 = ATemp1
        ATemp4 = ATemp2
        do ind=1, 10, 1
            ATemp3 = ATemp3 + 10._dl**ind
            ATemp4 = ATemp4 - 10._dl**ind
            BTemp3 = DesFR_BfuncA(ATemp3)
            BTemp4 = DesFR_BfuncA(ATemp4)
            if ((BTemp1-self%B0)*(BTemp3-self%B0)<0.or.(BTemp2-self%B0)*(BTemp4-self%B0)<0) exit
        end do
        if ((BTemp1-self%B0)*(BTemp3-self%B0)<0.and.(BTemp2-self%B0)*(BTemp4-self%B0)<0) then
            if ( feedback_level>1 ) write(*,'(a)') '   FAILURE as the root seems to be on both sides'
            return
        end if

        ! 7) Solve the equation B0(A)=B0_wanted.
        if ((BTemp1-self%B0)*(BTemp3-self%B0)<0) then
            realAp = zbrent(DesFR_BfuncA,ATemp1,ATemp3,1.d-50,self%B0,success)
            if (.not.success) then
                if ( feedback_level>1 ) write(*,'(a)') '   FAILURE right side solution not found'
                return
            end if
        else if ((BTemp2-self%B0)*(BTemp4-self%B0)<0) then
            realAp = zbrent(DesFR_BfuncA,ATemp2,ATemp4,1.d-50,self%B0,success)
            if (.not.success) then
                if ( feedback_level>1 ) write(*,'(a)') '   FAILURE left side solution not found'
                return
            end if
        else
            if ( feedback_level>1 ) write(*,'(a)') '   FAILURE the root was not on the right side nor the left one...'
            return
        end if

        if ( feedback_level>1 ) write(*,'(a,E13.4)') '   initial condition A  = ', realAp

        ! 8) Check if the result found is compatible with the requested one. This is required only for debug.
        BTemp1 = DesFR_BfuncA(realAp)
        if ( feedback_level>1 ) then
            write(*,'(a,E13.4)') '   B0 found = ', BTemp1
            write(*,'(a,E13.4)') '   B0 given = ', self%B0
        end if

        if (ABS(BTemp1-self%B0)/ABS(self%B0)>0.1_dl.and.ABS(BTemp1-self%B0)>1.d-8.and..false.) then
            if ( feedback_level>1 ) write(*,'(a)')  '   FAILURE designer code unable to find appropriate initial conditions'
            success = .false.
            return
        end if

        ! 9) output the value:
        A_ini = realAp

        return

    contains

        ! ---------------------------------------------------------------------------------------------
        !> Function that solves the designer f(R) equations and returns the present value of B0
        function DesFR_BfuncA(A)

            implicit none

            real(dl) :: A            !< input amplitude of the particular solution
            real(dl) :: DesFR_BfuncA !< value of B0
            real(dl) :: B0

            call self%solve_designer_equations( params_cache, A, B0, only_B0=.True., success=success )

            DesFR_BfuncA = B0

        end function

    end subroutine EFTCAMBDesignerFRFindInitialConditions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the number of parameters of the model.
    subroutine EFTCAMBDesignerFRComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_fR_designer)  :: self   !< the base class

        self%parameter_number = 1
        self%parameter_number = self%parameter_number +self%DesfRwDE%parameter_number

    end subroutine EFTCAMBDesignerFRComputeParametersNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints on the screen feedback information about the model.
    subroutine EFTCAMBDesignerFRFeedback( self, print_params )

        implicit none

        class(EFTCAMB_fR_designer)  :: self         !< the base class
        logical, optional           :: print_params !< optional flag that decised whether to print numerical values
                                                    !! of the parameters.

        write(*,*)
        write(*,'(a,a)')    '   Model               =  ', self%name
        write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number
        ! print model functions informations:
        if ( self%EFTwDE /= 0 ) then
            write(*,*)
            write(*,'(a,I3)')  '   EFTwDE              =', self%EFTwDE
        end if
        write(*,*)
        write(*,'(a24,F12.6)') '   B0                  =', self%B0

        call self%DesfRwDE%feedback( print_params )

    end subroutine EFTCAMBDesignerFRFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerFRParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_fR_designer)  :: self   !< the base class
        integer     , intent(in)    :: i      !< the index of the parameter
        character(*), intent(out)   :: name   !< the output name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! the first parameter is B0:
        else if ( i==1 ) then
            name = TRIM('B0')
            return
        ! the other parameters are the w_DE parameters:
        else
            call self%DesfRwDE%parameter_names( i-1, name )
            return
        end if

    end subroutine EFTCAMBDesignerFRParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerFRParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_fR_designer) :: self       !< the base class
        integer     , intent(in)    :: i         !< The index of the parameter
        character(*), intent(out)   :: latexname !< the output latex name of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! the first parameter is B0:
        else if ( i==1 ) then
            latexname = TRIM('B_0')
            return
        ! the other parameters are the w_DE parameters:
        else
            call self%DesfRwDE%parameter_names_latex( i-1, latexname )
            return
        end if

    end subroutine EFTCAMBDesignerFRParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name of the model
    subroutine EFTCAMBDesignerFRParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_fR_designer) :: self   !< the base class
        integer , intent(in)        :: i     !< The index of the parameter
        real(dl), intent(out)       :: value !< the output value of the i-th parameter

        ! check validity of input:
        if ( i<=0 .or. i>self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        ! the first parameter is B0:
        else if ( i==1 ) then
            value = self%B0
            return
        ! the other parameters are the w_DE parameters:
        else
            call self%DesfRwDE%parameter_value( i-1, value )
            return
        end if

    end subroutine EFTCAMBDesignerFRParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBDesignerFRBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                   :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: x, mu
        integer  :: ind

        x   = log(a)
        call self%EFTOmega%precompute(x, ind, mu )

        eft_cache%EFTOmegaV    = self%EFTOmega%value( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaP    = self%EFTOmega%first_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPP   = self%EFTOmega%second_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTOmegaPPP  = self%EFTOmega%third_derivative( x, index=ind, coeff=mu )
        eft_cache%EFTc         = 0._dl
        eft_cache%EFTLambda    = self%EFTLambda%value( x, index=ind, coeff=mu )
        eft_cache%EFTcdot      = 0._dl
        eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( x, index=ind, coeff=mu )

    end subroutine EFTCAMBDesignerFRBackgroundEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the value of the background EFT functions at a given time.
    subroutine EFTCAMBDesignerFRSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                   :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%EFTGamma1V  = 0._dl
        eft_cache%EFTGamma1P  = 0._dl
        eft_cache%EFTGamma2V  = 0._dl
        eft_cache%EFTGamma2P  = 0._dl
        eft_cache%EFTGamma3V  = 0._dl
        eft_cache%EFTGamma3P  = 0._dl
        eft_cache%EFTGamma4V  = 0._dl
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl
        eft_cache%EFTGamma5V  = 0._dl
        eft_cache%EFTGamma5P  = 0._dl
        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBDesignerFRSecondOrderEFTFunctions

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes dtauda = 1/sqrt(a^2H^2).For pure EFT std this has to be overridden
    !! for performance reasons.
    function EFTCAMBDesignerFRComputeDtauda( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                   :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: EFTCAMBDesignerFRComputeDtauda               !< the output dtauda

        real(dl) :: temp

        temp = eft_cache%grhoa2 +eft_par_cache%grhov*a*a*self%DesfRwDE%integral(a)
        EFTCAMBDesignerFRComputeDtauda = sqrt(3/temp)

    end function EFTCAMBDesignerFRComputeDtauda

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
    subroutine EFTCAMBDesignerFRComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                   :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%grhov_t = eft_par_cache%grhov*self%DesfRwDE%integral(a)
        eft_cache%adotoa  = sqrt( ( eft_cache%grhom_t +eft_cache%grhov_t )/3._dl )

    end subroutine EFTCAMBDesignerFRComputeAdotoa

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the two derivatives wrt conformal time of H.
    subroutine EFTCAMBDesignerFRComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                   :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        eft_cache%gpiv_t  = self%DesfRwDE%value(a)*eft_cache%grhov_t
        eft_cache%Hdot    = -0.5_dl*( eft_cache%adotoa**2 +eft_cache%gpresm_t +eft_cache%gpiv_t )
        eft_cache%Hdotdot = eft_cache%adotoa*( ( eft_cache%grhob_t +eft_cache%grhoc_t)/6._dl +2._dl*( eft_cache%grhor_t +eft_cache%grhog_t)/3._dl ) &
            & +eft_cache%adotoa*eft_cache%grhov_t*( 1._dl/6._dl +self%DesfRwDE%value(a) +1.5_dl*self%DesfRwDE%value(a)**2 -0.5_dl*a*self%DesfRwDE%first_derivative(a) ) &
            & +eft_cache%adotoa*eft_cache%grhonu_tot/6._dl -0.5_dl*eft_cache%adotoa*eft_cache%gpinu_tot -0.5_dl*eft_cache%gpinudot_tot

    end subroutine EFTCAMBDesignerFRComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBDesignerFRAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_fR_designer)                   :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBDesignerFRAdditionalModelStability          !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBDesignerFRAdditionalModelStability = .True.
        if ( self%DesfRwDE%value(a) > -1._dl/3._dl ) EFTCAMBDesignerFRAdditionalModelStability = .False.

    end function EFTCAMBDesignerFRAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_designer_fR

!----------------------------------------------------------------------------------------