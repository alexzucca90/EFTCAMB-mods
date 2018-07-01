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

!> @file 04p8_taylor_expansion_parametrizations_1D.f90
!! This file contains the definition of the Taylor expansion parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Taylor expansion parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_polynomial_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public polynomial_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: polynomial_parametrization_1D

        real(dl) :: alpha0
        real(dl) :: alpha1
        real(dl) :: alpha2
        real(dl) :: alpha3
        real(dl) :: alpha4
        real(dl) :: alpha5

    contains

        ! utility functions:
        procedure :: set_param_number      => PolynomialParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
        procedure :: init_parameters       => PolynomialParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => PolynomialParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => PolynomialParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => PolynomialParametrized1DValue               !< function that returns the value of the Taylor expansion.
        procedure :: first_derivative      => PolynomialParametrized1DFirstDerivative     !< function that returns the first derivative of the Taylor expansion.
        procedure :: second_derivative     => PolynomialParametrized1DSecondDerivative    !< function that returns the second derivative of the Taylor expansion.
        procedure :: third_derivative      => PolynomialParametrized1DThirdDerivative     !< function that returns the third derivative of the Taylor expansion.
        procedure :: integral              => PolynomialParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type polynomial_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor expansion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine PolynomialParametrized1DSetParamNumber( self )

        implicit none

        class(polynomial_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 6

    end subroutine PolynomialParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine PolynomialParametrized1DInitParams( self, array )

        implicit none

        class(polynomial_parametrization_1D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%alpha0 = array(1)
        self%alpha1 = array(2)
        self%alpha2 = array(3)
        self%alpha3 = array(4)
        self%alpha4 = array(5)
        self%alpha5 = array(6)

    end subroutine PolynomialParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine PolynomialParametrized1DParameterValues( self, i, value )

        implicit none

        class(polynomial_parametrization_1D)       :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%alpha0
            case(2)
                value = self%alpha1
            case(3)
                value = self%alpha2
            case(4)
                value = self%alpha3
            case(5)
                value = self%alpha4
            case(6)
                value = self%alpha5
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine PolynomialParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine PolynomialParametrized1DFeedback( self, print_params )

        implicit none

        class(polynomial_parametrization_1D) :: self         !< the base class
        logical, optional                :: print_params !< optional flag that decised whether to print numerical values
                                                         !! of the parameters.

        integer                                 :: i
        real(dl)                                :: param_value
        character(len=EFT_names_max_length)     :: param_name
        logical                                 :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,*)     'Taylor expansion parametrization: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine PolynomialParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function PolynomialParametrized1DValue( self, x, eft_cache )

        implicit none

        class(polynomial_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PolynomialParametrized1DValue                           !< the output value

        PolynomialParametrized1DValue = self%alpha0 +self%alpha1*x +self%alpha2*x**2 +self%alpha3*x**3 &
                & + self%alpha4*x**4 + self%alpha5*x**5

    end function PolynomialParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function PolynomialParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(polynomial_parametrization_1D)               :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PolynomialParametrized1DFirstDerivative                 !< the output value

        PolynomialParametrized1DFirstDerivative = self%alpha1+ 2._dl*self%alpha2*x + 3._dl*self%alpha3*x**2 &
                & + 4._dl * self%alpha4*x**3 + 5._dl*self%alpha5 * x**4


    end function PolynomialParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function PolynomialParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(polynomial_parametrization_1D)               :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PolynomialParametrized1DSecondDerivative                !< the output value

        PolynomialParametrized1DSecondDerivative = 2._dl*self%alpha2+6._dl*self%alpha3*x + 12._dl * self%alpha4 * x**2 &
                & + 20._dl * self%alpha5*x**3

    end function PolynomialParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function PolynomialParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(polynomial_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PolynomialParametrized1DThirdDerivative                 !< the output value

        PolynomialParametrized1DThirdDerivative = 6._dl * self%alpha3 + 24._dl * self%alpha4*x + 60._dl * self%alpha5*x**2

    end function PolynomialParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function PolynomialParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(polynomial_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: PolynomialParametrized1DIntegral                        !< the output value

        PolynomialParametrized1DIntegral = 0.d0
        write(*,*) "WARNING: integral of the 5th order polynomial called"

    end function PolynomialParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_polynomial_parametrizations_1D

!----------------------------------------------------------------------------------------
