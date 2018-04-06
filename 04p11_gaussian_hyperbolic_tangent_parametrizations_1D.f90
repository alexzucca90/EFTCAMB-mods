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

!> @file 04p10_reconstructed_fit_parametrizations_1D.f90
!! This file contains the definition of the reconstructed fit parametrization for the
!! dark energy density reconstruction, inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Taylor expansion parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

!> @author Alex Zucca: azucca@sfu.ca


module EFTCAMB_gaussian_hyperbolic_tangent_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public gaussian_hyperbolic_tangent_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: gaussian_hyperbolic_tangent_parametrization_1D

        real(dl) :: A
        real(dl) :: B
        real(dl) :: C
        real(dl) :: D
        real(dl) :: E

    contains

        ! utility functions:
        procedure :: set_param_number      => GaussianHyperbolicTangentParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
        procedure :: init_parameters       => GaussianHyperbolicTangentParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => GaussianHyperbolicTangentParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => GaussianHyperbolicTangentParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => GaussianHyperbolicTangentParametrized1DValue               !< function that returns the value of the Taylor expansion.
        procedure :: first_derivative      => GaussianHyperbolicTangentParametrized1DFirstDerivative     !< function that returns the first derivative of the Taylor expansion.
        procedure :: second_derivative     => GaussianHyperbolicTangentParametrized1DSecondDerivative    !< function that returns the second derivative of the Taylor expansion.
        procedure :: third_derivative      => GaussianHyperbolicTangentParametrized1DThirdDerivative     !< function that returns the third derivative of the Taylor expansion.
        procedure :: integral              => GaussianHyperbolicTangentParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type gaussian_hyperbolic_tangent_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor expansion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine GaussianHyperbolicTangentParametrized1DSetParamNumber( self )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 5

    end subroutine GaussianHyperbolicTangentParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine GaussianHyperbolicTangentParametrized1DInitParams( self, array )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D)          :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%A     = array(1)
        self%B     = array(2)
        self%C     = array(3)
        self%D     = array(4)
        self%E     = array(5)

    end subroutine GaussianHyperbolicTangentParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine GaussianHyperbolicTangentParametrized1DParameterValues( self, i, value )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D)       :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%A
            case(2)
                value = self%B
            case(3)
                value = self%C
            case(4)
                value = self%D
            case(5)
                value = self%E
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine GaussianHyperbolicTangentParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine GaussianHyperbolicTangentParametrized1DFeedback( self, print_params )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D)   :: self         !< the base class
        logical, optional                                       :: print_params !< optional flag that decised whether to print numerical values
                                                         !! of the parameters.

        integer                                     :: i
        real(dl)                                    :: param_value
        character(len=EFT_names_max_length)         :: param_name
        logical                                     :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,*)     'Reconstruction DE Fit parametrization: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine GaussianHyperbolicTangentParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function GaussianHyperbolicTangentParametrized1DValue( self, x, eft_cache )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D)       :: self      !< the base class
        real(dl), intent(in)                                        :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional          :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: GaussianHyperbolicTangentParametrized1DValue                 !< the output value

        real(dl) :: t1
        real(dl) :: t2
        real(dl) :: t4
        real(dl) :: t8
        real(dl) :: t10
        real(dl) :: t11
        real(dl) :: t13
        real(dl) :: t17

        t1 = x - self%D
        t2 = t1 ** 2
        t4 = exp(-self%C * t2)
        t8 = tanh(self%E * t1)
        t10 = 1._dl - self%D
        t11 = t10 ** 2
        t13 = exp(-self%C * t11)
        t17 = tanh(self%E * t10)

        !> value of the function
        GaussianHyperbolicTangentParametrized1DValue = (self%B * t4 + self%A) * t8 + 1._dl - (self%B * t13 + self%A) * t17



    end function GaussianHyperbolicTangentParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function GaussianHyperbolicTangentParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D)           :: self      !< the base class
        real(dl), intent(in)                                            :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional              :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: GaussianHyperbolicTangentParametrized1DFirstDerivative           !< the output value


        !> Maple output
        !> extra varibales
        real(dl) :: t2
        real(dl) :: t3
        real(dl) :: t5
        real(dl) :: t8
        real(dl) :: t15

        t2 = x - self%D
        t3 = t2 ** 2
        t5 = exp(-self%C * t3)
        t8 = tanh(self%E * t2)
        t15 = t8 ** 2

        !> First derivative
        GaussianHyperbolicTangentParametrized1DFirstDerivative = -2._dl * self%B * self%C * t2 * t5 * t8 + (self%B * t5 + self%A) * self%E * (1._dl - t15)


    end function GaussianHyperbolicTangentParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function GaussianHyperbolicTangentParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D)       :: self      !< the base class
        real(dl), intent(in)                                        :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional          :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: GaussianHyperbolicTangentParametrized1DSecondDerivative      !< the output value

        !> Maple output
        !> extra varibales
        real(dl) :: t1
        real(dl) :: t2
        real(dl) :: t3
        real(dl) :: t5
        real(dl) :: t7
        real(dl) :: t11
        real(dl) :: t19
        real(dl) :: t20
        real(dl) :: t26

        t1 = self%B * self%C
        t2 = x - self%D
        t3 = t2 ** 2
        t5 = exp(-self%C * t3)
        t7 = tanh(self%E * t2)
        t11 = self%C ** 2
        t19 = t7 ** 2
        t20 = 1._dl - t19
        t26 = self%E ** 2

        !> Second Derivative
        GaussianHyperbolicTangentParametrized1DSecondDerivative = -2._dl * t1 * t5 * t7 + 4._dl * self%B * t11 * t3 * t5 * t7 - 4._dl * t1 * t2 * t5 * self%E * t20 &
                    &- 2._dl * (self%B * t5 + self%A) * t26 * t7 * t20

    end function GaussianHyperbolicTangentParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function GaussianHyperbolicTangentParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D)     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: GaussianHyperbolicTangentParametrized1DThirdDerivative       !< the output value

        !> Maple output
        !> extra variables
        real(dl) :: t1
        real(dl) :: t2
        real(dl) :: t3
        real(dl) :: t4
        real(dl) :: t6
        real(dl) :: t9
        real(dl) :: t13
        real(dl) :: t15
        real(dl) :: t16
        real(dl) :: t17
        real(dl) :: t31
        real(dl) :: t40
        real(dl) :: t41

        t1 = self%C ** 2
        t2 = self%B * t1
        t3 = x - self%D
        t4 = t3 ** 2
        t6 = exp(-self%C * t4)
        t9 = tanh(self%E * t3)
        t13 = self%B * self%C
        t15 = t9 ** 2
        t16 = 1._dl - t15
        t17 = t6 * self%E * t16
        t31 = self%E ** 2
        t40 = (self%B * t6 + self%A) * t31 * self%E
        t41 = t16 ** 2

        !> third derivative
        GaussianHyperbolicTangentParametrized1DThirdDerivative = -8._dl * self%B * t1 * self%C * t4 * t3 * t6 * t9 + 12._dl* t13 * t3 * t6 * t31 * t9 * t16 &
                    & + 12._dl * t2 * t3 * t6 * t9 + 4._dl * t40 * t15 * t16 + 12._dl * t2 * t4 * t17 - 6._dl * t13 * t17 - 2._dl * t40 * t41

    end function GaussianHyperbolicTangentParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function GaussianHyperbolicTangentParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(gaussian_hyperbolic_tangent_parametrization_1D)       :: self      !< the base class
        real(dl), intent(in)                                        :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional          :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: GaussianHyperbolicTangentParametrized1DIntegral              !< the output value

        write(*,*) "WARNING: the function GaussianHyperbolicTangent%Integral is not supposed to be called. Returning 0"
        GaussianHyperbolicTangentParametrized1DIntegral = 0._dl

    end function GaussianHyperbolicTangentParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_gaussian_hyperbolic_tangent_parametrizations_1D

!----------------------------------------------------------------------------------------
