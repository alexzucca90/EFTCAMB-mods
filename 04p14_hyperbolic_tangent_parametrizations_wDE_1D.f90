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

!> @file 04p10_hyperbolic_tangent_wDE_parametrizations_1D.f90
!! This file contains the definition of the reconstructed fit parametrization for the
!! dark energy density reconstruction, inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Taylor expansion parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

!> @author Alex Zucca: azucca@sfu.ca


module EFTCAMB_hyperbolic_tangent_wDE_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public hyperbolic_tangent_wDE_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: hyperbolic_tangent_wDE_parametrization_1D

        real(dl) :: A
        real(dl) :: B
        real(dl) :: C



    contains

        ! utility functions:
        procedure :: set_param_number      => HyperbolicTangentWDEParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
        procedure :: init_parameters       => HyperbolicTangentWDEParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => HyperbolicTangentWDEParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => HyperbolicTangentWDEParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => HyperbolicTangentWDEParametrized1DValue               !< function that returns the value of the Taylor expansion.
        procedure :: first_derivative      => HyperbolicTangentWDEParametrized1DFirstDerivative     !< function that returns the first derivative of the Taylor expansion.
        procedure :: second_derivative     => HyperbolicTangentWDEParametrized1DSecondDerivative    !< function that returns the second derivative of the Taylor expansion.
        procedure :: third_derivative      => HyperbolicTangentWDEParametrized1DThirdDerivative     !< function that returns the third derivative of the Taylor expansion.
        procedure :: integral              => HyperbolicTangentWDEParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type hyperbolic_tangent_wDE_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor expansion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine HyperbolicTangentWDEParametrized1DSetParamNumber( self )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 3

    end subroutine HyperbolicTangentWDEParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine HyperbolicTangentWDEParametrized1DInitParams( self, array )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D)             :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%A     = array(1)
        self%B     = array(2)
        self%C     = array(3)

    end subroutine HyperbolicTangentWDEParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine HyperbolicTangentWDEParametrized1DParameterValues( self, i, value )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D)       :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%A
            case(2)
                value = self%B
            case(3)
                value = self%C
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop( 'EFTCAMB error' )
        end select

    end subroutine HyperbolicTangentWDEParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine HyperbolicTangentWDEParametrized1DFeedback( self, print_params )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D) :: self         !< the base class
        logical, optional                           :: print_params !< optional flag that decised whether to print numerical values
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

        write(*,*)     'Reconstruction Fit parametrization: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine HyperbolicTangentWDEParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function HyperbolicTangentWDEParametrized1DValue( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                                    :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional      :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentWDEParametrized1DValue                  !< the output value

        real(dl) :: t5
        real(dl) :: t6
        real(dl) :: t12

        t5 = tanh(self%B * (x - self%C))
        t6 = t5 ** 2
        t12 = tanh(self%B * (1._dl - self%C))

        !> value of the function
        HyperbolicTangentWDEParametrized1DValue = -0.3D1 / x * self%A * self%B * (1._dl - t6) / (-self%A * t12 + self%A * t5 + 1._dl) - 1._dl


    end function HyperbolicTangentWDEParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function HyperbolicTangentWDEParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D)                :: self      !< the base class
        real(dl), intent(in)                                            :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional              :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentWDEParametrized1DFirstDerivative                    !< the output value

        real(dl) :: t1
        real(dl) :: t6
        real(dl) :: t7
        real(dl) :: t8
        real(dl) :: t13
        real(dl) :: t15
        real(dl) :: t16
        real(dl) :: t20
        real(dl) :: t22
        real(dl) :: t28
        real(dl) :: t30
        real(dl) :: t32

        t1 = x ** 2
        t6 = tanh(self%B * (x - self%C))
        t7 = t6 ** 2
        t8 = 1._dl - t7
        t13 = tanh(self%B * (1._dl - self%C))
        t15 = -self%A * t13 + self%A * t6 + 1._dl
        t16 = 1._dl / t15
        t20 = 1._dl / x
        t22 = self%B ** 2
        t28 = self%A ** 2
        t30 = t8 ** 2
        t32 = t15 ** 2

        !> first derivative
        HyperbolicTangentWDEParametrized1DFirstDerivative = 3._dl / t1 * self%A * self%B * t8 * t16 + 6._dl * t20 * self%A * t22 * t6 * t8 * t16 &
                    &+ 3._dl * t20 * t28 * t22 * t30 / t32

    end function HyperbolicTangentWDEParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function HyperbolicTangentWDEParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D)    :: self         !< the base class
        real(dl), intent(in)                                :: x            !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional  :: eft_cache    !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentWDEParametrized1DSecondDerivative      !< the output value

        real(dl) :: t1
        real(dl) :: t7
        real(dl) :: t8
        real(dl) :: t9
        real(dl) :: t14
        real(dl) :: t16
        real(dl) :: t17
        real(dl) :: t21
        real(dl) :: t23
        real(dl) :: t29
        real(dl) :: t31
        real(dl) :: t33
        real(dl) :: t34
        real(dl) :: t38
        real(dl) :: t39
        real(dl) :: t40

        t1 = x ** 2
        t7 = tanh(self%B * (x - self%C))
        t8 = t7 ** 2
        t9 = 1._dl - t8
        t14 = tanh(self%B * (1._dl - self%C))
        t16 = -self%A * t14 + self%A * t7 + 1._dl
        t17 = 1._dl / t16
        t21 = 1._dl / t1
        t23 = self%B ** 2
        t29 = self%A ** 2
        t31 = t9 ** 2
        t33 = t16 ** 2
        t34 = 1._dl / t33
        t38 = 1._dl / x
        t39 = t38 * self%A
        t40 = t23 * self%B

        !> second derivative
         HyperbolicTangentWDEParametrized1DSecondDerivative = -6._dl / t1 / x * self%A * self%B * t9 * t17 - 0.12D2 * t21 * self%A * t23 * t7 * t9 * t17 &
                & - 6._dl * t21 * t29 * t23 * t31 * t34 + 6._dl * t39 * t40 * t31 * t17 - 12._dl * t39 * t40 * t8 * t9 * t17 &
                & - 18._dl * t38 * t29 * t40 * t7 * t31 * t34 - 6._dl * t38 * t29 * self%A * t40 * t31 * t9 / t33 / t16


    end function HyperbolicTangentWDEParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function HyperbolicTangentWDEParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D)    :: self      !< the base class
        real(dl), intent(in)                                :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional  :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentWDEParametrized1DThirdDerivative    !< the output value

        real(dl) :: t1
        real(dl) :: t2
        real(dl) :: t7
        real(dl) :: t8
        real(dl) :: t9
        real(dl) :: t14
        real(dl) :: t16
        real(dl) :: t17
        real(dl) :: t22
        real(dl) :: t24
        real(dl) :: t30
        real(dl) :: t32
        real(dl) :: t34
        real(dl) :: t35
        real(dl) :: t39
        real(dl) :: t40
        real(dl) :: t41
        real(dl) :: t57
        real(dl) :: t59
        real(dl) :: t62
        real(dl) :: t66
        real(dl) :: t68
        real(dl) :: t69
        real(dl) :: t74
        real(dl) :: t95
        real(dl) :: t97
        real(dl) :: t99

        t1 = x ** 2
        t2 = t1 ** 2
        t7 = tanh(self%B * (x - self%C))
        t8 = t7 ** 2
        t9 = 1._dl - t8
        t14 = tanh(self%B * (1._dl - self%C))
        t16 = -self%A * t14 + self%A * t7 + 1._dl
        t17 = 1._dl / t16
        t22 = 1._dl / t1 / x
        t24 = self%B ** 2
        t30 = self%A ** 2
        t32 = t9 ** 2
        t34 = t16 ** 2
        t35 = 1._dl / t34
        t39 = 1._dl / t1
        t40 = t39 * self%A
        t41 = t24 * self%B
        t57 = t30 * self%A
        t59 = t32 * t9
        t62 = 1._dl / t34 / t16
        t66 = 1._dl / x
        t68 = t24 ** 2
        t69 = t66 * self%A * t68
        t74 = t66 * t30
        t95 = t30 ** 2
        t97 = t32 ** 2
        t99 = t34 ** 2

        !> third derivative
        HyperbolicTangentWDEParametrized1DThirdDerivative = 18._dl / t2 * self%A * self%B * t9 * t17 + 36._dl * t22 * self%A * t24 * t7 * t9 * t17  &
                & + 18._dl * t22 * t30 * t24 * t32 * t35 - 18._dl * t40 * t41 * t32 * t17 + 36._dl * t40 * t41 * t8 * t9 * t17                      &
                & + 54._dl* t39 * t30 * t41 * t7 * t32 * t35 + 18._dl * t39 * t57* t41 * t59 * t62 - 48._dl * t69 * t32 * t17 * t7                  &
                & - 24._dl * t74 * t68 * t59 * t35 + 24._dl * t69 * t8 * t7 * t9 * t17 + 84._dl * t74 * t68 * t8 * t32 * t35                        &
                & + 72._dl * t66 * t57 * t68 * t7 * t59 * t62 + 18._dl * t66 * t95 * t68 * t97 / t99


    end function HyperbolicTangentWDEParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function HyperbolicTangentWDEParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(hyperbolic_tangent_wDE_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: HyperbolicTangentWDEParametrized1DIntegral              !< the output value

        !> Integral
        HyperbolicTangentWDEParametrized1DIntegral = (self%A*tanh(self%B*(x-self%C)) + 1._dl - self%A*tanh(self%B*(1._dl-self%C)))*x**2

    end function HyperbolicTangentWDEParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_hyperbolic_tangent_wDE_parametrizations_1D

!----------------------------------------------------------------------------------------
