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


module EFTCAMB_reconstructed_fit_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public reconstructed_fit_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: reconstructed_fit_parametrization_1D

        real(dl) :: P1
        real(dl) :: P2
        real(dl) :: P3
        real(dl) :: P4
        real(dl) :: P5
        real(dl) :: omegaL

    contains

        ! utility functions:
        procedure :: set_param_number      => ReconstructedFitParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
        procedure :: init_parameters       => ReconstructedFitParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => ReconstructedFitParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => ReconstructedFitParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => ReconstructedFitParametrized1DValue               !< function that returns the value of the Taylor expansion.
        procedure :: first_derivative      => ReconstructedFitParametrized1DFirstDerivative     !< function that returns the first derivative of the Taylor expansion.
        procedure :: second_derivative     => ReconstructedFitParametrized1DSecondDerivative    !< function that returns the second derivative of the Taylor expansion.
        procedure :: third_derivative      => ReconstructedFitParametrized1DThirdDerivative     !< function that returns the third derivative of the Taylor expansion.
        procedure :: integral              => ReconstructedFitParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type reconstructed_fit_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor expansion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine ReconstructedFitParametrized1DSetParamNumber( self )

        implicit none

        class(reconstructed_fit_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 6

    end subroutine ReconstructedFitParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine ReconstructedFitParametrized1DInitParams( self, array )

        implicit none

        class(reconstructed_fit_parametrization_1D)             :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%P1     = array(1)
        self%P2     = array(2)
        self%P3     = array(3)
        self%P4     = array(4)
        self%P5     = array(5)
        self%omegaL = array(6)

    end subroutine ReconstructedFitParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine ReconstructedFitParametrized1DParameterValues( self, i, value )

        implicit none

        class(reconstructed_fit_parametrization_1D)       :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%P1
            case(2)
                value = self%P2
            case(3)
                value = self%P3
            case(4)
                value = self%P4
            case(5)
                value = self%P5
            case(6)
                value = self%omegaL
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine ReconstructedFitParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine ReconstructedFitParametrized1DFeedback( self, print_params )

        implicit none

        class(reconstructed_fit_parametrization_1D) :: self         !< the base class
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

    end subroutine ReconstructedFitParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function ReconstructedFitParametrized1DValue( self, x, eft_cache )

        implicit none

        class(reconstructed_fit_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedFitParametrized1DValue                 !< the output value

        !> Maple output
        !> extra parameters
        real(dl) :: t13
        real(dl) :: t15
        real(dl) :: t2
        real(dl) :: t4
        real(dl) :: t5

        !> defining the extra parameters
        t2 = -x + self%P4
        t4 = tanh(self%P5 * t2)
        t5 = t4 ** 2
        t13 = t2 ** 2
        t15 = exp(-self%P3 * t13)

        !> value of the function
         ReconstructedFitParametrized1DValue = -0.2D1 / 0.3D1 * x * self%P2 * (-self%P5 * t5 / 0.2D1 - self%P3 * t2                     &
                                                * t4 + self%P5 / 0.2D1) * t15 + t5 * self%P1 * self%P5 * x / 0.3D1 - self%P1 * self%P5  &
                                                * x / 0.3D1 - 0.1D1




    end function ReconstructedFitParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function ReconstructedFitParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(reconstructed_fit_parametrization_1D)                     :: self      !< the base class
        real(dl), intent(in)                                            :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional              :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedFitParametrized1DFirstDerivative                    !< the output value

        !> Maple output
        !> some extra variables for optimization
        real(dl) :: t1
        real(dl) :: t12
        real(dl) :: t13
        real(dl) :: t15
        real(dl) :: t2
        real(dl) :: t23
        real(dl) :: t3
        real(dl) :: t36
        real(dl) :: t38
        real(dl) :: t5
        real(dl) :: t6

        !> Defining the extra variables
        t1 = self%P5 ** 2
        t2 = x * t1
        t3 = -x + self%P4
        t5 = tanh(self%P5 * t3)
        t6 = t5 ** 2
        t12 = x ** 2
        t13 = self%P3 * t12
        t15 = self%P5 * (-self%P3 * self%P4 * x + t13 - 0.1D1 / 0.4D1)
        t23 = self%P4 ** 2
        t36 = t3 ** 2
        t38 = exp(-self%P3 * t36)

        !> First Derivative
        ReconstructedFitParametrized1DFirstDerivative = (0.2D1 / 0.3D1 * t2 * t6 * t5 - 0.4D1 / 0.3D1                                       &
                                                        * t15 * t6 + 0.4D1 / 0.3D1 * (-t2 / 0.2D1 + (-0.2D1 * t13 * self%P4 +               &
                                                        t12 * x * self%P3 + (self%P3 * t23 - 0.1D1) * x + self%P4 / 0.2D1) * self%P3) * t5  &
                                                        + 0.4D1 / 0.3D1 * t15) * self%P2 * t38 + (0.2D1 / 0.3D1 * x * self%P5 * t5          &
                                                        + 0.1D1 / 0.3D1) * self%P1 * self%P5 * (t5 + 0.1D1) * (t5 - 0.1D1)



    end function ReconstructedFitParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function ReconstructedFitParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(reconstructed_fit_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedFitParametrized1DSecondDerivative      !< the output value

        !> Maple output
        !> extra varibales
        real(dl) :: t1
        real(dl) :: t11
        real(dl) :: t12
        real(dl) :: t13
        real(dl) :: t14
        real(dl) :: t20
        real(dl) :: t21
        real(dl) :: t25
        real(dl) :: t26
        real(dl) :: t32
        real(dl) :: t4
        real(dl) :: t40
        real(dl) :: t41
        real(dl) :: t6
        real(dl) :: t65
        real(dl) :: t67
        real(dl) :: t7
        real(dl) :: t70
        real(dl) :: t8

        !> defining the extra variables
        t1 = self%P5 ** 2
        t4 = -x + self%P4
        t6 = tanh(self%P5 * t4)
        t7 = t6 ** 2
        t8 = t7 ** 2
        t11 = self%P3 * self%P4
        t12 = t11 * x
        t13 = x ** 2
        t14 = self%P3 * t13
        t20 = x * t1
        t21 = t13 * x
        t25 = self%P4 ** 2
        t26 = self%P3 * t25
        t32 = 0.3D1 / 0.2D1 * (t21 * self%P3 - dble(2 * t14 * self%P4) + (t26 - &
                0.7D1 / 0.6D1) * x + 0.2D1 / 0.3D1 * self%P4) * self%P3
        t40 = t13 ** 2
        t41 = self%P3 ** 2
        t65 = t4 ** 2
        t67 = exp(-self%P3 * t65)
        t70 = x * self%P5

        !> second derivative
        ReconstructedFitParametrized1DSecondDerivative = -0.8D1 / 0.3D1 * self%P2 * (-0.3D1 / 0.4D1 * x *                           &
                                                        t1 * self%P5 * t8 + 0.3D1 / 0.2D1 * t1 * (-0.1D1 / 0.3D1 - t12 +            &
                                                        t14) * t7 * t6 - (-t20 + t32) * self%P5 * t7 + ((0.1D1 / 0.2D1 - 0.3D1      &
                                                        / 0.2D1 * t14 + 0.3D1 / 0.2D1 * t12) * t1 + (0.1D1 / 0.2D1                  &
                                                         + t40 * t41 - 0.3D1 * t21 * t41 * self%P4 + (-0.5D1 / 0.2D1 * self%P3 +    &
                                                        0.3D1 * t41 * t25) * t13 + (-t41 * t25 * self%P4 + 0.7D1 /                  &
                                                        0.2D1 * t11) * x - t26) * self%P3) * t6 + (-t20 / 0.4D1 + t32) * self%P5)   &
                                                        * t67 + (0.2D1 * t70 * t7 - 0.2D1 / 0.3D1 * t70 + 0.4D1 / 0.3D1 *           &
                                                        t6) * self%P1 * t1 * (t6 + 0.1D1) * (t6 - 0.1D1)


    end function ReconstructedFitParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function ReconstructedFitParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(reconstructed_fit_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedFitParametrized1DThirdDerivative       !< the output value

        !> Maple output
        !> extra variables
        real(dl) :: t1
        real(dl) :: t12
        real(dl) :: t14
        real(dl) :: t17
        real(dl) :: t2
        real(dl) :: t23
        real(dl) :: t25
        real(dl) :: t3
        real(dl) :: t32
        real(dl) :: t4
        real(dl) :: t40
        real(dl) :: t48
        real(dl) :: t55
        real(dl) :: t59
        real(dl) :: t6
        real(dl) :: t7
        real(dl) :: t8
        real(dl) :: t81
        real(dl) :: t83
        real(dl) :: t87

        !> defining the extra variables
        t1 = self%P5 ** 2
        t2 = t1 ** 2
        t3 = x * t2
        t4 = -x + self%P4
        t6 = tanh(self%P5 * t4)
        t7 = t6 ** 2
        t8 = t7 ** 2
        t12 = -t4
        t14 = x * t12 * self%P3
        t17 = t1 * self%P5
        t23 = t12 ** 2
        t25 = x * t23 * self%P3
        t32 = t7 * t6
        t40 = self%P3 ** 2
        t48 = (0.9D1 / 0.16D2 + x * t23 * t12 * t40 - (0.21D2 / 0.8D1 * &
            x - 0.9D1 / 0.8D1 * self%P4) * t12 * self%P3) * self%P3
        t55 = 0.9D1 / 0.4D1 * self%P4
        t59 = t23 ** 2
        t81 = t4 ** 2
        t83 = exp(-self%P3 * t81)
        t87 = x * self%P5

        !> Third derivative
        ReconstructedFitParametrized1DThirdDerivative = (0.8D1 * t3 * t8 * t6 - 0.16D2 / 0.3D1 * (-0.9D1                            &
                                                        / 0.8D1 + 0.3D1 * t14) * t17 * t8 + 0.16D2 * t1 * (-0.5D1 / 0.6D1           &
                                                        * x * t1 + (t25 - 0.5D1 / 0.4D1 * x + 0.3D1 / 0.4D1 * self%P4) * self%P3)   &
                                                        * t32 - 0.32D2 / 0.3D1 * self%P5 * ((0.3D1 / 0.4D1 - 0.2D1 * t14) *         &
                                                         t1 + t48) * t7 + 0.16D2 / 0.3D1 * (t3 - (0.3D1 * t25 - 0.15D2 / 0.4D1      &
                                                        * x + t55) * self%P3 * t1 + (x * t59 * t40 - (0.9D1 / 0.2D1 * x -           &
                                                         0.3D1 / 0.2D1 * self%P4) * t23 * self%P3 - t55 + 0.3D1 * x) * t40) * t6 +  &
                                                         0.32D2 / 0.3D1 * self%P5 * ((0.3D1 / 0.16D2 - t14 / 0.2D1) * t1 + t48))    &
                                                        * self%P2 * t83 + 0.8D1 * self%P1 * t17 * (t6 + 0.1D1) * (-0.1D1 / 0.4D1    &
                                                        + t87 * t32 - 0.2D1 / 0.3D1 * t87 * t6 + 0.3D1 / 0.4D1 * t7) * (t6 - 0.1D1)


    end function ReconstructedFitParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function ReconstructedFitParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(reconstructed_fit_parametrization_1D)        :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedFitParametrized1DIntegral              !< the output value

        !> Maple output
        !> extra variables
        real(dl) :: a
        real(dl) :: t1
        real(dl) :: t10
        real(dl) :: t12
        real(dl) :: t14
        real(dl) :: t22
        real(dl) :: t3
        real(dl) :: t4
        real(dl) :: t6
        real(dl) :: t9

        !> defining the extra variables
        t1 = -a + self%P4
        t3 = tanh(self%P5 * t1)
        t4 = t1 ** 2
        t6 = exp(-self%P3 * t4)
        t9 = -1 + self%P4
        t10 = t9 ** 2
        t12 = exp(-self%P3 * t10)
        t14 = tanh(self%P5 * t9)
        t22 = log((self%P2 * t12 * t14 - self%P2 * t3 * t6 + self%P1 * t14 - self%P1 * t3 + self%omegaL) / self%omegaL)

        !> Integral
        ReconstructedFitParametrized1DIntegral = -t22 / 0.3D1

    end function ReconstructedFitParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_reconstructed_fit_parametrizations_1D

!----------------------------------------------------------------------------------------
