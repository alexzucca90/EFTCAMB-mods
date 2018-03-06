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
        !real(dl) :: t13
        !real(dl) :: t15
        !real(dl) :: t2
        !real(dl) :: t4
        !real(dl) :: t5

        !> defining the extra parameters
        !t2 = -x + self%P4
        !t4 = tanh(self%P5 * t2)
        !t5 = t4 ** 2
        !t13 = t2 ** 2
        !t15 = exp(-self%P3 * t13)

        !> value of the function
        !ReconstructedFitParametrized1DValue = -0.2D1 / 0.3D1 * x * self%P2 * (-self%P5 * t5 / 0.2D1 - self%P3 * t2                     &
        !                                        * t4 + self%P5 / 0.2D1) * t15 + t5 * self%P1 * self%P5 * x / 0.3D1 - self%P1 * self%P5  &
        !                                        * x / 0.3D1 - 0.1D1


        !> MAPLE output.
        !> extra variables for optimization
        real(dl) :: t2
        real(dl) :: t3
        real(dl) :: t5
        real(dl) :: t8
        real(dl) :: t13
        real(dl) :: t15
        real(dl) :: t21
        real(dl) :: t22
        real(dl) :: t24
        real(dl) :: t28


        !> defining the extra variables
        t2 = x - self%P4
        t3 = t2 ** 2
        t5 = exp(-self%P3 * t3)
        t8 = tanh(self%P5 * t2)
        t13 = self%P2 * t5 + self%P1
        t15 = t8 ** 2
        t21 = 1._dl - self%P4
        t22 = t21 ** 2
        t24 = exp(-self%P3 * t22)
        t28 = tanh(self%P5 * t21)

        !> value of the function
        ReconstructedFitParametrized1DValue = -x * (-0.2D1 * self%P2 * self%P3 * t2 * t5 * t8 + t13 * self%P5 * (0.1D1 - t15))/&
                                                (t13 * t8 + self%omegaL - (self%P2 * t24 + self%P1) * t28) / 0.3D1 - 0.1D1

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
        !real(dl) :: t1
        !real(dl) :: t12
        !real(dl) :: t13
        !real(dl) :: t15
        !real(dl) :: t2
        !real(dl) :: t23
        !real(dl) :: t3
        !real(dl) :: t36
        !real(dl) :: t38
        !real(dl) :: t5
        !real(dl) :: t6

        !> Defining the extra variables
        !t1 = self%P5 ** 2
        !t2 = x * t1
        !t3 = -x + self%P4
        !t5 = tanh(self%P5 * t3)
        !t6 = t5 ** 2
        !t12 = x ** 2
        !t13 = self%P3 * t12
        !t15 = self%P5 * (-self%P3 * self%P4 * x + t13 - 0.1D1 / 0.4D1)
        !t23 = self%P4 ** 2
        !t36 = t3 ** 2
        !t38 = exp(-self%P3 * t36)

        !> First Derivative
        !ReconstructedFitParametrized1DFirstDerivative = (0.2D1 / 0.3D1 * t2 * t6 * t5 - 0.4D1 / 0.3D1                                       &
        !                                                * t15 * t6 + 0.4D1 / 0.3D1 * (-t2 / 0.2D1 + (-0.2D1 * t13 * self%P4 +               &
        !                                                t12 * x * self%P3 + (self%P3 * t23 - 0.1D1) * x + self%P4 / 0.2D1) * self%P3) * t5  &
        !                                                + 0.4D1 / 0.3D1 * t15) * self%P2 * t38 + (0.2D1 / 0.3D1 * x * self%P5 * t5          &
        !                                                + 0.1D1 / 0.3D1) * self%P1 * self%P5 * (t5 + 0.1D1) * (t5 - 0.1D1)


        !> Maple output
        !> extra parameters
        real(dl) :: t1
        real(dl) :: t2
        real(dl) :: t3
        real(dl) :: t5
        real(dl) :: t8
        real(dl) :: t13
        real(dl) :: t15
        real(dl) :: t16
        real(dl) :: t18
        real(dl) :: t20
        real(dl) :: t21
        real(dl) :: t23
        real(dl) :: t27
        real(dl) :: t29
        real(dl) :: t30
        real(dl) :: t35
        real(dl) :: t46
        real(dl) :: t54
        real(dl) :: t56

        !> Defining the extra variables
        t1 = self%P2 * self%P3
        t2 = x - self%P4
        t3 = t2 ** 2
        t5 = exp(-self%P3 * t3)
        t8 = tanh(self%P5 * t2)
        t13 = self%P2 * t5 + self%P1
        t15 = t8 ** 2
        t16 = 0.1D1 - t15
        t18 = -0.2D1 * t1 * t2 * t5 * t8 + t13 * self%P5 * t16
        t20 = 1._dl - self%P4
        t21 = t20 ** 2
        t23 = exp(-self%P3 * t21)
        t27 = tanh(self%P5 * t20)
        t29 = t13 * t8 + self%omegaL - (self%P2 * t23 + self%P1) * t27
        t30 = 0.1D1 / t29
        t35 = self%P3 ** 2
        t46 = self%P5 ** 2
        t54 = t18 ** 2
        t56 = t29 ** 2

        !> first derivative
        ReconstructedFitParametrized1DFirstDerivative = -t18 * t30 / 0.3D1 - x * (0.4D1 * self%P2 * t35 * t3 * t5 * t8 - 0.4D1 * t1 * t2 * t5 * self%P5 * t16 - &
                                                        0.2D1 * t13 * t46* t8 * t16 - 0.2D1 * t1 * t5 * t8) * t30 / 0.3D1 + x * t54 / t56 / 0.3D1



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
        !real(dl) :: t1
        !real(dl) :: t11
        !real(dl) :: t12
        !real(dl) :: t13
        !real(dl) :: t14
        !real(dl) :: t20
        !real(dl) :: t21
        !real(dl) :: t25
        !real(dl) :: t26
        !real(dl) :: t32
        !real(dl) :: t4
        !real(dl) :: t40
        !real(dl) :: t41
        !real(dl) :: t6
        !real(dl) :: t65
        !real(dl) :: t67
        !real(dl) :: t7
        !real(dl) :: t70
        !real(dl) :: t8

        !> defining the extra variables
        !t1 = self%P5 ** 2
        !t4 = -x + self%P4
        !t6 = tanh(self%P5 * t4)
        !t7 = t6 ** 2
        !t8 = t7 ** 2
        !t11 = self%P3 * self%P4
        !t12 = t11 * x
        !t13 = x ** 2
        !t14 = self%P3 * t13
        !t20 = x * t1
        !t21 = t13 * x
        !t25 = self%P4 ** 2
        !t26 = self%P3 * t25
        !t32 = 0.3D1 / 0.2D1 * (t21 * self%P3 - dble(2 * t14 * self%P4) + (t26 - &
        !        0.7D1 / 0.6D1) * x + 0.2D1 / 0.3D1 * self%P4) * self%P3
        !t40 = t13 ** 2
        !t41 = self%P3 ** 2
        !t65 = t4 ** 2
        !t67 = exp(-self%P3 * t65)
        !t70 = x * self%P5

        !> second derivative
        !ReconstructedFitParametrized1DSecondDerivative = -0.8D1 / 0.3D1 * self%P2 * (-0.3D1 / 0.4D1 * x *                           &
        !                                                t1 * self%P5 * t8 + 0.3D1 / 0.2D1 * t1 * (-0.1D1 / 0.3D1 - t12 +            &
        !                                                t14) * t7 * t6 - (-t20 + t32) * self%P5 * t7 + ((0.1D1 / 0.2D1 - 0.3D1      &
        !                                                / 0.2D1 * t14 + 0.3D1 / 0.2D1 * t12) * t1 + (0.1D1 / 0.2D1                  &
        !                                                 + t40 * t41 - 0.3D1 * t21 * t41 * self%P4 + (-0.5D1 / 0.2D1 * self%P3 +    &
        !                                                0.3D1 * t41 * t25) * t13 + (-t41 * t25 * self%P4 + 0.7D1 /                  &
        !                                                0.2D1 * t11) * x - t26) * self%P3) * t6 + (-t20 / 0.4D1 + t32) * self%P5)   &
        !                                                * t67 + (0.2D1 * t70 * t7 - 0.2D1 / 0.3D1 * t70 + 0.4D1 / 0.3D1 *           &
        !                                                t6) * self%P1 * t1 * (t6 + 0.1D1) * (t6 - 0.1D1)


        !> MAPLE output
        !> extra variables
        real(dl) :: t1
        real(dl) :: t2
        real(dl) :: t3
        real(dl) :: t5
        real(dl) :: t7
        real(dl) :: t11
        real(dl) :: t12
        real(dl) :: t17
        real(dl) :: t19
        real(dl) :: t20
        real(dl) :: t21
        real(dl) :: t25
        real(dl) :: t26
        real(dl) :: t28
        real(dl) :: t31
        real(dl) :: t33
        real(dl) :: t34
        real(dl) :: t36
        real(dl) :: t40
        real(dl) :: t42
        real(dl) :: t43
        real(dl) :: t47
        real(dl) :: t52
        real(dl) :: t53
        real(dl) :: t54
        real(dl) :: t55
        real(dl) :: t77
        real(dl) :: t78

        !> defining the extra variables
        t1 = self%P2 * self%P3
        t2 = x - self%P4
        t3 = t2 ** 2
        t5 = exp(-self%P3 * t3)
        t7 = tanh(self%P5 * t2)
        t11 = self%P3 ** 2
        t12 = self%P2 * t11
        t17 = t1 * t2
        t19 = t7 ** 2
        t20 = 0.1D1 - t19
        t21 = t5 * self%P5 * t20
        t25 = self%P2 * t5 + self%P1
        t26 = self%P5 ** 2
        t28 = t7 * t20
        t31 = 0.4D1 * t12 * t3 * t5 * t7 - 0.2D1 * t1 * t5 * t7 - 0.2D1 *t25 * t26 * t28 - 0.4D1 * t17 * t21
        t33 = 1._dl - self%P4
        t34 = t33 ** 2
        t36 = exp(-self%P3 * t34)
        t40 = tanh(self%P5 * t33)
        t42 = t25 * t7 + self%omegaL - (self%P2 * t36 + self%P1) * t40
        t43 = 0.1D1 / t42
        t47 = t2 * t5 * t7
        t52 = t25 * self%P5 * t20 - 0.2D1 * t1 * t47
        t53 = t52 ** 2
        t54 = t42 ** 2
        t55 = 0.1D1 / t54
        t77 = t25 * t26 * self%P5
        t78 = t20 ** 2

        !> second derivative
        ReconstructedFitParametrized1DSecondDerivative = -0.2D1 / 0.3D1 * t31 * t43 + 0.2D1 / 0.3D1 * t53 * t55 - x *       &
                                                        (-0.8D1 * self%P2 * t11 * self%P3 * t3 * t2 * t5 * t7 + 0.12D2      &
                                                        * t17 * t5 * t26 * t28 + 0.12D2 * t12 * t3 * t21 + 0.4D1 * t77*     &
                                                        t19 * t20 - 0.6D1 * t1 * t21 + 0.12D2 * t12 * t47 - 0.2D1 * t77*    &
                                                        t78) * t43 / 0.3D1 + x * t31 * t55 * t52 - 0.2D1 / 0.3D1 * x * t53  &
                                                        * t52 / t54 / t42


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
        !real(dl) :: t1
        !real(dl) :: t12
        !real(dl) :: t14
        !real(dl) :: t17
        !real(dl) :: t2
        !real(dl) :: t23
        !real(dl) :: t25
        !real(dl) :: t3
        !real(dl) :: t32
        !real(dl) :: t4
        !real(dl) :: t40
        !real(dl) :: t48
        !real(dl) :: t55
        !real(dl) :: t59
        !real(dl) :: t6
        !real(dl) :: t7
        !real(dl) :: t8
        !real(dl) :: t81
        !real(dl) :: t83
        !real(dl) :: t87

        !> defining the extra variables
        !t1 = self%P5 ** 2
        !t2 = t1 ** 2
        !t3 = x * t2
        !t4 = -x + self%P4
        !t6 = tanh(self%P5 * t4)
        !t7 = t6 ** 2
        !t8 = t7 ** 2
        !t12 = -t4
        !t14 = x * t12 * self%P3
        !t17 = t1 * self%P5
        !t23 = t12 ** 2
        !t25 = x * t23 * self%P3
        !t32 = t7 * t6
        !t40 = self%P3 ** 2
        !t48 = (0.9D1 / 0.16D2 + x * t23 * t12 * t40 - (0.21D2 / 0.8D1 * &
        !    x - 0.9D1 / 0.8D1 * self%P4) * t12 * self%P3) * self%P3
        !t55 = 0.9D1 / 0.4D1 * self%P4
        !t59 = t23 ** 2
        !t81 = t4 ** 2
        !t83 = exp(-self%P3 * t81)
        !t87 = x * self%P5

        !> Third derivative
        !ReconstructedFitParametrized1DThirdDerivative = (0.8D1 * t3 * t8 * t6 - 0.16D2 / 0.3D1 * (-0.9D1                            &
        !                                                / 0.8D1 + 0.3D1 * t14) * t17 * t8 + 0.16D2 * t1 * (-0.5D1 / 0.6D1           &
        !                                                * x * t1 + (t25 - 0.5D1 / 0.4D1 * x + 0.3D1 / 0.4D1 * self%P4) * self%P3)   &
        !                                                * t32 - 0.32D2 / 0.3D1 * self%P5 * ((0.3D1 / 0.4D1 - 0.2D1 * t14) *         &
        !                                                 t1 + t48) * t7 + 0.16D2 / 0.3D1 * (t3 - (0.3D1 * t25 - 0.15D2 / 0.4D1      &
        !                                                * x + t55) * self%P3 * t1 + (x * t59 * t40 - (0.9D1 / 0.2D1 * x -           &
        !                                                 0.3D1 / 0.2D1 * self%P4) * t23 * self%P3 - t55 + 0.3D1 * x) * t40) * t6 +  &
        !                                                 0.32D2 / 0.3D1 * self%P5 * ((0.3D1 / 0.16D2 - t14 / 0.2D1) * t1 + t48))    &
        !                                                * self%P2 * t83 + 0.8D1 * self%P1 * t17 * (t6 + 0.1D1) * (-0.1D1 / 0.4D1    &
        !                                                + t87 * t32 - 0.2D1 / 0.3D1 * t87 * t6 + 0.3D1 / 0.4D1 * t7) * (t6 - 0.1D1)

        !> MAPLE output
        !> extra variables
        real(dl) :: t1
        real(dl) :: t2
        real(dl) :: t3
        real(dl) :: t4
        real(dl) :: t6
        real(dl) :: t9
        real(dl) :: t10
        real(dl) :: t13
        real(dl) :: t15
        real(dl) :: t16
        real(dl) :: t17
        real(dl) :: t21
        real(dl) :: t22
        real(dl) :: t27
        real(dl) :: t30
        real(dl) :: t31
        real(dl) :: t33
        real(dl) :: t34
        real(dl) :: t38
        real(dl) :: t39
        real(dl) :: t40
        real(dl) :: t41
        real(dl) :: t44
        real(dl) :: t47
        real(dl) :: t49
        real(dl) :: t50
        real(dl) :: t52
        real(dl) :: t56
        real(dl) :: t58
        real(dl) :: t59
        real(dl) :: t61
        real(dl) :: t65
        real(dl) :: t73
        real(dl) :: t74
        real(dl) :: t75
        real(dl) :: t81
        real(dl) :: t84
        real(dl) :: t87
        real(dl) :: t102
        real(dl) :: t104
        real(dl) :: t114
        real(dl) :: t121
        real(dl) :: t122
        real(dl) :: t130
        real(dl) :: t142
        real(dl) :: t145
        real(dl) :: t147

        !> defining the extra variables
        t1 = self%P3 ** 2
        t2 = self%P2 * t1
        t3 = x - dble(self%P4)
        t4 = t3 ** 2
        t6 = exp(-self%P3 * t4)
        t9 = tanh(self%P5 * t3)
        t10 = t3 * t6 * t9
        t13 = self%P2 * self%P3
        t15 = t9 ** 2
        t16 = 0.1D1 - t15
        t17 = t6 * self%P5 * t16
        t21 = self%P2 * t1 * self%P3
        t22 = t4 * t3
        t27 = t2 * t4
        t30 = t13 * t3
        t31 = self%P5 ** 2
        t33 = t9 * t16
        t34 = t6 * t31 * t33
        t38 = self%P2 * t6 + self%P1
        t39 = t31 * self%P5
        t40 = t38 * t39
        t41 = t16 ** 2
        t44 = t15 * t16
        t47 = -0.8D1 * t21 * t22 * t6 * t9 + 0.12D2 * t2 * t10 - 0.6D1 * t13 * t17 + 0.12D2 * t27 * t17 + 0.12D2 * t30 * t34 - 0.2D1 * t40 * t41 + 0.4D1 * t40 * t44
        t49 = 1 - self%P4
        t50 = t49 ** 2
        t52 = exp(-self%P3 * dble(t50))
        t56 = tanh(self%P5 * dble(t49))
        t58 = t38 * t9 + self%omegaL - (self%P2 * t52 + self%P1) * t56
        t59 = 0.1D1 / t58
        t61 = t6 * t9
        t65 = t4 * t6 * t9
        t73 = -0.2D1 * t38 * t31 * t33 - 0.2D1 * t13 * t61 - 0.4D1 * t30 * t17 + 0.4D1 * t2 * t65
        t74 = t58 ** 2
        t75 = 0.1D1 / t74
        t81 = t38 * self%P5 * t16 - 0.2D1 * t13 * t10
        t84 = t81 ** 2
        t87 = 0.1D1 / t74 / t58
        t102 = t1 ** 2
        t104 = t4 ** 2
        t114 = t6 * t39
        t121 = t31 ** 2
        t122 = t38 * t121
        t130 = 0.16D2 * self%P2 * t102 * t104 * t6 * t9 + 0.24D2 * t13 * t6 *t31 * t9 * t16 - 0.8D1 * t122 * t15 * t9 * t16 + 0.16D2 * t30 * t114 * t41 &
                - 0.32D2 * t30 * t114 * t44 + 0.16D2 * t122 * t41 * t9 +0.48D2 * t2 * t3 * t17 - 0.32D2 * t21 * t22 * t17 + 0.12D2 * t2 *t61 - 0.48D2 * t21 * t65 - 0.48D2 * t27 * t34
        t142 = t73 ** 2
        t145 = t84 ** 2
        t147 = t74 ** 2

        ReconstructedFitParametrized1DThirdDerivative = -t47 * t59 + 0.3D1 * t73 * t75 * t81 - 0.2D1 * t84 * t81 * t87 - x * t130 * t59 / 0.3D1 + 0.4D1 / 0.3D1 &
                                                        * x * t47* t75 * t81 - 0.4D1 * x * t73 * t87 * t84 + x * t142 * t75 + 0.2D1 * x * t145 / t147


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
        real(dl) :: t1
        real(dl) :: t10
        real(dl) :: t12
        real(dl) :: t14
        real(dl) :: t2
        real(dl) :: t20
        real(dl) :: t4
        real(dl) :: t6
        real(dl) :: t9

        !> defining the extra variables
        t1 = -x + self%P4
        t2 = t1 ** 2
        t4 = exp(-self%P3 * t2)
        t6 = tanh(self%P5 * t1)
        t9 = -1.d0 + self%P4
        t10 = t9 ** 2
        t12 = exp(-self%P3 * t10)
        t14 = tanh(self%P5 * t9)
        t20 = x ** 2

        !> Integral
        ReconstructedFitParametrized1DIntegral = (self%P2 * t12 * t14 - self%P2 * t4 * t6 + self%P1 * t14 - self%P1 *   &
                                                    t6 + self%omegaL) * t20 / self%omegaL

    end function ReconstructedFitParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_reconstructed_fit_parametrizations_1D

!----------------------------------------------------------------------------------------
