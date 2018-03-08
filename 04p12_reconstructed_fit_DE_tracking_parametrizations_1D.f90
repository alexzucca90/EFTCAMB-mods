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


module EFTCAMB_reconstructed_DE_fit_tracking_parametrizations_1D

    use precision
    use AMLutils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public reconstructed_DE_fit_tracking_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: reconstructed_DE_fit_tracking_parametrization_1D

        real(dl) :: P1
        real(dl) :: P2
        real(dl) :: P3
        real(dl) :: P4
        real(dl) :: omegaL

    contains

        ! utility functions:
        procedure :: set_param_number      => ReconstructedDEFitTrackParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Taylor expansion parametrized function.
        procedure :: init_parameters       => ReconstructedDEFitTrackParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => ReconstructedDEFitTrackParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => ReconstructedDEFitTrackParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.

        ! evaluation procedures:
        procedure :: value                 => ReconstructedDEFitTrackParametrized1DValue               !< function that returns the value of the Taylor expansion.
        procedure :: first_derivative      => ReconstructedDEFitTrackParametrized1DFirstDerivative     !< function that returns the first derivative of the Taylor expansion.
        procedure :: second_derivative     => ReconstructedDEFitTrackParametrized1DSecondDerivative    !< function that returns the second derivative of the Taylor expansion.
        procedure :: third_derivative      => ReconstructedDEFitTrackParametrized1DThirdDerivative     !< function that returns the third derivative of the Taylor expansion.
        procedure :: integral              => ReconstructedDEFitTrackParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type reconstructed_DE_fit_tracking_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Taylor expansion parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Taylor expansion parametrized function.
    subroutine ReconstructedDEFitTrackParametrized1DSetParamNumber( self )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 5

    end subroutine ReconstructedDEFitTrackParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine ReconstructedDEFitTrackParametrized1DInitParams( self, array )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D)         :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)          :: array  !< input array with the values of the parameters.

        self%P1     = array(1)
        self%P2     = array(2)
        self%P3     = array(3)
        self%P4     = array(4)
        self%omegaL = array(5)

    end subroutine ReconstructedDEFitTrackParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine ReconstructedDEFitTrackParametrized1DParameterValues( self, i, value )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D)       :: self        !< the base class
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
                value = self%omegaL
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine ReconstructedDEFitTrackParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine ReconstructedDEFitTrackParametrized1DFeedback( self, print_params )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D)  :: self         !< the base class
        logical, optional                               :: print_params !< optional flag that decised whether to print numerical values
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

    end subroutine ReconstructedDEFitTrackParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function ReconstructedDEFitTrackParametrized1DValue( self, x, eft_cache )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D)     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedDEFitTrackParametrized1DValue                 !< the output value

        !> Maple output
        !> extra parameters
        real(dl) :: t3
        real(dl) :: t5

        !> defining the extra variables
        t3 = tanh(self%P3 * (x - self%P2))
        t5 = x ** 2

        !> value of the function
        ReconstructedDEFitTrackParametrized1DValue = self%P1 * t3 + self%P1 - self%P2 + self%P4 / t5 / x


    end function ReconstructedDEFitTrackParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function ReconstructedDEFitTrackParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D)                  :: self        !< the base class
        real(dl), intent(in)                                            :: x                    !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional              :: eft_cache            !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedDEFitTrackParametrized1DFirstDerivative                        !< the output value

        !> Maple output
        !> some extra variables for optimization
        real(dl) :: t4
        real(dl) :: t5
        real(dl) :: t8
        real(dl) :: t9

        !> defining the extra variables
        t4 = tanh(self%P3 * (x - self%P2))
        t5 = t4 ** 2
        t8 = x ** 2
        t9 = t8 ** 2

        !> firt derivative
        ReconstructedDEFitTrackParametrized1DFirstDerivative = self%P1 * self%P3 * (0.1D1 - t5) - 0.3D1 * self%P4 / t9


    end function ReconstructedDEFitTrackParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function ReconstructedDEFitTrackParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D)     :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedDEFitTrackParametrized1DSecondDerivative      !< the output value

        !> Maple output
        !> extra varibales
        real(dl) :: t1
        real(dl) :: t5
        real(dl) :: t6
        real(dl) :: t11
        real(dl) :: t12

        !> defining the extra variables
        t1 = self%P3 ** 2
        t5 = tanh(self%P3 * (x - self%P2))
        t6 = t5 ** 2
        t11 = x ** 2
        t12 = t11 ** 2

        !> second derivative
        ReconstructedDEFitTrackParametrized1DSecondDerivative = -0.2D1 * self%P1 * t1 * t5 * (0.1D1 - t6) + 0.12D2 * self%P4 / t12 / x


    end function ReconstructedDEFitTrackParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function ReconstructedDEFitTrackParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D)     :: self      !< the base class
        real(dl), intent(in)                                        :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional          :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedDEFitTrackParametrized1DThirdDerivative       !< the output value

        !> Maple output
        !> extra varibales
        real(dl) :: t1
        real(dl) :: t3
        real(dl) :: t6
        real(dl) :: t7
        real(dl) :: t8
        real(dl) :: t9
        real(dl) :: t15
        real(dl) :: t16

        !> defining the extra variables
        t1 = self%P3 ** 2
        t3 = self%P1 * t1 * self%P3
        t6 = tanh(self%P3 * (x - self%P2))
        t7 = t6 ** 2
        t8 = 0.1D1 - t7
        t9 = t8 ** 2
        t15 = x ** 2
        t16 = t15 ** 2

        !> third derivative
        ReconstructedDEFitTrackParametrized1DThirdDerivative = -0.2D1 * t3 * t9 + 0.4D1 * t3 * t7 * t8 - 0.60D2 * self%P4 / t16 / t15


    end function ReconstructedDEFitTrackParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function ReconstructedDEFitTrackParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(reconstructed_DE_fit_tracking_parametrization_1D)     :: self      !< the base class
        real(dl), intent(in)                                        :: x         !< the input scale factor
        type(EFTCAMB_timestep_cache), intent(in), optional          :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: ReconstructedDEFitTrackParametrized1DIntegral                !< the output value

        write(*,*) "WARNING: the function Reconstructed_DE_tracking_fit%Integral is not supposed to be called. Returning 0"
        ReconstructedDEFitTrackParametrized1DIntegral = 0._dl

    end function ReconstructedDEFitTrackParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_reconstructed_DE_fit_tracking_parametrizations_1D

!----------------------------------------------------------------------------------------
