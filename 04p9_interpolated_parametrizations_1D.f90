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

!> @file 04p9_interpolated_parametrizations_1D.f90
!! This file contains the interpolation of a function of the scale factor.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the interpolation parametrization,
!! inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone
!> @author Alex Zucca

module EFTCAMB_interpolated_function_1D

    use precision
    use AMLutils
    use GBD_Utils
    use IniFile
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    integer, parameter :: number_of_points = 100

    private

    public interpolated_function_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Taylor expansion parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: interpolated_function_1D

        ! Interpolation function variables
        character(len=50)                   :: file, file_name
        !integer                             :: number_of_points
        logical                             :: is_function_loaded = .false.

        ! the arrays for the interpolations.
        real(dl) :: InterpolationArray_a(number_of_points),      InterpolationArray_value(number_of_points), &
                    InterpolationArray_deriv1(number_of_points), InterpolationArray_deriv2(number_of_points), &
                    InterpolationArray_deriv3(number_of_points)

        real(dl) :: fake_param !< parameter inserted just for compatibility purposes.

    contains

        ! utility functions
        procedure :: set_param_number    => InterpolatedFunction1DSetParamNumber   !< subroutine that sets the number of parameters of the interpolated function.
        procedure :: init_parameters     => InterpolatedFunction1DInitParams       !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value     => InterpolatedFunction1DParameterValues  !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback            => InterpolatedFunction1DFeedback         !< subroutine that prints to screen the informations about the function.
        procedure :: init_from_file      => InterpolatedFunction1DInitFromFile     !< subroutine that initializes a few parameters from files
        procedure :: set_param_names     => InterpolatedFunction1DSetParamNames    !< subroutine that sets the file_name

        ! Initialization here
        procedure :: initialize_function => InterpolatedFunction1DInitialization   !< subroutine that initializes the interpolation

        ! evaluation procedures
        procedure :: value               => InterpolatedFunction1DValue            !< function that returns the value of the interpolation.
        procedure :: first_derivative    => InterpolatedFunction1DFirstDerivative  !< function that returns the first derivative of the interpolation.
        procedure :: second_derivative   => InterpolatedFunction1DSecondDerivative !< function that returns the second derivative of the interpolation.
        procedure :: third_derivative    => InterpolatedFunction1DThirdDerivative  !< function that returns the third derivative of the interpolation.
        procedure :: integral            => InterpolatedFunction1DIntegral         !< function that returns the strange integral that we need for w_DE.

        ! I think I should replace some stuff here... let's see

    end type interpolated_function_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the interpolated function
    ! ---------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
!> Subroutine the sets the number of parameter of the interpolation of the function
subroutine InterpolatedFunction1DSetParamNumber(self)
    implicit none
    class(interpolated_function_1D) :: self

    self%parameter_number = 0 ! Not sure about this.. I don't know what's their purpose then

end subroutine InterpolatedFunction1DSetParamNumber


! ---------------------------------------------------------------------------------------------
!> Subroutine that initializes the function parameters based on the values found in an input array.
subroutine InterpolatedFunction1DInitParams( self, array )

    implicit none

    class(interpolated_function_1D) :: self
    real(dl), dimension(self%parameter_number), intent(in) :: array

    ! That's it, nothing to do..
end subroutine InterpolatedFunction1DInitParams



! ---------------------------------------------------------------------------------------------
!> Subroutine that returns the value of the function i-th parameter.
subroutine InterpolatedFunction1DParameterValues( self, i, value )

    implicit none

    class(interpolated_function_1D) :: self
    integer,intent(in)              :: i
    real(dl), intent(out)           :: value

    value = 0.d0

end subroutine InterpolatedFunction1DParameterValues


! ---------------------------------------------------------------------------------------------
!> Subroutine that prints to screen the informations about the function.
subroutine InterpolatedFunction1DFeedback( self, print_params )

    implicit none

    class(interpolated_function_1D) :: self
    logical, optional               :: print_params
    logical                         :: print_params_temp

    if ( present(print_params) ) then
        print_params_temp = print_params
    else
        print_params_temp = .True.
    end if

    write(*,*)     'Interpolated Function: ', self%name
    if ( print_params_temp ) then
        write(*,*) self%file_name, '=', self%file
    end if

end subroutine InterpolatedFunction1DFeedback


! ---------------------------------------------------------------------------------------------
!> Subroutine that reads a Ini file looking for name of the file containing the function.
subroutine InterpolatedFunction1DInitFromFile( self, Ini )

    implicit none

    class(interpolated_function_1D)     :: self   !< the base class
    type(TIniFile)                      :: Ini    !< Input ini file
    character(len=50)                   :: file_name, file

    file_name = self%file_name
    file = Ini_Read_String_File(Ini, TRIM(file_name))

    self%file = file

    !initialize the function parameters from the vector:
    !call self%init_parameters( parameters )

end subroutine InterpolatedFunction1DInitFromFile


! ---------------------------------------------------------------------------------------------
!> subroutine that sets the name of the file that we want to load
subroutine InterpolatedFunction1DSetParamNames(self, param_names, param_names_latex)

    implicit none

    class(interpolated_function_1D)                         :: self            !< the base class
    character(len=50), intent(in), dimension(:)             :: param_names       !< the name of the file
    character(len=50), intent(in), dimension(:), optional   :: param_names_latex !< the name of the file in latex (useless)

    self%file_name = param_names(1)

end subroutine InterpolatedFunction1DSetParamNames




! ---------------------------------------------------------------------------------------------
!> subroutine that initializes the spline arrays (2nd derivative arrays) of the function
subroutine InterpolatedFunction1DInitialization(self)

    implicit none

    class(interpolated_function_1D) :: self
    integer                         :: i

    ! Open the file
    open(unit=7, file = self%file, status = "unknown")

    do i  = 1, number_of_points
        ! read a, f(a) from file and fill the array
        read(7,*) self%InterpolationArray_a(i), self%InterpolationArray_value(i), self%InterpolationArray_deriv1(i)
    end do

    ! close the file
    close(7)

    ! initialize the interpolations
    call GBD_spline_double(self%InterpolationArray_a, self%InterpolationArray_value,  number_of_points, self%InterpolationArray_deriv2)
    call GBD_spline_double(self%InterpolationArray_a, self%InterpolationArray_deriv1, number_of_points, self%InterpolationArray_deriv3)

    ! done
    self%is_function_loaded = .true.

end subroutine InterpolatedFunction1DInitialization


! ---------------------------------------------------------------------------------------------
!> Function that returns the value of the function in the scale factor.
function InterpolatedFunction1DValue( self, x, eft_cache )

    implicit none

    class(interpolated_function_1D)                    :: self        !< the base class
    real(dl), intent(in)                               :: x           !< the input scale factor
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache   !< the optional input EFTCAMB cache
    real(dl) :: InterpolatedFunction1DValue                           !< the output value

    if (.not. self%is_function_loaded) call self%initialize_function

    InterpolatedFunction1DValue = GBD_spline_val(x, self%InterpolationArray_a, self%InterpolationArray_value,&
                                                    self%InterpolationArray_deriv2, number_of_points)

end function InterpolatedFunction1DValue


! ---------------------------------------------------------------------------------------------
!> Function that returns the value of the first derivative, wrt scale factor, of the function.
function InterpolatedFunction1DFirstDerivative(self, x, eft_cache )

    implicit none

    class(interpolated_function_1D)                    :: self
    real(dl), intent(in)                               :: x
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache
    real(dl) :: InterpolatedFunction1DFirstDerivative

    if (.not. self%is_function_loaded) call self%initialize_function

    InterpolatedFunction1DFirstDerivative = GBD_spline_val(x, self%InterpolationArray_a, self%InterpolationArray_deriv1,&
                                                              self%InterpolationArray_deriv3, number_of_points)

end function InterpolatedFunction1DFirstDerivative


! ---------------------------------------------------------------------------------------------
!> Function that returns the second derivative of the function.
function InterpolatedFunction1DSecondDerivative( self, x, eft_cache )

    implicit none

    class(interpolated_function_1D)                    :: self
    real(dl), intent(in)                               :: x
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache
    real(dl) :: InterpolatedFunction1DSecondDerivative

    if (.not. self%is_function_loaded) call self%initialize_function

    InterpolatedFunction1DSecondDerivative = GBD_spline_1der(x, self%InterpolationArray_a,  self%InterpolationArray_deriv1,&
                                                                self%InterpolationArray_deriv3, number_of_points)

end function InterpolatedFunction1DSecondDerivative



! ---------------------------------------------------------------------------------------------
!> Function that returns the second derivative of the function.
function InterpolatedFunction1DThirdDerivative( self, x, eft_cache )

    implicit none

    class(interpolated_function_1D)                    :: self
    real(dl), intent(in)                               :: x
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache
    real(dl) :: InterpolatedFunction1DThirdDerivative

    if (.not. self%is_function_loaded) call self%initialize_function

    InterpolatedFunction1DThirdDerivative = GBD_spline_2der(x, self%InterpolationArray_a, self%InterpolationArray_deriv1,&
                                                            self%InterpolationArray_deriv3, number_of_points)

end function InterpolatedFunction1DThirdDerivative


function InterpolatedFunction1DIntegral( self, x, eft_cache )

    implicit none

    class(interpolated_function_1D) :: self
    real(dl), intent(in) :: x
    type(EFTCAMB_timestep_cache), intent(in), optional :: eft_cache
    real(dl) :: InterpolatedFunction1DIntegral

    ! Now what to do? rombint most likely
    if (.not. self%is_function_loaded) call self%initialize_function



end function InterpolatedFunction1DIntegral


    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_interpolated_function_1D
!----------------------------------------------------------------------------------------
