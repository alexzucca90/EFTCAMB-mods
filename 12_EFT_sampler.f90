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

!> @file 12_EFT_sampler.f90
!! This file contains the EFTCAMB sampler parameters


!----------------------------------------------------------------------------------------
!> This module contains the parameters for the sampling paramters

!> @author Alex Zucca

module EFT_sampler

    use Precision
    use IniFile
    use AMLutils
    use EFTCAMB_cache
    use EFT_def
    use random

    implicit none

    Type EFTSamplingParameters

        !> xDE MOD START: adding the parameters for the sampling of xDE
        !> xDE sampling variables
        integer :: xDE_n_bins
        real(dl), dimension(:), allocatable     :: xDE_scale_factor_means
        real(dl), dimension(:), allocatable     :: xDE_means
        real(dl), dimension(:), allocatable     :: xDE_sample
        real(dl), dimension(:,:), allocatable   :: xDE_covmat
        real(dl), dimension(:), allocatable     :: xDE_covmat_lin
        real(dl), dimension(:), allocatable     :: xDE_bare_sample
        real(dl), dimension(:), allocatable     :: xDE_covmat_low_tri

        !> Input file for sampling xDE
        character(LEN=Ini_max_string_len)       :: xDE_means_filename
        character(LEN=Ini_max_string_len)       :: xDE_covmat_filename
        logical                                 :: first = .true.
        !> xDE MOD END

    contains

        procedure :: initialize_xDE_distribution    => EFTSamplingParameters_InitlizexDEDistribution
        procedure :: draw_sample_xDE                => EFTSamplingParameters_DrawSamplexDE



    end type EFTSamplingParameters



contains
!-----------------------------------------------------------
!> This subroutine draw a sample for the DE density xDE(a)
subroutine EFTSamplingParameters_InitlizexDEDistribution( self, Ini )

    implicit none

    class(EFTSamplingParameters)    :: self     !< the base class
    type(TIniFile)                  :: Ini      !< Input ini file
    integer :: i,j

    !> Read the parameters from file
    self%xDE_means_filename     = Ini_Read_String_File( Ini, 'xDE_means_filename' )
    self%xDE_covmat_filename    = Ini_Read_String_File( Ini, 'xDE_covmat_filename' )
    self%xDE_n_bins             = Ini_Read_Int_File( Ini, 'xDE_n_bins', 40 )

    !> Allocate the xDE distributions array
    allocate(self%xDE_scale_factor_means(self%xDE_n_bins))
    allocate(self%xDE_means(self%xDE_n_bins-1))
    allocate(self%xDE_covmat(self%xDE_n_bins-1,self%xDE_n_bins-1))
    !allocate(Array_a(xDE_n_bins), Array_xDE(xDE_n_bins))
    !allocate(Array_xDE1(xDE_n_bins), Array_xDE2(xDE_n_bins))
    allocate(self%xDE_covmat_low_tri((self%xDE_n_bins-1)*((self%xDE_n_bins)/2)))

    !> read the arrays
    open(unit = 1, file = self%xDE_means_filename,  status = 'unknown')
    open(unit = 2, file = self%xDE_covmat_filename, status = 'unknown')

    write(*,*) "    Reading the covariance matrix"
    do i = 1, self%xDE_n_bins - 1
        !write(*,*) i
        read(1,*) self%xDE_means(i)
        read(2,*) self%xDE_covmat(i,:)
    end do

    close(1);close(2)

    !> building the covariance matrix for the random subroutine
    allocate(self%xDE_covmat_lin((self%xDE_n_bins - 1) *(self%xDE_n_bins)/2))
    allocate(self%xDE_sample(self%xDE_n_bins))

    write(*,*) "    Generating the reshhaped covariance matrix"
    do i = 1, self%xDE_n_bins-1
        do j = 1, self%xDE_n_bins-1
            if(j .ge. i) then
                self%xDE_covmat_lin(j*(j-1)/2+i) = self%xDE_covmat(i,j)
            end if
        end do
    end do

    !> here I need to generate the scale factor vector
    do i = 0,self%xDE_n_bins-1
        self%xDE_scale_factor_means(i+1) = 1.d-3 + i*(1.d0 - 1.d-3)/(self%xDE_n_bins-1)
    end do

end subroutine EFTSamplingParameters_InitlizexDEDistribution


!-----------------------------------------------------------
!> This subroutine draw a sample for the DE density xDE(a)
subroutine EFTSamplingParameters_DrawSamplexDE( self )

    implicit none

    class(EFTSamplingParameters)    :: self     !< the base class
    integer :: ier
    integer :: i
    real(dl), dimension(:), allocatable :: xDE_bare_sample

    allocate(xDE_bare_sample(self%xDE_n_bins - 1))

    !> generating the xDE from the distribution
    call  random_mvnorm(self%xDE_n_bins-1, self%xDE_means, self%xDE_covmat_lin, self%xDE_covmat_low_tri, self%first, xDE_bare_sample, ier)

    !> filling the array for the interpolations
    self%xDE_sample(self%xDE_n_bins) = 1.d0
    do i = 1, self%xDE_n_bins-1
        self%xDE_sample(i) = xDE_bare_sample(self%xDE_n_bins - i)
        !write(*,*) i, self%xDE_sample(i)
    end do

    self%first = .false.

end subroutine EFTSamplingParameters_DrawSamplexDE

end module EFT_sampler

!----------------------------------------------------------------------------------------
