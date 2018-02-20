subroutine gl10(n,fcn,y, dt)
use Precision
use AMLUtils
implicit none
!---------------------------------------------------------------
! Gauss-Legendre 10th order
! written by Alex Zucca: azucca@sfu.ca
! paramters:
!            n = number of equations
!            fnc = external subroutines that computes the derivatives
!                   fcn has the form
!
!                   subroutine fcn(n,x,y,yprime)
!                       integer :: n
!                       real    :: y(n), yprime(n)
!
!                          ***************
!
!                   end subroutine fcn
!
!            y = array of variables, y(1) is the independent variable
!            dt = interval of independent variable
!
!   if you use this subroutine for your projects, please remove
!
!          use AMLUtils
!          use Precision
!
!      and uncomment the line
!          integer,parameter :: dl=kind(1.d0)
!
!
!
!---------------------------------------------------------------
    integer, parameter :: s = 5
    integer :: n

    !> uncomment the following line if you use this subroutine for your projects
    !integer,parameter :: dl=kind(1.d0)

    real(dl) :: y(n), dt
    real(dl) :: g(n,s), dx
    integer :: i, k
    external fcn

!> Butcher tableau for 8th order Gauss-Legendre method
real(dl), parameter :: a(s,s) = reshape((/ &
0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
0.5923172126404727187856601017997934066Q-1 /),(/s,s/))
real, parameter ::   b(s) = (/ &
1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
1.1846344252809454375713202035995868132Q-1 /)

    !> some debug info
    !write(*,*) " y = ", y

    !> iterate trial steps
    g=0.0
    do k = 1,16
        g = matmul(g,a)
        do i = 1,s
            !> Debug info
            !write(*,*) "y+g(:,i)*dt = ", y + g(:,i)*dt

            call fcn(n, y + g(:,i)*dt, g(:,i))
        end do
        !> Debug info
        !write(*,*) g(1,:)
        !pause
    end do

    !> update the solution
    y = y + matmul(g,b)*dt

    !> Debug info
    !write(*,*) " y = ", y
    !write(*,*)
    !pause
end subroutine gl10

