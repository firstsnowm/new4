program sod_shock_tube
    implicit none
   
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: delta_x = 0.001_dp
    real(dp), parameter :: delta_t = 0.0002_dp
    real(dp), parameter :: length_l = 0.5_dp
    real(dp), parameter :: length_r = 0.5_dp
    integer, parameter :: n_l = int(length_l / delta_x)
    integer, parameter :: n_r = int(length_r / delta_x)
    integer, parameter :: time_step = int(1 / delta_t)
    integer :: i, j, m

    real(dp), dimension(n_l + n_r) :: p, psi, u, mass, e, c
    real(dp), dimension(3, n_l + n_r) :: F_plus, F_minus
    real(dp), dimension(3, n_l + n_r) :: F_plus_minus_1, F_minus_minus_1
    real(dp), dimension(3, n_l + n_r) :: F_plus_plus_1, F_minus_plus_1
    real(dp), dimension(3, n_l + n_r) :: F_plus_minus_3, F_minus_minus_3
    real(dp), dimension(3, n_l + n_r) :: F_plus_plus_3, F_minus_plus_3
    REAL(dp), dimension(3, n_l + n_r) :: x1, x2, x3, x4

    real(dp), dimension(n_l + n_r) :: lamda1, lamda2, lamda3
    real(dp), parameter :: gama = 1.4_dp

    ! Initialize conditions
    p(1:n_l) = 1.0_dp
    p(n_l+1:n_l+n_r) = 0.1_dp
    psi(1:n_l) = 1.0_dp
    psi(n_l+1:n_l+n_r) = 0.125_dp
    u = 0.0_dp
    mass = psi * u
    e = p / (gama - 1.0_dp) + 0.5_dp * psi * u**2

    ! Time iteration
    
    do i = 1, 500
        c = sqrt(gama * p / psi)
        lamda1 = u
        lamda2 = u + c
        lamda3 = u - c

        ! Calculate fluxes using Steger-Warming flux vector splitting
        F_plus(1,:) = psi / (2.0_dp * gama) * (2.0_dp * (gama - 1.0_dp) &
        * max(lamda1, 0.0_dp) + max(lamda2, 0.0_dp) + max(lamda3, 0.0_dp))
        F_plus(2,:) = psi / (2.0_dp * gama) * (2.0_dp * (gama - 1.0_dp) &
        * max(lamda1, 0.0_dp) * u + max(lamda2, 0.0_dp) * (u + c) + max(lamda3, 0.0_dp) * (u - c))
        F_plus(3,:) = psi / (2.0_dp * gama) * &
        ((gama - 1.0_dp) * max(lamda1, 0.0_dp) * u**2 + &
        (3.0_dp - gama) / (2.0_dp * (gama - 1.0_dp)) * &
        (max(lamda2, 0.0_dp) + max(lamda3, 0.0_dp)) * c**2 + &
        0.5_dp * max(lamda2, 0.0_dp) * (u + c)**2 + &
        0.5_dp * max(lamda3, 0.0_dp) * (u - c)**2)


        F_minus(1,:) = psi / (2.0_dp * gama) * (2.0_dp * (gama - 1.0_dp)&
         * min(lamda1, 0.0_dp) + min(lamda2, 0.0_dp) + min(lamda3, 0.0_dp))
        F_minus(2,:) = psi / (2.0_dp * gama) * (2.0_dp * (gama - 1.0_dp) &
        * min(lamda1, 0.0_dp) * u + min(lamda2, 0.0_dp) * (u + c) + min(lamda3, 0.0_dp) * (u - c))
        F_minus(3,:) = psi / (2.0_dp * gama) * &
               ((gama - 1.0_dp) * min(lamda1, 0.0_dp) * u**2 + &
               (3.0_dp - gama) / (2.0_dp * (gama - 1.0_dp)) * &
               (min(lamda2, 0.0_dp) + min(lamda3, 0.0_dp)) * c**2 + &
               0.5_dp * min(lamda2, 0.0_dp) * (u + c)**2 + &
               0.5_dp * min(lamda3, 0.0_dp) * (u - c)**2)


        ! Calculate flux differences
        do m = 3, n_l + n_r - 2
            do j = 1, 3
                F_plus_minus_1(j,m) = F_plus(j,m) - F_plus(j,m-1)
                F_minus_minus_1(j,m) = F_minus(j,m) - F_minus(j,m-1)

                F_plus_plus_1(j,m) = F_plus(j,m+1) - F_plus(j,m)
                F_minus_plus_1(j,m) = F_minus(j,m+1) - F_minus(j,m)

                F_plus_minus_3(j,m) = F_plus(j,m-1) - F_plus(j,m-2)
                F_minus_minus_3(j,m) = F_minus(j,m-1) - F_minus(j,m-2)

                F_plus_plus_3(j,m) = F_plus(j,m+2) - F_plus(j,m+1)
                F_minus_plus_3(j,m) = F_minus(j,m+2) - F_minus(j,m+1)
            end do
        end do

        ! Calculate min mod coefficients
        do m = 3, n_l + n_r - 2
            do j = 1, 3
                x1(j,m) = REAL(logical_to_int(signum(F_plus_minus_1(j,m)) == signum(F_plus_plus_1(j,m))), dp)
x2(j,m) = REAL(logical_to_int(signum(F_minus_plus_1(j,m)) == signum(F_minus_plus_3(j,m))), dp)
x3(j,m) = REAL(logical_to_int(signum(F_plus_minus_3(j,m)) == signum(F_plus_minus_1(j,m))), dp)
x4(j,m) = REAL(logical_to_int(signum(F_minus_minus_1(j,m)) == signum(F_minus_plus_1(j,m))), dp)
            end do
        end do

        ! Update conserved quantities
        do m = 3, n_l + n_r - 2
            psi(m) = psi(m) - delta_t / delta_x * (F_plus(1,m) + 0.5_dp * x1(1,m) * signum(F_plus_minus_1(1,m)) &
                 * min(abs(F_plus_minus_1(1,m)), abs(F_plus_plus_1(1,m))) &
                 + F_minus(1,m+1) - 0.5_dp * x2(1,m) * signum(F_minus_plus_1(1,m)) &
                 * min(abs(F_minus_plus_1(1,m)), abs(F_minus_plus_3(1,m))) &
                 - F_plus(1,m-1) - 0.5_dp * x3(1,m) * signum(F_plus_minus_3(1,m))&
                  * min(abs(F_plus_minus_3(1,m)), abs(F_plus_minus_1(1,m))) &
                 - F_minus(1,m) + 0.5_dp * x4(1,m) * signum(F_minus_minus_1(1,m)) &
                 * min(abs(F_minus_minus_1(1,m)), abs(F_minus_plus_1(1,m))))

            mass(m) = mass(m) - delta_t / delta_x * (F_plus(2,m) + 0.5_dp * x1(2,m) * signum(F_plus_minus_1(2,m))&
                  * min(abs(F_plus_minus_1(2,m)), abs(F_plus_plus_1(2,m))) &
                 + F_minus(2,m+1) - 0.5_dp * x2(2,m) * signum(F_minus_plus_1(2,m))&
                  * min(abs(F_minus_plus_1(2,m)), abs(F_minus_plus_3(2,m))) &
                 - F_plus(2,m-1) - 0.5_dp * x3(2,m) * signum(F_plus_minus_3(2,m)) &
                 * min(abs(F_plus_minus_3(2,m)), abs(F_plus_minus_1(2,m))) &
                 - F_minus(2,m) + 0.5_dp * x4(2,m) * signum(F_minus_minus_1(2,m))&
                  * min(abs(F_minus_minus_1(2,m)), abs(F_minus_plus_1(2,m))))
                    u(m) = mass(m) / psi(m)

        e(m) = e(m) - delta_t / delta_x * (F_plus(3,m) + 0.5_dp * x1(3,m) * signum(F_plus_minus_1(3,m)) &
             * min(abs(F_plus_minus_1(3,m)), abs(F_plus_plus_1(3,m))) &
             + F_minus(3,m+1) - 0.5_dp * x2(3,m) * signum(F_minus_plus_1(3,m))&
              * min(abs(F_minus_plus_1(3,m)), abs(F_minus_plus_3(3,m))) &
             - F_plus(3,m-1) - 0.5_dp * x3(3,m) * signum(F_plus_minus_3(3,m)) &
             * min(abs(F_plus_minus_3(3,m)), abs(F_plus_minus_1(3,m))) &
             - F_minus(3,m) + 0.5_dp * x4(3,m) * signum(F_minus_minus_1(3,m)) &
             * min(abs(F_minus_minus_1(3,m)), abs(F_minus_plus_1(3,m))))

        p(m) = (gama - 1.0_dp) * (e(m) - 0.5_dp * psi(m) * u(m)**2)
        end do
! 输出当前时间步长的信息
        open (11,file='3.plt')
write(*, '(A,F6.4,A)') 'Time step: t = ', delta_t * i, ' s'
do m = 3, n_l + n_r - 2
! 输出压力、密度、速度和能量等信息
write(11, 10)delta_x*m ,p(m), psi(m), u(m), e(m)
10 FORMAT (1x,f20.14,2x,f20.14,3x,f20.14,4x,f20.14,5x,f20.14,6x,f20.14)
enddo
close(11)

! 可以根据需要添加其他变量的输出语句

    !! Plot at each time step (optional)
    !write(*, '(A,F6.4,A)') 'Time step: t = ', delta_t * i, ' s'
    !! Plotting code here if needed
    end do
    contains
    function signum(x)
    real(dp), intent(in) :: x
    real(dp) :: signum

    if (x > 0.0_dp) then
        signum = 1.0_dp
    else if (x < 0.0_dp) then
        signum = -1.0_dp
    else
        signum = 0.0_dp
    end if
end function signum
INTEGER FUNCTION logical_to_int(logical_value)
    LOGICAL, INTENT(IN) :: logical_value
    IF (logical_value) THEN
        logical_to_int = 1
    ELSE
        logical_to_int = 0
    END IF
END FUNCTION logical_to_int

end program sod_shock_tube
