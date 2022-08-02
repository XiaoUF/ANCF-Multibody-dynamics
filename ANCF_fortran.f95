! The purpose of this project is to simulate multiple flexible beams
! with ANCF beam model
program flexible_beam

    implicit none

    ! Define termination time, time marching and output intervel
    double precision, parameter :: tend = 10.0d1  ! End time of simulation
    double precision, parameter :: h = 1.0d-5     ! Simulation time step
    double precision, parameter :: hout = 5.0d-3  ! Time step to write to output file (result folder)
    double precision, parameter :: tlog = 0.1d0   ! Time step to write simulation progress (out.log)

    integer :: i, ios, j
    logical :: dirExist

    ! Input parameters for beam
    double precision, parameter :: E = 8.0d8  ! Young's modulus (Pa)
    double precision, parameter :: rho = 8.5d2  ! Density (kg/m^3)
    integer, parameter :: NE = 20   ! Number of ANCF elements per beam
    integer, parameter :: NDOF = 126  ! Number of nodal coordinates, gradient=deficient
    double precision, parameter :: L = 100.0d-2/20.0d0 ! Length of one ANCF element (m)
    double precision, parameter :: thick = 5.0d-4 ! Thickness of blade (m)
    double precision, parameter :: width = 1.0d-2 ! Width of blade (m)
    double precision, parameter :: gx = 0.0d0 ! X direction gravitational acceleration
    double precision, parameter :: gy = 0.0d0 ! Y direction gravitational acceleration
    double precision, parameter :: gz = 1.729d0  ! Z direction gravitational acceleration
    double precision :: A = thick * width ! Cross-section area (m^2)
    double precision :: Inertia = width * thick**3 / 1.2d1  ! Second moment of area (m^4)
    double precision, parameter :: dx0 = 0.10d0 ! Distance between blades
    ! Specify constraints at the ends of beam
    ! 0 : free
    ! 1 : roller constraint
    ! 2 : fixed constraint
    integer, parameter :: botCnstr = 2  ! Fixed constraint at bottom of beam
    integer, parameter :: topCnstr = 0  ! Free constraint at top of beam

    double precision, parameter :: Uv = 0.3d0
    double precision, parameter :: Cv = 0.6d0
    double precision, parameter :: kv = 1.0d0/0.125d0
    double precision, parameter :: zv = 0.5125d0
    double precision, parameter :: tml = 0.125d0

    ! Number of actual beam and output
    integer, parameter :: NBEAM = 5

    ! Define folder name of beams
    character(len=8), dimension(NBEAM) :: beam_folder
    character(len=3) :: char_i

    ! Arries for current status
    double precision, dimension(NDOF, NBEAM) :: disp_cur, vel_cur, acc_cur

    ! Mass matrix
    double precision, dimension(NDOF, NDOF) :: M
    double precision, dimension(NDOF, NDOF) :: M_L, M_U
    double precision, dimension(NDOF, NDOF) :: Mtemp
    integer, dimension(NDOF) :: pivot
    integer :: ok

    ! Gravitational force, total force, elastic force, and damping force
    double precision, dimension(NDOF, NBEAM) :: Qg, Q, Qe, Qp, Ql
    double precision, dimension(2, NBEAM) :: Qdrag, Qam, Qimp

    ! Recorded beam status
    double precision, dimension(NDOF, NBEAM) :: disp_prev, vel_prev, acc_prev

    ! ! Velocity mode for the calculation of damping force
    ! integer, parameter :: nmode = 9   ! Number of mode needs to be smaller than 9
    ! double precision, dimension(NDOF, nmode) :: vel_mode_hor, vel_mode_ver

    ! Handler of input and output files
    integer, parameter :: INUNIT = 77
    integer, parameter :: OUTUNIT = 78
    integer, parameter :: LOGUNIT = 79

    ! Computed time step size
    double precision :: deltaT

    ! Current time and output number
    double precision :: t
    integer :: iout, ilog

    ! CPU time and date
    character(len=8) :: cpu_date
    character(len=10) :: cpu_time
    character(len=1000) :: outchar

    ! Collision counter
    integer, dimension(NBEAM) :: col_count

    ! Elastic and kinetic energy of beam
    double precision, dimension(NBEAM) :: energy_e, energy_t
    double precision, dimension(NDOF) :: vec_tmp
    double precision :: DDOT
    integer :: count

    ! Program start
    open (unit = LOGUNIT, file = 'out.log', iostat = ios, status = 'unknown',&
      action = 'write')
    call date_and_time(cpu_date, cpu_time)
    write (LOGUNIT, *) 'Simulation starts at ',cpu_time(1:2),':',cpu_time(3:4),&
        ':', cpu_time(5:6),' ',cpu_date(1:4),'/',cpu_date(5:6),'/',cpu_date(7:8)
    close (LOGUNIT, status = 'keep')

    ! Read input file
    ! Check if result folder already exists
    inquire(file='result/.', exist = dirExist)
    if (dirExist) then
        ! In case result already exists, read data and continue time marching
        print *, 'Folder exists'
    else
        ! Exam time step size
        call time_step(L, thick, rho, E, deltaT)
        if (h > deltaT) then
            write (LOGUNIT, *) 'Reduce time step size to ', deltaT
        end if

        ! Assign initial status
        t = 0.0d0
        call init_cond(disp_cur, vel_cur, NDOF, NE, L, NBEAM)
        do i = 1, NBEAM
            ! call init_cond(disp_cur(:,i), vel_cur(:,i), NDOF, NE, L)
            call gravitational_force(Qg(:,i), rho, A, L, gx, gy, gz, NDOF, NE)
        end do
        call mass_matrix(M, rho, A, width, L, NDOF, NE, botCnstr, topCnstr)
        Mtemp = M
        M_L(:,:) = 0.0d0
        M_U(:,:) = 0.0d0
        call DGETRF(NDOF, NDOF, Mtemp, NDOF, pivot, ok)
        do i = 1, NDOF
            M_L(i,i) = 1.0d0
            M_U(1:i,i) = Mtemp(1:i,i)
            M_L(i+1:NDOF,i) = Mtemp(i+1:NDOF,i)
        end do

        ! The total force vector and acceleration vector
        call point_load(disp_cur, NDOF, NE, E, dx0, L, thick, NBEAM, Qp, Qimp, col_count)
        do i = 1, NBEAM
            call elastic_force(disp_cur(:,i), NDOF, NE, E, A, Inertia, L, Qe(:,i),&
                energy_e(i))
            call line_load_single(disp_cur(:,i), vel_cur(:,i), NDOF, NE, L, width,&
                dx0, i, Ql(:,i), Qdrag(:,i), Qam(:,i), t, Uv, Cv, kv, zv, tml)
            Q(:,i) = Qg(:,i) + Ql(:,i) + Qp(:,i) - Qe(:,i)
            call add_constraint_load(NDOF, botCnstr, topCnstr, Q(:,i))
            call linear_system(M_L, M_U, Q(:,i), acc_cur(:,i), pivot, NDOF)
        end do

        ! Initialize beam status
        disp_prev(:,:) = 0.0d0
        vel_prev(:,:) = 0.0d0
        acc_prev(:,:) = 0.0d0

        ! Update history
        disp_prev = disp_cur
        vel_prev = vel_cur
        acc_prev = acc_cur

        ! Assign current time and output number
        iout = 1
        ilog = 1
        energy_t(:) = 0.0d0
        col_count(:) = 0

        ! Create folder and save beam parameters
        call system('mkdir result')
        do i = 1, NBEAM
            write(char_i, '(I3.3)') i
            beam_folder(i) = 'beam_'//trim(adjustl(char_i))
            call system('mkdir result/'//beam_folder(i))
            open (unit = OUTUNIT, file = 'result/'//beam_folder(i)//'/parameters.txt',&
                iostat = ios, status = 'new', action = 'write')
            write (OUTUNIT, *) 'ne,n,L,E,A,I,rho,x0,d,t,gx,gy,gz,botCnstr,topCnstr'
            write (OUTUNIT, *) NE,',',NDOF,',',L,',',E,',',A,',',Inertia,',',rho,',',&
                DBLE(i-1)*dx0,',',thick,',',width,',',gx,',',gy,',',gz,',',&
                botCnstr,',',topCnstr
            close (OUTUNIT, status = 'keep')

            ! Create the string for DOFs
            outchar = 't'
            do j = 1, NDOF
                write (char_i, '(I3.1)') j
                outchar = trim(outchar)//',dof_'//trim(adjustl(char_i))
            end do
            ! Create other files
            open (unit = OUTUNIT, file = 'result/'//beam_folder(i)//'/position_result.txt',&
                iostat = ios, status = 'new', action = 'write')
            write (OUTUNIT, *) trim(outchar)
            close (OUTUNIT, status = 'keep')
            open (unit = OUTUNIT, file = 'result/'//beam_folder(i)//'/velocity_result.txt',&
                iostat = ios, status = 'new', action = 'write')
            write (OUTUNIT, *) trim(outchar)
            close (OUTUNIT, status = 'keep')
            open (unit = OUTUNIT, file = 'result/'//beam_folder(i)//'/acceleration_result.txt',&
                iostat = ios, status = 'new', action = 'write')
            write (OUTUNIT, *) trim(outchar)
            close (OUTUNIT, status = 'keep')
            outchar = 't,collision_count,elastic_energy,kinetic_energy'
            open (unit = OUTUNIT, file = 'result/'//beam_folder(i)//'/collision_and_energy.txt',&
                iostat = ios, status = 'new', action = 'write')
            write (OUTUNIT, *) trim(outchar)
            close (OUTUNIT, status = 'keep')
            outchar = 't,drag_x,drag_z'
            open (unit = OUTUNIT, file = 'result/'//beam_folder(i)//'/drag_force.txt',&
                iostat = ios, status = 'new', action = 'write')
            write (OUTUNIT, *) trim(outchar)
            close (OUTUNIT, status = 'keep')
            outchar = 't,fam_x,fam_z'
            open (unit = OUTUNIT, file = 'result/'//beam_folder(i)//'/added_mass_force.txt',&
                iostat = ios, status = 'new', action = 'write')
            write (OUTUNIT, *) trim(outchar)
            close (OUTUNIT, status = 'keep')
            outchar = 't,fimp_x,fimp_z'
            open (unit = OUTUNIT, file = 'result/'//beam_folder(i)//'/impact_force.txt',&
                iostat = ios, status = 'new', action = 'write')
            write (OUTUNIT, *) trim(outchar)
            close (OUTUNIT, status = 'keep')
            call write_output(beam_folder(i), OUTUNIT, t, disp_cur(:,i),&
                vel_cur(:,i), acc_cur(:,i), col_count(i), energy_e(i),&
                energy_t(i), Qdrag(:,i), Qam(:,i), Qimp(:,i), NDOF)
        end do
    end if

    count = 0

    ! Start time marching
    do while (t <= tend)

        ! Determine impact force from current status of beams
        call point_load(disp_prev, NDOF, NE, E, dx0, L, thick, NBEAM, Qp, Qimp, col_count)


        ! Calculate the status for next time step

        do i = 1, NBEAM
            call runge_kutta(h, disp_prev(:,i), vel_prev(:,i), acc_prev(:,i),&
                NDOF, NE, botCnstr, topCnstr, E, A, Inertia, L, M_L, M_U,&
                Qg(:,i), Qp(:,i), disp_cur(:,i), vel_cur(:,i), acc_cur(:,i),&
                Q(:,i), Qdrag(:,i), Qam(:,i), energy_e(i), pivot, width, dx0,&
                i, t, Uv, Cv, kv, zv, tml)
        end do

        ! Store solution in output file
        t = t + h

        if (t >= iout * hout) then
            ! Compute kinetic energy at time t
            energy_t(:) = 0.0d0
            do i = 1, NBEAM
                vec_tmp(:) = 0.0d0
                call DGEMV('n',NDOF,NDOF,1.0d0,M,NDOF,vel_cur(:,i),1,0.0d0,vec_tmp,1)
                energy_t(i) = DDOT(NDOF,vel_cur(:,i),1,vec_tmp,1)
            end do
            ! Save output
            iout = iout + 1
            do i = 1, NBEAM
                call write_output(beam_folder(i), OUTUNIT, t, disp_cur(:,i),&
                    vel_cur(:,i), acc_cur(:,i), col_count(i), energy_e(i),&
                    energy_t(i), Qdrag(:,i), Qam(:,i), Qimp(:,i), NDOF)
            end do
        end if

        ! Update status of beam
        disp_prev = disp_cur
        vel_prev = vel_cur
        acc_prev = acc_cur

        ! Write simulation progress to out.log file
        if (t >= ilog * tlog) then
            ilog = ilog + 1
            call date_and_time(cpu_date, cpu_time)
            open (unit = LOGUNIT, file = 'out.log', iostat = ios, status = 'old',&
                action = 'write', position = 'append')
            write (LOGUNIT, *) 'Simulation proceeds to ', t, ' at ',&
                cpu_time(1:2),':', cpu_time(3:4),':',cpu_time(5:6),' ',&
                cpu_date(1:4),'/',cpu_date(5:6),'/',cpu_date(7:8)
            close (LOGUNIT, status = 'keep')
        end if
    end do

    ! Save result to output folder
    call date_and_time(cpu_date, cpu_time)
    open (unit = LOGUNIT, file = 'out.log', iostat = ios, status = 'old',&
        action = 'write', position = 'append')
    write (LOGUNIT, *) 'Simulation ends at ',cpu_time(1:2),':',cpu_time(3:4),&
        ':',cpu_time(5:6),' ',cpu_date(1:4),'/',cpu_date(5:6),'/',cpu_date(7:8)
    close (LOGUNIT, status = 'keep')

end program flexible_beam

subroutine time_step(L, t, rho, E, deltaT)
    ! This subroutine compute the required time step for a beam
    implicit none

    double precision, intent(in) :: L, t, rho, E
    double precision, intent(out) :: deltaT

    deltaT = DSQRT(3.0d0) * L**2 / (4.73d0)**2 / t * DSQRT(rho / E)

end subroutine time_step

subroutine init_cond(disp_0, vel_0, ndof, ne, L, nbeam)
    ! This subroutine assign initial condition to a beam
    implicit none

    integer, intent(in) :: ndof, ne, nbeam
    double precision, intent(in) :: L
    double precision, dimension(ndof,nbeam), intent(out) :: disp_0, vel_0
    integer :: ie
    logical :: fileExist
    integer :: ibeam

    inquire(file='init_coord', exist = fileExist)
    if (fileExist) then
      open(unit = 101,file='init_coord',action = 'read')
      open(unit = 102,file='init_vel',action = 'read')
      read(101,*) disp_0(:,1)
      read(102,*) vel_0(:,1)
      close(101)
      close(102)
      do ibeam = 2,nbeam
        disp_0(:,ibeam) = disp_0(:,1)
        vel_0(:,ibeam) = vel_0(:,1)
      end do
    else
      disp_0(:,:) = 0.0d0
      vel_0(:,:) = 0.0d0
      do ibeam = 1,nbeam
        do ie = 1, (ne + 1)
            disp_0(6*ie-3,ibeam) = (ie-1) * L
            disp_0(6*ie,ibeam) = 1.0d0
        end do
      end do
    end if

end subroutine init_cond

subroutine mass_matrix(M, rho, A, width, L, ndof, ne, botCnstr, topCnstr)
    ! This subroutine calcualte the generalized mass matrix for a beam
    implicit none

    integer, intent(in) :: ndof, ne, botCnstr, topCnstr
    double precision, intent(in) :: rho, A, L, width
    double precision, dimension(ndof, ndof), intent(out) :: M
    double precision :: m11, m12, m13, m14, m22, m23, m24, m33, m34, m44, factor
    double precision, dimension(12, 12) :: Me
    double precision, dimension(6, 6) :: Me00, Me01, Me10, Me11
    double precision, parameter :: pi = 4.0d0 * DATAN(1.0d0)
    double precision :: rhof = 1000.0d0
    integer :: i, j

    ! Evalute the element level mass matrix (local matrix)
    factor = rho * A * L + rhof * pi / 4.0d0 * width**2 * L
    m11 = factor * 1.3d1 / 3.5d1
    m12 = factor * L * 1.1d1 / 2.1d2
    m13 = factor * 9.0d0 / 7.0d1
    m14 = - factor * L * 1.3d1 / 4.2d2
    m22 = factor * L**2 / 1.05d2
    m23 = - m14
    m24 = - factor * L**2 / 1.4d2
    m33 = m11
    m34 = - m12
    m44 = m22
    Me(:,:) = 0.0d0
    do i = 1, 3
        Me(i,i) = m11
        Me(i+3,i) = m12
        Me(i+6,i) = m13
        Me(i+9,i) = m14
    end do
    do i = 4, 6
        Me(i-3,i) = m12
        Me(i,i) = m22
        Me(i+3,i) = m23
        Me(i+6,i) = m24
    end do
    do i = 7, 9
        Me(i-6,i) = m13
        Me(i-3,i) = m23
        Me(i,i) = m33
        Me(i+3,i) = m34
    end do
    do i = 10, 12
        Me(i-9,i) = m14
        Me(i-6,i) = m24
        Me(i-3,i) = m34
        Me(i,i) = m44
    end do

    ! Assemble to the global mass matrix
    Me00 = Me(1:6, 1:6)
    Me01 = Me(1:6, 7:12)
    Me10 = Me(7:12, 1:6)
    Me11 = Me(7:12, 7:12)
    M(:,:) = 0.0d0
    do i = 1, ne
        j = (i-1) * 6 + 1
        M(j:j+5,j:j+5) = M(j:j+5,j:j+5) + Me00
        M(j+6:j+11,j:j+5) = Me10
        M(j:j+5,j+6:j+11) = Me01
        M(j+6:j+11,j+6:j+11) = M(j+6:j+11,j+6:j+11) + Me11
    end do

    ! Adjust mass matrix according to constraint
    if (botCnstr == 1) then
        M(1:3,:) = 0.0d0
        do i = 1,3
            M(i,i) = 1.0d0
        end do
    else if (botCnstr == 2) then
        M(1:6,:) = 0.0d0
        do i = 1,6
            M(i,i) = 1.0d0
        end do
    end if

    if (topCnstr == 1) then
        M((ndof-5):(ndof-3),:) = 0.0d0
        do i = 1,3
            M(ndof-6+i,ndof-6+i) = 1.0d0
        end do
    else if (topCnstr == 2) then
        M((ndof-5):ndof,:) = 0.0d0
        do i = 1,6
            M(ndof-6+i,ndof-6+i) = 1.0d0
        end do
    end if

end subroutine mass_matrix

subroutine gravitational_force(Q, rho, A, L, gx, gy, gz, ndof, ne)
    ! This subroutine calculate the beam generalized gravitational force
    implicit none

    double precision, intent(in) :: rho, A, L, gx, gy, gz
    integer, intent(in) :: ndof, ne
    double precision, dimension(ndof), intent(out) :: Q
    double precision, dimension(4) :: factor
    double precision, dimension(12) :: Qg
    integer :: i

    ! Evaluate the element gravitational force
    factor(1) = rho * A * L/2.0d0
    factor(2) = rho * A * L**2/1.2d1
    factor(3) = factor(1)
    factor(4) = - factor(2)
    do i = 1, 4
        Qg(3*(i-1)+1) = factor(i) * gx
        Qg(3*(i-1)+2) = factor(i) * gy
        Qg(3*(i-1)+3) = factor(i) * gz
    end do

    ! Assemble generalized gravitational force for a beam with 'ne' elements
    Q(:) = 0.0d0
    do i = 1, ne
        Q((6*i-5):(6*i+6)) = Q((6*i-5):(6*i+6)) + Qg
    end do

end subroutine gravitational_force

! subroutine mode_velocity_vector(nmode, vel_mode, ndof, ne, L, dir)
!     ! Horizontal velocity vector for different mode
!     ! dir = 0 : horizontal velocity vector
!     ! dir = 1 : vertical velocity vector
!     implicit none
!
!     integer, intent(in) :: ndof, ne, nmode, dir
!     double precision, intent(in) :: L
!     double precision, dimension(ndof, nmode), intent(out) :: vel_mode
!     double precision, dimension(ndof) :: vel_mode_cur
!     integer :: i
!     double precision :: DNRM2
!
!     select case (dir)
!         case (0)
!             do i = 1, nmode
!                 call mode_vector(vel_mode_cur, ndof, ne, L, i, 1)
!                 vel_mode_cur = vel_mode_cur / DNRM2(ndof,vel_mode_cur,1)
!                 vel_mode(:, i) = vel_mode_cur
!             end do
!         case (1)
!             do i = 1, nmode
!                 call mode_vector_ver(vel_mode_cur, ndof, ne, L, i, 1)
!                 vel_mode_cur = vel_mode_cur / DNRM2(ndof,vel_mode_cur,1)
!                 vel_mode(:, i) = vel_mode_cur
!             end do
!     end select
!
! end subroutine mode_velocity_vector
!
! subroutine mode_vector(v_mode, ndof, ne, L, mode, term)
!     ! Generate the beam mode vector based on the number of element
!     ! term = 0 : return position vector
!     ! term = 1 : return velocity vector
!     implicit none
!
!     integer, intent(in) :: ndof, ne, mode, term
!     double precision, dimension(ndof), intent(out) :: v_mode
!     double precision, intent(in) :: L
!     double precision :: L_beam
!     integer :: i
!
!     v_mode(:) = 0.0d0
!     L_beam = ne * L
!
!     do i = 1, ne
!         call mode_shape(v_mode(6*i+1), mode, L_beam, L*i, 0)
!         call mode_shape(v_mode(6*i+4), mode, L_beam, L*i, 1)
!     end do
!
!     select case (term)
!         case (0)
!             v_mode(6) = 1
!             do i = 1, ne
!                 v_mode(6*i+3) = L*i
!                 v_mode(6*i+6) = 1
!             end do
!     end select
!
! end subroutine mode_vector
!
! subroutine mode_vector_ver(v_mode, ndof, ne, L, mode, term)
!     ! Generate the beam mode vector based on the number of element
!     ! term = 0 : return position vector
!     ! term = 1 : return velocity vector
!     implicit none
!
!     integer, intent(in) :: ndof, ne, mode, term
!     double precision, dimension(ndof), intent(out) :: v_mode
!     double precision, intent(in) :: L
!     double precision :: L_beam
!     integer :: i
!
!     v_mode(:) = 0.0d0
!     L_beam = ne * L
!
!     do i = 1, ne
!         call mode_shape(v_mode(6*i+3), mode, L_beam, L*i, 0)
!     end do
!
!     select case (term)
!         case (0)
!             v_mode(6) = 1.0d0
!             do i = 1, ne
!                 v_mode(6*i+6) = 1.0d0
!             end do
!     end select
!
! end subroutine mode_vector_ver
!
! subroutine mode_shape(res, mode, L, x, der)
!     ! The mode shape function of cantilever beam (Euler beam theory)
!     ! der = 0 : mode shape function
!     ! der = 1 : derivative of mode shape function with respect to x
!
!     double precision, intent(out) :: res
!     integer, intent(in) :: mode, der
!     double precision, intent(in) :: L, x
!     double precision, dimension(9) :: kn_list
!     double precision :: kn
!
!     data kn_list / 1.875104068711961d0, 4.694091132974174d0, 7.854757438237613d0,&
!      10.995540734875467d0, 14.137168391046471d0, 17.278759532088237d0,&
!      20.420352251041251d0, 23.561944901806445d0, 26.703537555518299d0 /
!     kn = kn_list(mode) / L
!
!     select case (der)
!         case (0)
!             res = (DSIN(kn*L)+DSINH(kn*L))*(DCOS(kn*x)-DCOSH(kn*x))&
!              - (DCOS(kn*L)+DCOSH(kn*L))*(DSIN(kn*x)-DSINH(kn*x))
!         case (1)
!             res = -kn*(DSIN(kn*L)+DSINH(kn*L))*(DSIN(kn*x)+DSINH(kn*x))&
!              - kn*(DCOS(kn*L)+DCOSH(kn*L))*(DCOS(kn*x)-DCOSH(kn*x))
!     end select
!
! end subroutine mode_shape

subroutine gauss_points(n, x, w)
    ! The nodes and weights for Gauss quadrature integration in [-1, 1]
    ! x is nodes (in [-1, 1]) and w is weight to evaluate the integrals
    ! using Gauss quadrature of order n

    implicit none

    integer, intent(in) :: n
    double precision, dimension(n), intent(out) :: x, w

    select case (n)
        case (2)
            x = [-DSQRT(1.0d0/3.0d0), DSQRT(1.0d0/3.0d0)]
            w = [1.0d0, 1.0d0]
        case (3)
            x = [-DSQRT(3.0d0/5.0d0), 0.0d0, DSQRT(3.0d0/5.0d0)]
            w = [5.0d0/9.0d0, 8.0d0/9.0d0, 5.0d0/9.0d0]
        case (4)
            x = [-DSQRT(3.0d0/7.0d0+DSQRT(1.2d2)/3.5d1),&
                 -DSQRT(3.0d0/7.0d0-DSQRT(1.2d2)/3.5d1),&
                 DSQRT(3.0d0/7.0d0-DSQRT(1.2d2)/3.5d1),&
                 DSQRT(3.0d0/7.0d0+DSQRT(1.2d2)/3.5d1)]
            w = [1.0d0/2.0d0-5.0d0/(3.0d0*DSQRT(1.2d2)),&
                 1.0d0/2.0d0+5.0d0/(3.0d0*DSQRT(1.2d2)),&
                 1.0d0/2.0d0+5.0d0/(3.0d0*DSQRT(1.2d2)),&
                 1.0d0/2.0d0-5.0d0/(3.0d0*DSQRT(1.2d2))]
        case (5)
            x = [-(DSQRT(5.0d0+2.0d0*DSQRT(1.0d1/7.0d0)))/3.0d0,&
                 -(DSQRT(5.0d0-2.0d0*DSQRT(1.0d1/7.0d0)))/3.0d0,&
                 0.0d0, (DSQRT(5.0d0-2.0d0*DSQRT(1.0d1/7.0d0)))/3.0d0,&
                 (DSQRT(5.0d0+2.0d0*DSQRT(1.0d1/7.0d0)))/3.0d0]
            w = [(3.22d2-1.3d1*DSQRT(7.0d1))/9.0d2,&
                 (3.22d2+1.3d1*DSQRT(7.0d1))/9.0d2, 1.28d2/2.25d2,&
                 (3.22d2+1.3d1*DSQRT(7.0d1))/9.0d2,&
                 (3.22d2-1.3d1*DSQRT(7.0d1))/9.0d2]
    end select

end subroutine gauss_points

subroutine elastic_force(disp_cur, ndof, ne, E, A, Inertia, L, Qe, energy_e)
    ! Calculate the generalized elastic force for beam
    implicit none

    integer, intent(in) :: ndof, ne
    double precision, intent(in) :: E, A, Inertia, L
    double precision, dimension(ndof), intent(in) :: disp_cur
    double precision, dimension(ndof), intent(out) :: Qe
    double precision, intent(out) :: energy_e
    double precision, dimension(3) :: x3, w3   ! Gauss points for 3 node case
    double precision, dimension(5) :: x5, w5   ! Gauss points for 5 node case
    integer :: ie, istart, iend, k    ! Iterative integers
    double precision :: x                      ! Position
    double precision, dimension(12) :: Qn, Qf, en_int, ef_int  ! Element level force
    double precision :: energy_n, energy_f, energy_n_int, energy_f_int

    ! Calculate the nodes and weights for Gauss quadrature on the interval
    ! [-1, 1] using 5 nodes (for the axial deformation) and 3 nodes (for the
    ! flexural deformation)
    call gauss_points(3, x3, w3)
    call gauss_points(5, x5, w5)

    ! Initialize elastic force vector
    Qe(:) = 0.0d0

    ! Loop over all beam elements
    do ie = 1, ne
        ! Find range of nodal coordinates for the current element
        istart = 6 * ie - 5
        iend = 6 * ie + 6

        ! Calculate elastic force due to axial deformation
        Qn(:) = 0.0d0
        energy_n = 0.0d0
        do k = 1, 5
            x = (1.0d0 + x5(k)) * L / 2.0d0
            call axial_force_integrand(disp_cur(istart:iend), x, L, en_int, energy_n_int)
            Qn = Qn + w5(k) * en_int
            energy_n = energy_n + w5(k) * energy_n_int
        end do
        Qn = E * A * L / 2.0d0 * Qn
        energy_n = E * A * L / 4.0d0 * energy_n

        ! Calculate elastic force due to flexural deformation
        Qf(:) = 0.0d0
        energy_f = 0.0d0
        do k = 1, 3
            x = (1.0d0 + x3(k)) * L / 2.0d0
            call flexural_force_integrand(disp_cur(istart:iend), x, L, ef_int, energy_f_int)
            Qf = Qf + w3(k) * ef_int
            energy_f = energy_f + w3(k) * energy_f_int
        end do
        Qf = E * Inertia * L / 2.0d0 * Qf
        energy_f = E * Inertia * L / 4.0d0 * energy_f

        ! Add the two sources of elastic forces
        Qe(istart:iend) = Qe(istart:iend) + Qn + Qf
        energy_e = energy_n + energy_f
    end do

end subroutine elastic_force

subroutine axial_force_integrand(disp, x, L, en_int, energy_n_int)
    ! Calculate the kernel of integral of elastic force due to axial deformation

    implicit none

    double precision, dimension(12), intent(in) :: disp      ! Element displacement
    double precision, intent(in) :: x, L                     ! Position and length of element
    double precision, dimension(12), intent(out) :: en_int   ! Integrand
    double precision, intent(out) :: energy_n_int            ! Integrand for energy
    double precision, dimension(4) :: Sx                     ! Shape function
    double precision, dimension(3, 4) :: disp_rs             ! Reshape of displacement
    double precision, dimension(3) :: rx                     ! Position vector derivative
    double precision :: epsilon                              ! Axial strain
    double precision, dimension(12) :: epsilon_e             ! Axial strain derivative
    integer :: i, j
    double precision :: DDOT


    ! First derivative of shape function with respect to x
    call shape_fun(x, L, 1, Sx)

    ! First derivative of position vector r with respect to x
    disp_rs = RESHAPE(disp, (/ 3, 4 /))
    rx(:) = 0.0d0
    call DGEMV('n',3,4,1.0d0,disp_rs,3,Sx,1,1.0d0,rx,1)

    ! Axial strain of element, calculate from rx
    epsilon = (DDOT(3,rx,1,rx,1) - 1.0d0) / 2.0d0

    ! Partial derivative of axial strain to base vector e
    do i = 1, 4
        do j = 1, 3
            epsilon_e(3*(i-1)+j) = Sx(i) * rx(j)
        end do
    end do

    ! Integrand
    en_int = epsilon * epsilon_e
    energy_n_int = epsilon * epsilon

end subroutine axial_force_integrand

subroutine flexural_force_integrand(disp, x, L, ef_int, energy_f_int)
    ! Calculate the kernel of integral of elastic force due to flexural deformation

    implicit none

    double precision, dimension(12), intent(in) :: disp      ! Element displacement
    double precision, intent(in) :: x, L                     ! Position and length of element
    double precision, dimension(12), intent(out) :: ef_int   ! Integrand
    double precision, intent(out) :: energy_f_int            ! Integrand for energy
    double precision, dimension(4) :: Sx, Sxx                ! Shape function
    double precision, dimension(3, 4) :: disp_rs             ! Reshape of displacement
    double precision, dimension(3) :: rx, rxx, v             ! Position vector derivative
    double precision :: rx_mag, v_mag                        ! Magnitude of rx
    double precision, dimension(3,3) :: rx_tilde, rxx_tilde
    double precision, dimension(3, 12) :: Sx_mat, Sxx_mat
    integer :: i, j
    double precision :: DNRM2
    double precision, dimension(3,12) :: temp_1, temp_2
    double precision, dimension(12) :: temp_3, temp_4


    ! First and second derivative of shape function with respect to x
    call shape_fun(x, L, 1, Sx)
    call shape_fun(x, L, 2, Sxx)

    ! Shape function matrix
    Sx_mat(:,:) = 0.0d0
    Sxx_mat(:,:) = 0.0d0
    do i = 1, 4
        do j = 1, 3
            Sx_mat(j, 3*(i-1)+j) = Sx(i)
            Sxx_mat(j, 3*(i-1)+j) = Sxx(i)
        end do
    end do

    ! First and second derivative of position vector r with respect to x
    disp_rs = RESHAPE(disp, (/3, 4/))
    rx(:) = 0.0d0
    rxx(:) = 0.0d0
    call DGEMV('n',3,4,1.0d0,disp_rs,3,Sx,1,0.0d0,rx,1)
    call DGEMV('n',3,4,1.0d0,disp_rs,3,Sxx,1,0.0d0,rxx,1)
    rx_mag = DNRM2(3,rx,1)

    ! Calculate v = rx x rxx
    call cross(rx, rxx, v)
    v_mag = DNRM2(3,v,1)

    ! Evaluate the kernel of integral
    call tilde(rx, rx_tilde)
    call tilde(rxx, rxx_tilde)
    temp_1(:,:) = 0.0d0
    temp_2(:,:) = 0.0d0
    temp_3(:) = 0.0d0
    temp_4(:) = 0.0d0
    call DGEMM('n','n',3,12,3,1.0d0,rx_tilde,3,Sxx_mat,3,0.0d0,temp_1,3)
    call DGEMM('n','n',3,12,3,1.0d0,rxx_tilde,3,Sx_mat,3,0.0d0,temp_2,3)
    temp_1 = temp_1 - temp_2
    call DGEMV('t',3,12,1.0d0,temp_1,3,v,1,0.0d0,temp_3,1)
    call DGEMV('t',3,12,1.0d0,Sx_mat,3,rx,1,0.0d0,temp_4,1)
    ef_int = temp_3/rx_mag**6 - 3.0d0*v_mag**2/rx_mag**8*temp_4
    energy_f_int = v_mag**2 / rx_mag**6

end subroutine flexural_force_integrand

subroutine shape_fun(x, L, der, S)
    ! Generate the shape function matrix, and also the derivatives of shape
    ! function based on the third input of function 'der'
    ! Note 1 : input value of x should be between 0 and L
    ! Note 2 :
    !         In case der = 0, return S
    !         In case der = 1, return S_x
    !         In case der = 2, return S_xx

    implicit none

    double precision, intent(in) :: x, L
    integer, intent(in) :: der
    double precision, dimension(4), intent(out) :: S
    double precision :: xi

    xi = x / L

    select case (der)
        case (0)
            S(1) = 1.0d0 - 3.0d0*xi**2 + 2.0d0*xi**3
            S(2) = L * (xi - 2.0d0*xi**2 + xi**3)
            S(3) = 3.0d0*xi**2 - 2.0d0*xi**3
            S(4) = L * (-xi**2 + xi**3)
        case (1)
            S(1) = (6.0d0*xi**2 - 6.0d0*xi) / L
            S(2) = 1.0d0 - 4.0d0*xi + 3.0d0*xi**2
            S(3) = (-6.0d0*xi**2 + 6.0d0*xi) / L
            S(4) = -2.0d0*xi + 3.0d0*xi**2
        Case (2)
            S(1) = (1.2d1*xi - 6.0d0) / L**2
            S(2) = (-4.0d0 + 6.0*xi) / L
            S(3) = (6.0d0 - 1.2d1*xi) / L**2
            S(4) = (-2.0d0 + 6.0d0*xi) / L
    end select

end subroutine shape_fun

subroutine cross(v1, v2, res)
    ! Compute the cross product of two vectors

    implicit none

    double precision, dimension(3), intent(in) :: v1, v2
    double precision, dimension(3), intent(out) :: res

    res(1) = v1(2)*v2(3) - v1(3)*v2(2)
    res(2) = v1(3)*v2(1) - v1(1)*v2(3)
    res(3) = v1(1)*v2(2) - v1(2)*v2(1)

end subroutine cross

subroutine tilde(vec, vec_mat)
    ! Provide the matrix form of vector to help compute the cross product

    implicit none

    double precision, dimension(3), intent(in) :: vec
    double precision, dimension(3,3), intent(out) :: vec_mat

    vec_mat(1,1) = 0.0d0
    vec_mat(2,1) = vec(3)
    vec_mat(3,1) = -vec(2)
    vec_mat(1,2) = -vec(3)
    vec_mat(2,2) = 0.0d0
    vec_mat(3,2) = vec(1)
    vec_mat(1,3) = vec(2)
    vec_mat(2,3) = -vec(1)
    vec_mat(3,3) = 0.0d0

end subroutine tilde

! subroutine damping_force(ndof, M, vel_cur, vel_mode_hor, vel_mode_ver, Qd)
!     ! Assign damping force, which consider different damping ratio for
!     ! different frequency vibration mode
!     implicit none
!
!     integer, intent(in) :: ndof
!     double precision, dimension(ndof, ndof), intent(in) :: M
!     double precision, dimension(ndof), intent(out) :: Qd
!     double precision, dimension(ndof), intent(in) :: vel_cur
!     double precision, dimension(ndof) :: vel_temp_hor, vel_temp_ver
!     double precision, dimension(ndof) :: vel_temp_tot, vel_sum, vel_rem
!     integer, parameter :: nmode = 9
!     double precision, dimension(ndof, nmode), intent(in) :: vel_mode_hor, vel_mode_ver
!     double precision, parameter :: ratio = 1.0d1
!     double precision, dimension(nmode) :: alpha
!     ! double precision :: alpha_v
!     ! integer :: i
!     ! double precision :: DDOT
!
!     alpha(:) = ratio
!     alpha(1) = 5.0d4
!     ! alpha_v = 1.0d4
!     Qd(:) = 0.0d0
!     vel_sum(:) = 0.0d0
!
!     call DGEMV('n',ndof,ndof,ratio,M,ndof,vel_cur,1,1.0d0,Qd,1)
!
!     ! do i = 1, nmode
!     !     vel_temp_hor = DDOT(ndof, vel_cur, 1, vel_mode_hor(:,i), 1) * vel_mode_hor(:,i)
!     !     vel_temp_ver = DDOT(ndof, vel_cur, 1, vel_mode_ver(:,i), 1) * vel_mode_ver(:,i)
!     !     vel_temp_tot = vel_temp_hor + vel_temp_ver
!     !     vel_sum = vel_sum + vel_temp_tot
!     !     call DGEMV('n',ndof,ndof,alpha(i),M,ndof,vel_temp_tot,1,1.0d0,Qd,1)
!     !     ! call DGEMV('n',ndof,ndof,alpha_v,M,ndof,vel_temp_ver,1,1.0d0,Qd,1)
!     !     ! vel_sum = vel_sum + vel_temp_hor + vel_temp_ver
!     ! end do
!     ! vel_rem = vel_cur - vel_sum
!     ! call DGEMV('n',ndof,ndof,alpha(1),M,ndof,vel_rem,1,1.0d0,Qd,1)
!
! end subroutine damping_force

subroutine add_constraint_load(ndof, botCnstr, topCnstr, Q)
    ! Add constraint condition to load vector
    implicit none

    integer, intent(in) :: ndof, botCnstr, topCnstr
    double precision, dimension(ndof), intent(inout) :: Q

    if (botCnstr == 1) then
        Q(1:3) = 0.0d0
    else if (botCnstr == 2) then
        Q(1:6) = 0.0d0
    end if

    if (topCnstr == 1) then
        Q((ndof-5):(ndof-3)) = 0.0d0
    else if (topCnstr == 2) then
        Q((ndof-5):ndof) = 0.0d0
    end if

end subroutine add_constraint_load

! subroutine update_history(res_cur, res_prev, res_neg1, res_neg2, res_neg3, ndof, nbeam)
!     ! This subroutine update history for the purpose to time marching
!
!     implicit none
!
!     integer, intent(in) :: ndof, nbeam
!     double precision, dimension(ndof, nbeam), intent(in) :: res_cur
!     double precision, dimension(ndof, nbeam), intent(inout) :: res_prev, res_neg1, res_neg2, res_neg3
!
!     res_neg3 = res_neg2
!     res_neg2 = res_neg1
!     res_neg1 = res_prev
!     res_prev = res_cur
!
! end subroutine update_history

! subroutine point_load_profile(disp, params, tmaxdof, nbeam, f_out, ratio_out, nf_out, Qp)
!     ! This subroutine generate the point load profile based on current status of beams
!
!     use data_params
!     implicit none
!
!     integer, intent(in) :: tmaxdof, nbeam
!     ! integer, dimension(nbeam), intent(out) :: count
!     type(beam_param), dimension(nbeam), intent(in) :: params
!     double precision, dimension(tmaxdof, nbeam), intent(in) :: disp
!     double precision, dimension(tmaxdof, nbeam), intent(out) :: Qp
!     integer :: ie1, ie2, istart1, iend1, istart2, iend2
!     logical :: res
!     integer, parameter :: MAXFE = 250, MAXF = 1000
!     double precision, dimension(3, MAXFE) :: f
!     double precision, dimension(MAXFE) :: ratio1, ratio2
!     double precision, dimension(3, MAXF, nbeam) :: f_out
!     double precision, dimension(MAXF, nbeam) :: ratio_out
!     integer, dimension(nbeam) :: nf_out
!     integer :: i, nf, j, ii, jj
!     double precision :: pos_local
!     double precision, dimension(12) :: Qele
!     double precision, dimension(4) :: S
!     double precision, dimension(3, 12) :: S_mat
!
!
!     ! Initialize number of force
!     nf_out(:) = 0
!
!     ! Initialize force vector
!     Qp(:,:) = 0.0d0
!
!     ! Loop over all the beams for collision force
!     do i = 1, (nbeam-1)
!         ! Loop over all the elements from current beam and next beam
!         do ie1 = 1, params(i)%ne
!             ! Start and end number for current element in current beam
!             istart1 = 6 * (ie1 - 1) + 1
!             iend1 = 6 * (ie1 + 1)
!             do ie2 = 1, params(i+1)%ne
!                 istart2 = 6 * (ie2 - 1) + 1
!                 iend2 = 6 * (ie2 + 1)
!                 ! Use AABB to check if collision is possible between two beam elements
!                 call collision_tendency(disp(istart1:iend1,i), params(i)%x0,&
!                     params(i)%d, disp(istart2:iend2,i+1), params(i+1)%x0,&
!                     params(i+1)%d, res)
!                 if (res) then
!                     ! When collision tendency exists, determine if actual collision
!                     ! exists, and what will be impact force and position
!                     ! print *, ie1
!                     ! print *, ie2
!                     call collision_particle(disp(istart1:iend1,i), params(i),&
!                         disp(istart2:iend2,i+1), params(i+1), f, ratio1,&
!                         ratio2, nf)
!                     if (nf == 0) then
!                         CYCLE
!                     else
!                         f_out(:,(nf_out(i)+1):(nf_out(i)+nf),i) = f(:,1:nf)
!                         f_out(:,(nf_out(i+1)+1):(nf_out(i+1)+nf),i+1) = -f(:,1:nf)
!                         ratio_out((nf_out(i)+1):(nf_out(i)+nf),i) = (DBLE(ie1)&
!                             - 1.0d0 + ratio1(1:nf)) / DBLE(params(i)%ne)
!                         ratio_out((nf_out(i+1)+1):(nf_out(i+1)+nf),i+1) =&
!                             (DBLE(ie2) - 1.0d0 + ratio2(1:nf)) / DBLE(params(i+1)%ne)
!                         nf_out(i:i+1) = nf_out(i:i+1) + nf
!                     end if
!                 end if
!             end do
!         end do
!     end do
!
!     ! Loop over all the beams
!     do i = 1, nbeam
!         ! If there exist contact forces
!         if (nf_out(i) .ne. 0) then
!             ! Loop over forces
!             do j = 1, nf_out(i)
!                 ie1 = CEILING(ratio_out(j,i) * DBLE(params(i)%ne))
!                 pos_local = (ratio_out(j,i) * DBLE(params(i)%ne) - DBLE(ie1) + 1.0d0) * params(i)%L
!                 call shape_fun(pos_local, params(i)%L, 0, S)
!                 ! Shape function matrix
!                 S_mat(:,:) = 0.0d0
!                 do ii = 1, 4
!                     do jj = 1, 3
!                         S_mat(jj, 3*(ii-1)+jj) = S(ii)
!                     end do
!                 end do
!                 Qele(:) = 0.0d0
!                 call DGEMV('t',3,12,1.0d0,S_mat,3,f_out(:,j,i),1,1.0d0,Qele,1)
!                 ! Construct the global force vector
!                 Qp(6*ie1-5:6*ie1+6,i) = Qp(6*ie1-5:6*ie1+6,i) + Qele
!             end do
!         end if
!     end do
!
! end subroutine point_load_profile
!
! subroutine collision_tendency(disp1_in, x01, d1, disp2_in, x02, d2, res)
!     ! First layer of collision detection (AABB)
!
!     implicit none
!
!     double precision, dimension(12), intent(in) :: disp1_in, disp2_in
!     double precision, dimension(12) :: disp1, disp2
!     double precision, intent(in) :: x01, d1, x02, d2
!     double precision, dimension(4) :: b1, b2
!     logical, intent(out) :: res
!
!     disp1 = disp1_in
!     disp2 = disp2_in
!
!     ! Correct the value of position vector with initial position
!     if ((disp1(1) + x01) < (disp2(1) + x02)) then
!         disp1(1) = disp1(1) + x01 + d1/2.0d0
!         disp1(7) = disp1(7) + x01 + d1/2.0d0
!         disp2(1) = disp2(1) + x02 - d2/2.0d0
!         disp2(7) = disp2(7) + x02 - d2/2.0d0
!     else
!         disp1(1) = disp1(1) + x01 - d1/2.0d0
!         disp1(7) = disp1(7) + x01 - d1/2.0d0
!         disp2(1) = disp2(1) + x02 + d2/2.0d0
!         disp2(7) = disp2(7) + x02 + d2/2.0d0
!     end if
!
!     ! Get boundary of AABB box
!     call box_boundary(disp1, b1)
!     call box_boundary(disp2, b2)
!
!     if ((b1(2)>b2(1)) .and. (b1(1)<b2(2)) .and. (b1(4)>b2(3)) .and. (b1(3)<b2(4))) then
!         res = .true.
!     else
!         res = .false.
!     end if
!
! end subroutine collision_tendency
!
! subroutine box_boundary(disp, b)
!     ! Determine the boundary nodes of AABB box
!
!     implicit none
!
!     double precision, dimension(12), intent(in) :: disp
!     double precision, dimension(4), intent(out) :: b
!
!     b(1) = min(disp(1), disp(7))
!     b(2) = max(disp(1), disp(7))
!     b(3) = min(disp(3), disp(9))
!     b(4) = max(disp(3), disp(9))
!
! end subroutine box_boundary
!
! subroutine collision_particle(disp1_in, params1, disp2_in, params2, f, ratio1, ratio2, nf)
!     ! Determine the collision force between two beam elements by simulating beam
!     ! with particles
!
!     use data_params
!     implicit none
!
!     type(beam_param), intent(in) :: params1, params2
!     double precision, dimension(12), intent(in) :: disp1_in, disp2_in
!     double precision, dimension(12) :: disp1, disp2
!     integer, parameter :: MAXF = 250
!     integer, parameter :: MAXPAR = 1000
!     double precision, dimension(3, MAXF), intent(out) :: f
!     double precision, dimension(MAXF), intent(out) :: ratio1, ratio2
!     integer, intent(out) :: nf
!     double precision, parameter :: nu = 0.3d0
!     double precision :: E_star, R, K, dlim, dist, f_mag
!     integer :: n_par_1, n_par_2
!     double precision, dimension(MAXPAR) :: x_1
!     double precision, dimension(MAXPAR) :: x_2
!     integer :: i, j
!     double precision, dimension(3, MAXPAR) :: coord_1
!     double precision, dimension(3, MAXPAR) :: coord_2
!     double precision, dimension(3, 4) :: disp1_rs, disp2_rs
!     double precision, dimension(4) :: S
!     double precision, dimension(3) :: vec
!     double precision :: DNRM2
!
!     ! Correct the value of position vector with initial position
!     disp1 = disp1_in
!     disp2 = disp2_in
!     disp1(1) = disp1(1) + params1%x0
!     disp1(7) = disp1(7) + params1%x0
!     disp2(1) = disp2(1) + params2%x0
!     disp2(7) = disp2(7) + params2%x0
!
!     ! Separate beam element into particles
!     n_par_1 = CEILING(params1%L/params1%d)
!     n_par_2 = CEILING(params2%L/params2%d)
!     ! Position of particles in the beam element
!     x_1(n_par_1) = params1%L
!     x_2(n_par_2) = params2%L
!     do i = 1, (n_par_1-1)
!         x_1(i) = i * params1%d
!     end do
!     do i = 1, (n_par_2-1)
!         x_2(i) = i * params2%d
!     end do
!
!     ! Coordinate of particles
!     coord_1(:,:) = 0.0d0
!     coord_2(:,:) = 0.0d0
!     disp1_rs = RESHAPE(disp1, (/ 3, 4 /))
!     disp2_rs = RESHAPE(disp2, (/ 3, 4 /))
!     do i = 1, n_par_1
!         call shape_fun(x_1(i), params1%L, 0, S)
!         call DGEMV('n',3,4,1.0d0,disp1_rs,3,S,1,1.0d0,coord_1(:,i),1)
!     end do
!     do i = 1, n_par_2
!         call shape_fun(x_2(i), params2%L, 0, S)
!         call DGEMV('n',3,4,1.0d0,disp2_rs,3,S,1,1.0d0,coord_2(:,i),1)
!     end do
!
!     ! Contact stiffness of sphere
!     E_star = 1.0d-2 * (1.0d0/((1.0d0-nu**2)/params1%E + (1.0d0-nu**2)/params2%E))
!     R = 1.0d0 / (2.0d0/params1%d + 2.0d0/params2%d)
!     K = 4.0d0 / 3.0d0 * E_star * DSQRT(R)
!
!     ! Looping over particles, check collision
!     f(:,:) = 0.0d0
!     ratio1(:) = 0.0d0
!     ratio2(:) = 0.0d0
!     dlim = params1%d/2.0d0 + params2%d/2.0d0
!     nf = 0
!     do i = 1, n_par_1
!         do j = 1, n_par_2
!             dist = DNRM2(3,coord_1(:,i)-coord_2(:,j),1)
!             if (dist < dlim) then
!                 nf = nf + 1
!                 f_mag = K * (dlim - dist)**1.5
!                 vec = (coord_1(:,i) - coord_2(:,j)) / dist
!                 f(:,nf) = f_mag * vec
!                 ratio1(nf) = x_1(i)/params1%L
!                 ratio2(nf) = x_2(j)/params2%L
!                 ! print *, i
!                 ! print *, j
!                 ! print *, dist
!                 ! print *, f_mag
!                 ! print *, vec
!                 ! print *, DNRM2(3,vec,1)
!                 ! print *, nf
!                 ! print *, f(:,nf)
!             end if
!         end do
!     end do
!
! end subroutine collision_particle

subroutine point_load(disp, ndof, ne, E, dx0, L, thick, nbeam, Qp, Qimp, col_count)
    ! This subroutine generate the point load based on current status of beams
    implicit none

    integer, intent(in) :: ndof, ne, nbeam
    double precision, intent(in) :: E, dx0, L, thick
    integer, dimension(nbeam), intent(inout) :: col_count
    double precision, dimension(ndof, nbeam), intent(in) :: disp
    double precision, dimension(ndof, nbeam), intent(out) :: Qp
    double precision, dimension(2, nbeam), intent(out) :: Qimp
    integer :: i, ie1, ie2
    logical :: res
    double precision, dimension(2, ne) :: mu_cur, v1_cur, v2_cur
    double precision, dimension(2, ne) :: mu_next, v1_next, v2_next
    double precision, dimension(2, ne) :: v1_range_cur, v2_range_cur
    double precision, dimension(2, ne) :: v1_range_next, v2_range_next

    Qp(:,:) = 0.0d0
    Qimp(:,:) = 0.0d0

    ! Loop over all the beams for collision force
    do i = 1, (nbeam-1)
        ! Compute OBB information
        if (i == 1) then
            call OBB_info(disp(:,i), ndof, ne, dx0, L, 1, mu_cur, v1_cur,&
              v2_cur, v1_range_cur, v2_range_cur)
        else
            mu_cur = mu_next
            v1_cur = v1_next
            v2_cur = v2_next
            v1_range_cur = v1_range_next
            v2_range_cur = v2_range_next
        end if
        call OBB_info(disp(:,i+1), ndof, ne, dx0, L, i+1, mu_next, v1_next,&
          v2_next, v1_range_next, v2_range_next)
        ! Loop over all the elements from current beam and next beam
        do ie1 = 1, ne
            do ie2 = 1, ne
                ! Use OBB to check if collision iis possible
                call OBB_disjoint(mu_cur(:,ie1), v1_cur(:,ie1), v2_cur(:,ie1), v1_range_cur(:,ie1),&
                  v2_range_cur(:,ie1), mu_next(:,ie2), v1_next(:,ie2), v2_next(:,ie2),&
                  v1_range_next(:,ie2), v2_range_next(:,ie2), thick, res)
                if (.not. res) then
                    call collision_particle_level(disp(6*(ie1-1)+1:6*(ie1+1),i),&
                      disp(6*(ie2-1)+1:6*(ie2+1),i+1), Qp(6*ie1-5:6*ie1+6,i),&
                      Qp(6*ie2-5:6*ie2+6,i+1), Qimp(:,i), Qimp(:,i+1),&
                      col_count(i), col_count(i+1), E, thick, L, dx0, i)
                end if
            end do
        end do
    end do

end subroutine point_load

subroutine OBB_info(disp, ndof, ne, dx0, L, ibeam, mu, v1, v2, v1_range, v2_range)
    ! Compute OBB information
    implicit none

    integer, intent(in) :: ndof, ne, ibeam
    double precision, intent(in) :: dx0, L
    double precision, dimension(ndof), intent(in) :: disp
    double precision, dimension(2, ne), intent(out) :: mu, v1, v2
    double precision, dimension(2, ne), intent(out) :: v1_range, v2_range
    integer :: ie, icoord
    double precision, dimension(3, 4) :: disp_rs
    double precision, dimension(2, 5) :: coord, diff
    double precision, dimension(3) :: coord_cur
    double precision, dimension(5) :: cov
    double precision :: DDOT, prod
    double precision, dimension(4) :: S

    ! Loop over elements
    do ie = 1, ne
        ! Compute coordinates for selected nodes
        disp_rs = RESHAPE(disp(6*(ie-1)+1:6*(ie+1)), (/ 3, 4 /))
        disp_rs(1, 1) = disp_rs(1, 1) + DBLE(ibeam-1) * dx0
        disp_rs(1, 3) = disp_rs(1, 3) + DBLE(ibeam-1) * dx0
        coord(1, 1) = disp_rs(1, 1)
        coord(2, 1) = disp_rs(3, 1)
        coord(1, 5) = disp_rs(1, 3)
        coord(2, 5) = disp_rs(3, 3)
        do icoord = 1, 3
            coord_cur(:) = 0.0d0
            call shape_fun(0.25d0*DBLE(icoord)*L, L, 0, S)
            call DGEMV('n',3,4,1.0d0,disp_rs,3,S,1,1.0d0,coord_cur,1)
            coord(1, icoord+1) = coord_cur(1)
            coord(2, icoord+1) = coord_cur(3)
        end do
        ! Compute mean point and two axes for oriented box
        mu(:,ie) = SUM(coord, dim=2) / 5.0d0
        cov(:) = 0.0d0
        do icoord = 1, 5
            diff(:,icoord) = coord(:,icoord) - mu(:,ie)
            cov(1) = cov(1) + diff(1,icoord) * diff(1,icoord)
            cov(2) = cov(2) + diff(1,icoord) * diff(2,icoord)
            cov(3) = cov(3) + diff(2,icoord) * diff(2,icoord)
        end do
        cov = cov / 5.0d0
        cov(4) = ((cov(3)-cov(1)) + DSQRT((cov(3)-cov(1))**2 + 4.0d0*cov(2)**2)) / 2.0d0
        cov(5) = ((cov(3)-cov(1)) - DSQRT((cov(3)-cov(1))**2 + 4.0d0*cov(2)**2)) / 2.0d0
        if (cov(2) == 0.0d0) then
            v1(1,ie) = 1.0d0
            v1(2,ie) = 0.0d0
            v2(1,ie) = 0.0d0
            v2(2,ie) = 1.0d0
        else
            v1(1,ie) = cov(2) / DSQRT(cov(4)**2 + cov(2)**2)
            v1(2,ie) = cov(4) / DSQRT(cov(4)**2 + cov(2)**2)
            v2(1,ie) = cov(2) / DSQRT(cov(5)**2 + cov(2)**2)
            v2(2,ie) = cov(5) / DSQRT(cov(5)**2 + cov(2)**2)
        end if
        v1_range(:,ie) = 0.0d0
        v2_range(:,ie) = 0.0d0
        ! Compute the extreme range of box
        do icoord = 1, 5
            prod = DDOT(2,v1(:,ie),1,(coord(:,icoord)-mu(:,ie)),1)
            v1_range(1,ie) = MIN(v1_range(1,ie), prod)
            v1_range(2,ie) = MAX(v1_range(2,ie), prod)
            prod = DDOT(2,v2(:,ie),1,(coord(:,icoord)-mu(:,ie)),1)
            v2_range(1,ie) = MIN(v2_range(1,ie), prod)
            v2_range(2,ie) = MAX(v2_range(2,ie), prod)
        end do
    end do

end subroutine OBB_info

subroutine OBB_disjoint(mu1, v11, v21, v11_range, v21_range, mu2, v12, v22,&
    v12_range, v22_range, thick, res)

    implicit none

    double precision, dimension(2), intent(in) :: mu1, v11, v21, v11_range,&
        v21_range, mu2, v12, v22, v12_range, v22_range
    double precision, intent(in) :: thick
    logical, intent(out) ::  res
    double precision, dimension(2) :: v_mu

    res = .False.
    v_mu = mu2 - mu1
    call check_disjoint(v_mu, v11, v11_range, v12, v12_range, v22, v22_range, thick, res)
    if (res) then
      RETURN
    end if
    call check_disjoint(v_mu, v21, v21_range, v12, v12_range, v22, v22_range, thick, res)
    if (res) then
      RETURN
    end if
    call check_disjoint(v_mu, v12, v12_range, v11, v11_range, v21, v21_range, thick, res)
    if (res) then
      RETURN
    end if
    call check_disjoint(v_mu, v22, v22_range, v11, v11_range, v21, v21_range, thick, res)

end subroutine OBB_disjoint

subroutine check_disjoint(v_mu, v_base, range_base, v1, range1, v2, range2, thick, res)

    implicit none

    double precision, dimension(2), intent(in) :: v_mu, v_base, range_base,&
        v1, range1, v2, range2
    double precision, intent(in) :: thick
    logical, intent(out) :: res
    double precision :: d, r_base, r, DDOT, prod

    res = .False.
    d = DDOT(2,v_base,1,v_mu,1)
    r = 0.0d0
    if (d > 0) then
        r_base = range_base(2)
        prod = DDOT(2,v_base,1,v1,1)
        if (prod > 0) then
          r = r - range1(1) * prod
        else
          r = r - range1(2) * prod
        end if
        prod = DDOT(2,v_base,1,v2,1)
        if (prod > 0) then
          r = r - range2(1) * prod
        else
          r = r - range2(2) * prod
        end if
    else
        d = -d
        r_base = -range_base(1)
        prod = DDOT(2,v_base,1,v1,1)
        if (prod > 0) then
          r = r + range1(2) * prod
        else
          r = r + range1(1) * prod
        end if
        prod = DDOT(2,v_base,1,v2,1)
        if (prod > 0) then
          r = r + range2(2) * prod
        else
          r = r + range2(1) * prod
        end if
    end if

    if (d - r_base - r - thick > 0) then
        res = .True.
    end if

end subroutine check_disjoint

subroutine collision_particle_level(disp1_in, disp2_in, Qp1, Qp2, Qimp1, Qimp2,&
  count1, count2, E, thick, L, dx0, ibeam)
    ! Determine the collision force by simulating beam with particles
    implicit none

    double precision, intent(in) :: E, thick, L, dx0
    integer, intent(in) :: ibeam
    double precision, dimension(12), intent(in) :: disp1_in, disp2_in
    integer, intent(inout) :: count1, count2
    double precision, dimension(3, 4) :: disp1_rs, disp2_rs
    double precision, parameter :: nu = 0.3d0
    double precision :: E_star, R, K, dlim_L1, dlim_L2, dist, f_mag
    integer :: i1, j1, i2, j2, i3, j3, ii, jj
    integer, parameter :: ncenter = 5
    integer, parameter :: npar = 4
    double precision, dimension(3, ncenter) :: coord_L1_1, coord_L2_1
    double precision, dimension(3, ncenter) :: coord_L1_2, coord_L2_2
    double precision, dimension(3, npar) :: coord_L3_1, coord_L3_2
    double precision, dimension(4) :: S
    double precision, dimension(3) :: vec
    double precision :: DNRM2
    double precision, dimension(12), intent(inout) :: Qp1, Qp2
    double precision, dimension(2), intent(inout) :: Qimp1, Qimp2
    double precision :: pos
    double precision, dimension(3, 12) :: S_mat

    ! Correct the value of position vector with initial position
    disp1_rs = RESHAPE(disp1_in, (/ 3, 4 /))
    disp2_rs = RESHAPE(disp2_in, (/ 3, 4 /))
    disp1_rs(1,1) = disp1_rs(1,1) + DBLE(ibeam-1) * dx0
    disp1_rs(1,3) = disp1_rs(1,3) + DBLE(ibeam-1) * dx0
    disp2_rs(1,1) = disp2_rs(1,1) + DBLE(ibeam) * dx0
    disp2_rs(1,3) = disp2_rs(1,3) + DBLE(ibeam) * dx0

    ! Contact stiffness of sphere
    ! E_star = 1.0d-2 * (1.0d0/((1.0d0-nu**2)/E + (1.0d0-nu**2)/E))
    E_star = 1.0d-3 * E / 2.0d0 / (1.0d0-nu**2)
    ! R = 1.0d0 / (2.0d0/thick + 2.0d0/thick)
    R = thick / 4.0d0
    K = 4.0d0 / 3.0d0 * E_star * DSQRT(R)

    ! First level of overlap check
    coord_L1_1(:,:) = 0.0d0
    coord_L1_2(:,:) = 0.0d0
    do i1 = 1, ncenter
        call shape_fun((DBLE(i1)-0.5d0)*L/DBLE(ncenter), L, 0, S)
        call DGEMV('n',3,4,1.0d0,disp1_rs,3,S,1,1.0d0,coord_L1_1(:,i1),1)
        call shape_fun((DBLE(i1)-0.5d0)*L/DBLE(ncenter), L, 0, S)
        call DGEMV('n',3,4,1.0d0,disp2_rs,3,S,1,1.0d0,coord_L1_2(:,i1),1)
    end do
    ! print *, 'L1'
    ! print *, coord_L1_1
    ! print *, coord_L1_2
    dlim_L1 = L/DBLE(ncenter) + thick
    dlim_L2 = L/DBLE(ncenter*ncenter) + thick
    do i1 = 1, ncenter
      do j1 = 1, ncenter
        dist = DNRM2(3,coord_L1_1(:,i1)-coord_L1_2(:,j1),1)
          if (dist < dlim_L1) then
            ! Second level of overlap check
            coord_L2_1(:,:) = 0.0d0
            coord_L2_2(:,:) = 0.0d0
            do i2 = 1, ncenter
              call shape_fun(((DBLE(i2)-0.5d0)/DBLE(ncenter*ncenter)+DBLE(i1-1)/DBLE(ncenter))*L, L, 0, S)
              call DGEMV('n',3,4,1.0d0,disp1_rs,3,S,1,1.0d0,coord_L2_1(:,i2),1)
              call shape_fun(((DBLE(i2)-0.5d0)/DBLE(ncenter*ncenter)+DBLE(j1-1)/DBLE(ncenter))*L, L, 0, S)
              call DGEMV('n',3,4,1.0d0,disp2_rs,3,S,1,1.0d0,coord_L2_2(:,i2),1)
            end do
            do i2 = 1, ncenter
              do j2 = 1, ncenter
                dist = DNRM2(3,coord_L2_1(:,i2)-coord_L2_2(:,j2),1)
                if (dist < dlim_L2) then
                  ! Third level of overlap check (particle level)
                  coord_L3_1(:,:) = 0.0d0
                  coord_L3_2(:,:) = 0.0d0
                  do i3 = 1, npar
                    call shape_fun((DBLE(i3)-0.5d0)*thick+(DBLE(i2-1)/DBLE(ncenter*ncenter)+DBLE(i1-1)/DBLE(ncenter))*L, L, 0, S)
                    call DGEMV('n',3,4,1.0d0,disp1_rs,3,S,1,1.0d0,coord_L3_1(:,i3),1)
                    call shape_fun((DBLE(i3)-0.5d0)*thick+(DBLE(j2-1)/DBLE(ncenter*ncenter)+DBLE(j1-1)/DBLE(ncenter))*L, L, 0, S)
                    call DGEMV('n',3,4,1.0d0,disp2_rs,3,S,1,1.0d0,coord_L3_2(:,i3),1)
                  end do
                  ! print *, 'L3'
                  ! print *, i2, j2
                  ! print *, coord_L3_1
                  ! print *, coord_L3_2
                  do i3 = 1, npar
                    do j3 = 1, npar
                      dist = DNRM2(3,coord_L3_1(:,i3)-coord_L3_2(:,j3),1)
                      if (dist < thick) then
                        count1 = count1 + 1
                        count2 = count2 + 1
                        f_mag = K * (thick - dist)**1.5
                        vec = (coord_L3_1(:,i3) - coord_L3_2(:,j3)) / dist
                        pos = (DBLE(i3)-0.5d0)*thick+(DBLE(i2-1)/DBLE(ncenter*ncenter)+DBLE(i1-1)/DBLE(ncenter))*L
                        call shape_fun(pos, L, 0, S)
                        S_mat(:,:) = 0.0d0
                        do ii = 1, 4
                          do jj = 1, 3
                            S_mat(jj,3*(ii-1)+jj) = S(ii)
                          end do
                        end do
                        call DGEMV('t',3,12,f_mag,S_mat,3,vec,1,1.0d0,Qp1,1)
                        pos = (DBLE(j3)-0.5d0)*thick+(DBLE(j2-1)/DBLE(ncenter*ncenter)+DBLE(j1-1)/DBLE(ncenter))*L
                        call shape_fun(pos, L, 0, S)
                        S_mat(:,:) = 0.0d0
                        do ii = 1, 4
                          do jj = 1, 3
                            S_mat(jj,3*(ii-1)+jj) = S(ii)
                          end do
                        end do
                        call DGEMV('t',3,12,-f_mag,S_mat,3,vec,1,1.0d0,Qp2,1)
                        Qimp1(1) = Qimp1(1) + f_mag * vec(1)
                        Qimp1(2) = Qimp1(2) + f_mag * vec(3)
                        Qimp2(1) = Qimp2(1) - f_mag * vec(1)
                        Qimp2(2) = Qimp2(2) - f_mag * vec(3)
                      end if
                    end do
                  end do
                end if
              end do
            end do
          end if
        end do
      end do

end subroutine collision_particle_level

subroutine line_load_single(disp, vel, ndof, ne, L, width, dx0, ibeam, Ql, Qdrag,&
    Qam, t, Uv, Cv, kv, zv, tml)
    ! This subroutine generate the line load on the beam based on position
    ! and velocity of beam, also from fluid velocity, fluid density, drag
    ! coefficient
    implicit none

    integer, intent(in) :: ndof, ne, ibeam
    double precision, intent(in) :: L, width, dx0, Uv, Cv, kv, zv, tml
    double precision, dimension(ndof), intent(in) :: disp, vel
    double precision, dimension(ndof), intent(out) :: Ql
    double precision, dimension(12) :: Qele, Qele_temp
    double precision, dimension(2), intent(out) :: Qdrag, Qam
    double precision, dimension(2) :: Qdrag_ele, Qam_ele
    double precision, parameter :: rhof = 1000.0d0, Cm = 1.0d0, kw = 1.9904d0, h = 0.30d0
    double precision :: Cd
    integer, parameter :: np = 5
    integer :: i, ie, istart, iend, ii, jj
    double precision, dimension(4, np) :: S, Sx
    double precision, dimension(3, 12, np) :: S_mat, Sx_mat
    double precision :: x
    double precision, dimension(3) :: r_c, v_c, rx_c, v_f, v_rel, fl_c, v_reln, a_fn
    double precision, dimension(3) :: f_am, f_tot, a_f
    double precision, parameter :: pi = 4.0d0 * DATAN(1.0d0)
    double precision, dimension(5) :: x5, w5
    double precision :: DNRM2, DDOT
    double precision, parameter :: v_vortex = 0.5d0
    double precision, intent(in) :: t

    ! omega = 2.0d0 * pi / Tw
    ! Uw = 9.8d0 * kw * aw / omega
    ! KC = Uw * Tw / width
    ! Cd = MAX(10.0d0 * KC**(-1.0d0/3.0d0), 1.95d0)
    Cd = 3.6d0
    ! The nodes and weights for Gauss quadrature on the interval [-1, 1]
    call gauss_points(np, x5, w5)

    ! Compute shape function and derivative of shape function
    S_mat(:,:,:) = 0.0d0
    Sx_mat(:,:,:) = 0.0d0
    do i = 1, np
        x = (1.0d0 + x5(i)) / 2.0d0
        call shape_fun(x * L, L, 0, S(:,i))
        call shape_fun(x * L, L, 1, Sx(:,i))
        do ii = 1, 4
            do jj = 1, 3
                S_mat(jj, 3*(ii-1)+jj, i) = S(ii,i)
                Sx_mat(jj, 3*(ii-1)+jj, i) = Sx(ii,i)
            end do
        end do
    end do

    ! Initialize external force vector
    Ql(:) = 0.0d0
    Qdrag(:) = 0.0d0
    Qam(:) = 0.0d0

    ! Loop over all beam elements
    do ie = 1, ne
        istart = 6 * (ie - 1) + 1
        iend = 6 * (ie + 1)
        ! Line force of element
        Qele(:) = 0.0d0
        Qdrag_ele(:) = 0.0d0
        Qam_ele(:) = 0.0d0
        ! Loop over all element points
        do i = 1, np
            ! Position vector and derivative
            r_c(:) = 0.0d0
            call DGEMV('N',3,12,1.0d0,S_mat(:,:,i),3,disp(istart:iend),&
                1,0.0d0,r_c,1)
            r_c(1) = r_c(1) + DBLE(ibeam-1) * dx0
            ! Velocity vector
            v_c(:) = 0.0d0
            rx_c(:) = 0.0d0
            call DGEMV('N',3,12,1.0d0,S_mat(:,:,i),3,vel(istart:iend),&
                1,0.0d0,v_c,1)
            call DGEMV('N',3,12,1.0d0,Sx_mat(:,:,i),3,disp(istart:iend),&
                1,0.0d0,rx_c,1)
            ! Fluid velocity
            v_f = 0.0d0
            a_f = 0.0d0
            if (r_c(3) > (zv-0.5d0*tml) .and. r_c(1) - Cv * t <= 0.0d0) then
              v_f(1) = Uv * DCOS(kv*pi*(r_c(1)-Cv*t)) * DSIN(kv*pi*(r_c(3)-zv))
              v_f(3) = -Uv * DSIN(kv*pi*(r_c(1)-Cv*t)) * DCOS(kv*pi*(r_c(3)-zv))
              a_f(1) = Cv*Uv*DSIN(kv*pi*(r_c(1)-Cv*t)) * DSIN(kv*pi*(r_c(3)-zv))
              a_f(3) = Cv*Uv*DCOS(kv*pi*(r_c(1)-Cv*t)) * DCOS(kv*pi*(r_c(3)-zv))
            end if
            if (r_c(1) - Cv * t > 0.0d0) then
              v_f(1) = 0.05d0
            end if

            ! Relative velocity
            v_rel = v_f - v_c
            v_reln = v_rel - DDOT(3,v_rel,1,rx_c,1) / DDOT(3,rx_c,1,rx_c,1) * rx_c
            fl_c = 0.5d0 * rhof * Cd * width * DNRM2(3,v_reln,1) * v_reln
            ! Added mass
            a_fn = a_f - DDOT(3,a_f,1,rx_c,1) / DDOT(3,rx_c,1,rx_c,1) * rx_c
            f_am = rhof * Cm * pi * width**2 / 4.0d0 * a_fn
            f_tot = fl_c + f_am
            if (t <= 1.0d0) then
              f_tot = t * f_tot
            end if
            ! Compute line force of element using Gauss integration
            Qele_temp(:) = 0.0d0
            call DGEMV('T',3,12,1.0d0,S_mat(:,:,i),3,f_tot,1,0.0d0,Qele_temp,1)
            Qele = Qele + w5(i) * Qele_temp
            Qdrag_ele(1) = Qdrag_ele(1) + w5(i) * fl_c(1)
            Qdrag_ele(2) = Qdrag_ele(2) + w5(i) * fl_c(3)
            Qam_ele(1) = Qam_ele(1) + w5(i) * f_am(1)
            Qam_ele(2) = Qam_ele(2) + w5(i) * f_am(3)
        end do
        ! Rescale line force
        Qele = L / 2.0d0 * Qele
        Qdrag_ele = L / 2.0d0 * Qdrag_ele
        Qam_ele = L / 2.0d0 * Qam_ele
        ! Construct the global force vector
        Ql(istart:iend) = Ql(istart:iend) + Qele
        Qdrag = Qdrag + Qdrag_ele
        Qam = Qam + Qam_ele
    end do

end subroutine line_load_single


subroutine runge_kutta(h, disp_prev, vel_prev, acc_prev, ndof, ne, botCnstr, topCnstr, E,&
    A, Inertia, L, M_L, M_U, Qg, Qp, disp, vel,&
    acc, Q, Qdrag, Qam, energy_e, pivot, width, dx0, ibeam, t, Uv, Cv, kv, zv, tml)
    ! Compute new system status with Runge Kutta method
    ! Using four stage, fourth order accurate method
    implicit none

    integer, intent(in) :: ndof, ne, botCnstr, topCnstr, ibeam
    double precision, intent(in) :: E, A, Inertia, L, h, width, dx0, t, Uv, Cv, kv, zv, tml
    double precision, dimension(ndof), intent(in) :: disp_prev, vel_prev, acc_prev, Qg, Qp
    double precision, dimension(ndof, ndof), intent(in) :: M_L, M_U
    double precision, dimension(ndof), intent(out) :: disp, vel, acc, Q
    double precision, dimension(2), intent(out) :: Qdrag, Qam
    double precision, intent(out) :: energy_e
    double precision, dimension(ndof) :: disp_xi_1, vel_xi_1, acc_xi_1, Q_xi_1
    double precision, dimension(ndof) :: disp_xi_2, vel_xi_2, acc_xi_2, Q_xi_2
    double precision, dimension(ndof) :: disp_xi_3, vel_xi_3, acc_xi_3, Q_xi_3
    double precision, dimension(ndof) :: disp_xi_4, vel_xi_4
    double precision, dimension(ndof) :: Qe_temp, Ql_temp
    double precision :: energy_e_temp
    integer, dimension(ndof), intent(in) :: pivot

    ! Compute the status at 'xi_1' step
    disp_xi_1 = h * vel_prev
    vel_xi_1 = h * acc_prev
    call elastic_force(disp_prev+0.5d0*disp_xi_1, ndof, ne, E, A, Inertia, L, Qe_temp, energy_e_temp)
    call line_load_single(disp_prev+0.5d0*disp_xi_1, vel_prev+0.5d0*vel_xi_1, ndof, ne, L, width,&
     dx0, ibeam, Ql_temp, Qdrag, Qam, t+0.5d0*h, Uv, Cv, kv, zv, tml)
    Q_xi_1 = Qg + Qp + Ql_temp - Qe_temp
    call add_constraint_load(ndof, botCnstr, topCnstr, Q_xi_1)
    call linear_system(M_L, M_U, Q_xi_1, acc_xi_1, pivot, ndof)

    ! Compute the status at 'xi_2' step
    disp_xi_2 = h * (vel_prev + 0.5d0 * vel_xi_1)
    vel_xi_2 = h * acc_xi_1
    call elastic_force(disp_prev+0.5d0*disp_xi_2, ndof, ne, E, A, Inertia, L, Qe_temp, energy_e_temp)
    call line_load_single(disp_prev+0.5d0*disp_xi_2, vel_prev+0.5d0*vel_xi_2, ndof, ne, L, width,&
     dx0, ibeam, Ql_temp, Qdrag, Qam, t+0.5d0*h, Uv, Cv, kv, zv, tml)
    Q_xi_2 = Qg + Qp + Ql_temp - Qe_temp
    call add_constraint_load(ndof, botCnstr, topCnstr, Q_xi_2)
    call linear_system(M_L, M_U, Q_xi_2, acc_xi_2, pivot, ndof)

    ! Compute the status at 'xi_3' step
    disp_xi_3 = h * (vel_prev + 0.5d0 * vel_xi_2)
    vel_xi_3 = h * acc_xi_2
    call elastic_force(disp_prev+disp_xi_3, ndof, ne, E, A, Inertia, L, Qe_temp, energy_e_temp)
    call line_load_single(disp_prev+disp_xi_3, vel_prev+vel_xi_3, ndof, ne, L, width,&
     dx0, ibeam, Ql_temp, Qdrag, Qam, t+h, Uv, Cv, kv, zv, tml)
    Q_xi_3 = Qg + Qp + Ql_temp - Qe_temp
    call add_constraint_load(ndof, botCnstr, topCnstr, Q_xi_3)
    call linear_system(M_L, M_U, Q_xi_3, acc_xi_3, pivot, ndof)

    ! Compute the status at 'xi_4' step
    disp_xi_4 = h * (vel_prev + vel_xi_3)
    vel_xi_4 = h * acc_xi_3

    ! Compute the status at next time step
    disp = disp_prev + 1.0d0/6.0d0*(disp_xi_1+2.0d0*disp_xi_2+2.0d0*disp_xi_3+disp_xi_4)
    vel = vel_prev + 1.0d0/6.0d0*(vel_xi_1+2.0d0*vel_xi_2+2.0d0*vel_xi_3+vel_xi_4)
    call elastic_force(disp, ndof, ne, E, A, Inertia, L, Qe_temp, energy_e)
    call line_load_single(disp, vel, ndof, ne, L, width, dx0, ibeam, Ql_temp,&
     Qdrag, Qam, t+h, Uv, Cv, kv, zv, tml)
    Q = Qg + Qp + Ql_temp - Qe_temp
    call add_constraint_load(ndof, botCnstr, topCnstr, Q)
    call linear_system(M_L, M_U, Q, acc, pivot, ndof)

end subroutine runge_kutta

! subroutine ABAM4(h, disp_prev, vel_neg3, vel_neg2, vel_neg1, vel_prev, Q_neg3,&
!     Q_neg2, Q_neg1, Q_prev, ndof, ne, botCnstr, topCnstr, E, A, Inertia, L, M, M_L,&
!     M_U, Qg, Qp, Ql, vel_mode_hor, vel_mode_ver, disp, vel, acc, Q, energy_e, pivot)
!     ! Compute new system status with 4th order Adams-Bashforth Adams-Moulton scheme
!     implicit none
!
!     integer, intent(in) :: ndof, ne, botCnstr, topCnstr
!     double precision, intent(in) :: E, A, Inertia, L
!     double precision, dimension(ndof), intent(in) :: disp_prev
!     double precision, dimension(ndof), intent(in) :: vel_neg3, vel_neg2, vel_neg1, vel_prev
!     double precision, dimension(ndof), intent(in) :: Q_neg3, Q_neg2, Q_neg1, Q_prev
!     double precision, dimension(ndof), intent(in) :: Qg, Qp, Ql
!     double precision, intent(in) :: h
!     double precision, dimension(ndof, ndof), intent(in) :: M, M_L, M_U
!     double precision, dimension(ndof, 9), intent(in) :: vel_mode_hor, vel_mode_ver
!     double precision, dimension(ndof), intent(out) :: disp, vel, acc, Q
!     double precision, dimension(ndof) :: disp_next_bar, vel_next_bar, Q_next_bar, Qtemp, Qe, Qd, acc_temp
!     integer, dimension(ndof), intent(in) :: pivot
!     double precision, intent(out) :: energy_e
!
!     disp_next_bar = disp_prev + h/2.4d1 * (5.5d1*vel_prev - 5.9d1*vel_neg1 + 3.7d1*vel_neg2 - 9.0d0*vel_neg3)
!     Qtemp = (5.5d1*Q_prev - 5.9d1*Q_neg1 + 3.7d1*Q_neg2 - 9.0d0*Q_neg3)/2.4d1
!     call linear_system(M_L, M_U, Qtemp, acc_temp, pivot, ndof)
!     ! acc_temp = Qtemp
!     vel_next_bar = vel_prev + h * acc_temp
!     call elastic_force(disp_next_bar, ndof, ne, E, A, Inertia, L, Qe, energy_e)
!     ! call damping_force(ndof, M, vel_next_bar, vel_mode_hor, vel_mode_ver, Qd)
!     Qd = 0.0d0
!     Q_next_bar = Qg + Qp + Ql - Qe - Qd
!     call add_constraint_load(ndof, botCnstr, topCnstr, Q_next_bar)
!
!     disp = disp_prev + h/2.4d1 * (9.0d0*vel_next_bar + 1.9d1*vel_prev - 5.0d0*vel_neg1 + vel_neg2)
!     Qtemp = (9.0d0*Q_next_bar + 1.9d1*Q_prev - 5.0d0*Q_neg1 + Q_neg2) / 2.4d1
!     call linear_system(M_L, M_U, Qtemp, acc_temp, pivot, ndof)
!     ! acc_temp = Qtemp
!     vel = vel_prev + h * acc_temp
!     call elastic_force(disp, ndof, ne, E, A, Inertia, L, Qe, energy_e)
!     ! call damping_force(ndof, M, vel, vel_mode_hor, vel_mode_ver, Qd)
!     Qd = 0.0d0
!     Q = Qg + Qp + Ql - Qe - Qd
!     call add_constraint_load(ndof, botCnstr, topCnstr, Q)
!     call linear_system(M_L, M_U, Q, acc, pivot, ndof)
!
! end subroutine ABAM4

subroutine write_output(beam_folder, OUTUNIT, t, disp, vel, acc, col_count, ene_e,&
   ene_t, qdrag, qam, qimp, n)
    ! This subroutine write data to output files

    implicit none

    character(len=8), intent(in) :: beam_folder
    integer, intent(in) :: OUTUNIT
    double precision, intent(in) :: t, ene_e, ene_t
    ! double precision, intent(in) :: t
    integer, intent(in) :: col_count
    double precision, dimension(2), intent(in) :: qdrag, qam, qimp
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: disp, vel, acc
    integer :: idof, ios

    ! Open file to save displacement information
    open (unit = OUTUNIT, file = 'result/'//beam_folder//'/position_result.txt',&
        iostat = ios, status = 'old', action = 'write', position = 'append')
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') t,','
    do idof = 1, n-1
        write (OUTUNIT, '(E20.12E2, A)', advance = 'no') disp(idof),','
    end do
    write (OUTUNIT, '(E20.12E2)') disp(n)
    close (OUTUNIT, status = 'keep')

    ! Open file to save velocity information
    open (unit = OUTUNIT, file = 'result/'//beam_folder//'/velocity_result.txt',&
        iostat = ios, status = 'old', action = 'write', position = 'append')
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') t,','
    do idof = 1, n-1
        write (OUTUNIT, '(E20.12E2, A)', advance = 'no') vel(idof),','
    end do
    write (OUTUNIT, '(E20.12E2)') vel(n)
    close (OUTUNIT, status = 'keep')

    ! Open file to save acceleration information
    open (unit = OUTUNIT, file = 'result/'//beam_folder//'/acceleration_result.txt',&
        iostat = ios, status = 'old', action = 'write', position = 'append')
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') t,','
    do idof = 1, n-1
        write (OUTUNIT, '(E20.12E2, A)', advance = 'no') acc(idof),','
    end do
    write (OUTUNIT, '(E20.12E2)') acc(n)
    close (OUTUNIT, status = 'keep')

    ! Open file to save collision count, elastic and kinetic energy
    open (unit = OUTUNIT, file = 'result/'//beam_folder//'/collision_and_energy.txt',&
        iostat = ios, status = 'old', action = 'write', position = 'append')
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') t,','
    write (OUTUNIT, '(I20, A)', advance = 'no') col_count,','
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') ene_e,','
    write (OUTUNIT, '(E20.12E2)') ene_t
    close (OUTUNIT, status = 'keep')

    ! Open file to save drag force
    open (unit = OUTUNIT, file = 'result/'//beam_folder//'/drag_force.txt',&
        iostat = ios, status = 'old', action = 'write', position = 'append')
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') t,','
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') qdrag(1),','
    write (OUTUNIT, '(E20.12E2)') qdrag(2)
    close (OUTUNIT, status = 'keep')

    ! Open file to save added mass force
    open (unit = OUTUNIT, file = 'result/'//beam_folder//'/added_mass_force.txt',&
        iostat = ios, status = 'old', action = 'write', position = 'append')
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') t,','
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') qam(1),','
    write (OUTUNIT, '(E20.12E2)') qam(2)
    close (OUTUNIT, status = 'keep')

    ! Open file to save added mass force
    open (unit = OUTUNIT, file = 'result/'//beam_folder//'/impact_force.txt',&
        iostat = ios, status = 'old', action = 'write', position = 'append')
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') t,','
    write (OUTUNIT, '(E20.12E2, A)', advance = 'no') qimp(1),','
    write (OUTUNIT, '(E20.12E2)') qimp(2)
    close (OUTUNIT, status = 'keep')

end subroutine write_output

subroutine linear_system(L, U, b, X, piv, n)
    ! Solving linear system
    implicit none

    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: piv
    double precision, dimension(n, n), intent(in) :: L, U
    double precision, dimension(n), intent(in) :: b
    double precision, dimension(n), intent(out) :: X
    double precision, dimension(n) :: btr, Y
    double precision :: temp
    integer :: i

    ! Transfer the vector based on pivot array
    btr = b
    do i = 1, n
        temp = btr(i)
        btr(i) = btr(piv(i))
        btr(piv(i)) = temp
    end do

    ! Compute array of Y
    Y = btr
    call DTRSV('L','N','U',n,L,n,Y,1)

    ! Compute array of X
    X = Y
    call DTRSV('U','N','N',n,U,n,X,1)

end subroutine linear_system
