module md_pro 
use machina_basic, only: i4, f8
use basis 
use read_input
use cint
use scf
use gradient
use output
implicit none 

real(kind=f8), parameter :: pi = 4.0_f8*atan(1.0_f8)
real(kind=f8), parameter :: kB = 1.38064852D-23
real(kind=f8), parameter :: J2au = 2.2937122783963248D+17
real(kind=f8), parameter :: fs2au = 4.1341373335182112D+01
real(kind=f8), parameter :: dalton2au = 1822.888486209
real(kind=f8) :: scf_start_time
real(kind=f8) :: scf_end_time
real(kind=f8) :: grad_start_time
real(kind=f8) :: grad_end_time
real(kind=f8) :: scf_time 
real(kind=f8) :: grad_time
type(cint_info) :: cint_in_md
type(atom_basis_info),dimension(:),allocatable :: bas_in_md
type(scf_info) :: scf_in_md
type(mole_input_info) :: mole_in_md
type(grad_info) :: grad_in_md

type,public :: md 
    integer :: nstep
    character(len=200) :: job
    real(kind=f8) :: dt
    real(kind=f8) :: init_temp
    real(kind=f8) :: bath_temp
    real(kind=f8) :: con_time
    real(kind=f8) :: temp_cal
    real(kind=f8) :: kin
    real(kind=f8) :: pot
    real(kind=f8) :: zeta_old
    real(kind=f8) :: zeta_new
    real(kind=f8),allocatable,dimension(:,:) :: pos_old
    real(kind=f8),allocatable,dimension(:,:) :: pos_new
    real(kind=f8),allocatable,dimension(:,:) :: vel_old
    real(kind=f8),allocatable,dimension(:,:) :: vel_new
    real(kind=f8),allocatable,dimension(:,:) :: vel_half
    real(kind=f8),allocatable,dimension(:,:) :: acc_old
    real(kind=f8),allocatable,dimension(:,:) :: acc_new
    real(kind=f8),allocatable,dimension(:,:) :: bond_length
    real(kind=f8),allocatable,dimension(:,:,:) :: bond_angle

    contains
        procedure,public :: set_md_info => set_md_info_sub
        procedure,public :: init_md => init_md_sub
        procedure,public :: random_guess => random_guess_sub
        procedure,public :: md_int => md_int_sub
        procedure,public :: cal_kin_pot => cal_kin_pot_sub
        procedure,public :: cal_temp => cal_temp_sub
        procedure,public :: berendsen => berendsen_sub
        procedure,public :: xyz_out => xyz_out_sub
        procedure,public :: E_out => E_out_sub
        procedure,public :: calc_angle => calc_angle_sub
        procedure,public :: pop_out => pop_out_sub
end type 

contains

subroutine set_md_info_sub(this,mole_in_md)
implicit none
    class(md) :: this
    type(mole_input_info),intent(in) :: mole_in_md
    integer :: status

    this%nstep = mole_in_md%nstep
    this%dt = mole_in_md%dt*fs2au
    this%init_temp = mole_in_md%init_temp
    this%bath_temp = mole_in_md%bath_temp
    this%con_time = mole_in_md%con_time*fs2au
    this%job = mole_in_md%md_job
    allocate(this%pos_old(mole_in_md%natom,3),stat=status)
    allocate(this%pos_new(mole_in_md%natom,3),stat=status)
    allocate(this%vel_old(mole_in_md%natom,3),stat=status)
    allocate(this%vel_new(mole_in_md%natom,3),stat=status)
    allocate(this%vel_half(mole_in_md%natom,3),stat=status)
    allocate(this%acc_old(mole_in_md%natom,3),stat=status)
    allocate(this%acc_new(mole_in_md%natom,3),stat=status)
    allocate(this%bond_length(mole_in_md%natom,mole_in_md%natom))
    allocate(this%bond_angle(mole_in_md%natom,mole_in_md%natom,mole_in_md%natom))

end subroutine

subroutine random_guess_sub(this,random)
implicit none
    class(md) :: this
    real(kind=f8),intent(inout) :: random
    real(kind=f8) :: z1,z2

    call random_seed()
    call random_number(z1)
    call random_seed()
    call random_number(z2)
    random = sqrt(-2.0_f8*log(z1))*sin(2.0_f8*pi*z2)
end subroutine

subroutine init_md_sub(this,mole_in_md,cint_in_md)
implicit none
    class(md) :: this
    type(mole_input_info),intent(in) :: mole_in_md
    type(cint_info),intent(in) :: cint_in_md
    real(kind=f8),dimension(3) :: avg_vel
    real(kind=f8),allocatable,dimension(:,:) :: vel_random
    real(kind=f8) :: vtemp
    real(kind=f8) :: factor
    integer :: i,status

    allocate(vel_random(mole_in_md%natom,3),stat=status)
    avg_vel = 0.0_f8
    do i = 1,mole_in_md%natom    
        call this%random_guess(vel_random(i,1))
        call this%random_guess(vel_random(i,2))
        call this%random_guess(vel_random(i,3))
        avg_vel(:) = avg_vel(:) + vel_random(i,:)
    end do 
    avg_vel(:) = avg_vel(:)/mole_in_md%natom
    !subracting the average velocity so the center of mass does not move
    do i = 1, mole_in_md%natom
        vel_random(i,1) = vel_random(i,1) - avg_vel(1)
        vel_random(i,2) = vel_random(i,2) - avg_vel(2)
        vel_random(i,3) = vel_random(i,3) - avg_vel(3)
    end do  
    do i = 1, mole_in_md%natom
        this%vel_old(i,:) = sqrt((this%init_temp*kB*J2au)/(cint_in_md%mass(i)*dalton2au))*vel_random(i,:)
        this%pos_old(i,1) = mole_in_md%x_list(i)/0.529177249_f8
        this%pos_old(i,2) = mole_in_md%y_list(i)/0.529177249_f8
        this%pos_old(i,3) = mole_in_md%z_list(i)/0.529177249_f8
    end do 
    deallocate(vel_random)
end subroutine

subroutine cal_kin_pot_sub(this,cint_in_md,scf_in_md,vel)
implicit none
    class(md) :: this
    type(cint_info),intent(in) :: cint_in_md
    type(scf_info),intent(in) :: scf_in_md
    real(kind=f8),allocatable,dimension(:,:),intent(in) :: vel
    real(kind=f8) :: vtemp
    integer :: i 

    vtemp = 0.0_f8
    do i=1,cint_in_md%natom 
        vtemp = vtemp + dalton2au*cint_in_md%mass(i)*((vel(i,1))**2.0 &
                                                     +(vel(i,2))**2.0 &
                                                     +(vel(i,3))**2.0)
    end do 
    this%kin = 0.5_f8 * vtemp 
    this%pot = scf_in_md%Ent - this%kin 
end subroutine

subroutine cal_temp_sub(this,cint_in_md)
implicit none
    class(md) :: this
    type(cint_info),intent(in) :: cint_in_md

    this%temp_cal = (2.0*this%kin)/(3.0*kB*J2au*cint_in_md%natom)
end subroutine

subroutine berendsen_sub(this,cint_in_md)
implicit none
    class(md) :: this
    type(cint_info),intent(in) :: cint_in_md
    real(kind=f8) :: f

    call this%cal_temp(cint_in_md)
    f = sqrt(1.0_f8+(this%dt*(this%bath_temp/this%temp_cal-1.0_f8)/this%con_time))
    if (f > 1.1_f8) then
        f = 1.1_f8
    elseif (f < 0.9_f8) then 
        f = 0.9_f8
    end if
    this%vel_new(:,:) = f * this%vel_new(:,:)
end subroutine

subroutine md_int_sub(this,mole_in_md,cint_in_md,bas_in_md,scf_in_md,grad_in_md)
implicit none
    class(md) :: this
    type(mole_input_info),intent(inout) :: mole_in_md
    type(cint_info),intent(inout) :: cint_in_md
    type(atom_basis_info),dimension(:),allocatable,intent(in) :: bas_in_md
    type(scf_info),intent(inout) :: scf_in_md
    type(grad_info),intent(in) :: grad_in_md
    real(kind=f8),dimension(3) :: vel_cent
    real(kind=f8),dimension(3) :: mv
    real(kind=f8) :: mtot
    real(kind=f8) :: Etemp1
    real(kind=f8) :: Etemp2
    integer :: i,n

    write(*,"(/,a,i5,1x,a,/)") '===== ',0,'AIMD-step  ====='
    ! set basis and cint infomation
    call cint_in_md%set_mole_bas_info(bas_in_md,mole_in_md)
    call cint_in_md%set_cint_info(mole_in_md,1)
    ! init the md
    call this%set_md_info(mole_in_md)
    write(*,"(1x,a)") 'Initialize MD...'
    call this%init_md(mole_in_md,cint_in_md)
    ! start SCF
    write(*,"(1x,a)") 'Start the ab initio SCF...'
    call cpu_time(scf_start_time)
    if (trim(mole_in_md%QC_method)=='hf') then
        call cint_in_md%cal_cint()
        call scf_in_md%scf_pro(cint_in_md,mole_in_md)
        call energy_out(scf_in_md%Ent,scf_in_md%iter_cong)
        Etemp1 = scf_in_md%Ent 
    else if (trim(mole_in_md%QC_method)=='mp2(fc)' .or. trim(mole_in_md%QC_method)=='mp2(full)') then 
        call cint_in_md%cal_cint()
        call scf_in_md%scf_pro(cint_in_md,mole_in_md)
        call scf_in_md%init_mp2(cint_in_md)
        call scf_in_md%cal_mp2(cint_in_md)
        call out_mp2(scf_in_md%Emp2,scf_in_md%Ent)
        Etemp1 = scf_in_md%Ent 
        Etemp2 = scf_in_md%Emp2 
    end if 
    call cpu_time(scf_end_time)
    scf_time = scf_end_time - scf_start_time
    ! start gradient
    write(*,"(1x,a)") 'Start calculate the first gradinet...'
    call cpu_time(grad_start_time)
    if (trim(mole_in_md%grad_job)=='analy') then 
        call cint_in_md%cal_gradint()
        call grad_in_md%cal_analy_grad(cint_in_md,bas_in_md,scf_in_md,mole_in_md)
    else
        call grad_in_md%cal_num_grad(cint_in_md,scf_in_md,mole_in_md)
    end if
    call cpu_time(grad_end_time)
    grad_time = grad_end_time - grad_start_time
    mole_in_md%guess = 'read'
    ! cal the init_temp
    call this%cal_kin_pot(cint_in_md,scf_in_md,this%vel_old)
    call this%cal_temp(cint_in_md)
    write(*,"(/,1x,a,f8.2,a)") 'Temperature: ',this%temp_cal,' K.'
    if (mole_in_md%nbond_angle /= 0) then
        call this%calc_angle(mole_in_md,scf_in_md)
    end if 
    if (mole_in_md%nbond_angle /= 0 .or. mole_in_md%nbond_length /= 0) then
        open(15,file='mole_geom.dat',status='replace')
        call this%pop_out(mole_in_md,scf_in_md,0)    
    end if 
    open(13,file='pos.xyz',status='replace')
    call this%xyz_out(mole_in_md,0)
    call out_md(this%nstep,this%dt*this%nstep/fs2au,this%job,this%init_temp,'Berendsen',this%bath_temp,this%con_time/fs2au)
    write(99,"(A8,3x,A16,3x,A16,3x,A16,3x,A8,3x,A10,3x,A10)") 'time(fs)','  Etot(Hartree) ','  Ekin(Hartree) ',&
                                                              '  Epot(Hartree) ','  T(K)  ',' t_scf(s) ',' t_grad(s)'
    write(99,"(A8,3x,A16,3x,A16,3x,A16,3x,A8,3x,A10,3x,A10)") '========','================','================',&
                                                              '================','========','==========','=========='
    if (trim(mole_in_md%QC_method)=='hf') then
        scf_in_md%Ent = Etemp1
    else if (trim(mole_in_md%QC_method)=='mp2(fc)' .or. trim(mole_in_md%QC_method)=='mp2(full)') then 
        scf_in_md%Ent = Etemp1
        scf_in_md%Emp2 = Etemp2
    end if 
    call this%E_out(scf_in_md,mole_in_md,0)
    write(*,"(/,1x,a)") 'Position in Born: '
    write(*,'(/,2x,a)') ' Atom           X                 Y                 Z       '
    write(*,'(2x,a)')   '------  ----------------  ----------------  ----------------'
    do i=1,cint_in_md%natom 
        write(*,"(4x,a2,5x,F15.10,3x,F15.10,3x,F15.10)") mole_in_md%elem_list(i),this%pos_old(i,1)&
                                                        ,this%pos_old(i,2),this%pos_old(i,3)
    end do 
    write(*,"(/,1x,a)") 'Velocity: '
    write(*,'(/,2x,a)') ' Atom           X                 Y                 Z       '
    write(*,'(2x,a)')   '------  ----------------  ----------------  ----------------'
    do i=1,cint_in_md%natom 
        write(*,"(4x,a2,5x,F15.10,3x,F15.10,3x,F15.10)") mole_in_md%elem_list(i),this%vel_old(i,1)&
                                                        ,this%vel_old(i,2),this%vel_old(i,3)
    end do 
    ! start md
    do n = 1,this%nstep
        write(*,"(/,a,i5,1x,a,/)") '===== ',n,'AIMD-step  ====='
        do i=1,cint_in_md%natom
            this%acc_old(i,:) = grad_in_md%force(i,:)/(cint_in_md%mass(i)*dalton2au)
            if (mole_in_md%md_job=='Velocity_Verlet') then
                this%pos_new(i,:) = this%pos_old(i,:) + this%vel_old(i,:)*this%dt &
                                    + 0.5_f8*this%acc_old(i,:)*(this%dt**2.0)
            else if(mole_in_md%md_job == 'leapfrog') then 
                this%vel_half(i,:) = this%vel_old(i,:)+0.5*this%dt*this%acc_old(i,:)
                this%pos_new(i,:) = this%pos_old(i,:) + this%vel_half(i,:)*this%dt
            end if 
            mole_in_md%x_list(i) = this%pos_new(i,1)*0.529177249_f8
            mole_in_md%y_list(i) = this%pos_new(i,2)*0.529177249_f8
            mole_in_md%z_list(i) = this%pos_new(i,3)*0.529177249_f8
        end do 
        call this%xyz_out(mole_in_md,n)
        call cint_in_md%set_mole_bas_info(bas_in_md,mole_in_md)
        call cint_in_md%set_cint_info(mole_in_md,1)
        call cpu_time(scf_start_time)
        call cint_in_md%cal_cint()
        write(*,"(1x,a)") 'Start the ab initio SCF...'
        call cpu_time(scf_start_time)
        if (trim(mole_in_md%QC_method)=='hf') then
            call cint_in_md%cal_cint()
            call scf_in_md%scf_pro(cint_in_md,mole_in_md)
            call energy_out(scf_in_md%Ent,scf_in_md%iter_cong)
            Etemp1 = scf_in_md%Ent
        else if (trim(mole_in_md%QC_method)=='mp2(fc)' .or. trim(mole_in_md%QC_method)=='mp2(full)') then 
            call cint_in_md%cal_cint()
            call scf_in_md%scf_pro(cint_in_md,mole_in_md)
            call scf_in_md%init_mp2(cint_in_md)
            call scf_in_md%cal_mp2(cint_in_md)
            call out_mp2(scf_in_md%Emp2,scf_in_md%Ent)
            Etemp1 = scf_in_md%Ent 
            Etemp2 = scf_in_md%Emp2
        end if 
        call cpu_time(scf_end_time)
        scf_time = scf_end_time - scf_start_time
        call cpu_time(grad_start_time)
        write(*,"(1x,a)") 'Start calculate the first gradinet...'
        if (trim(mole_in_md%grad_job)=='analy') then 
            call cint_in_md%cal_gradint()
            call grad_in_md%cal_analy_grad(cint_in_md,bas_in_md,scf_in_md,mole_in_md)
        else
            call grad_in_md%cal_num_grad(cint_in_md,scf_in_md,mole_in_md)
        end if 
        call cpu_time(grad_end_time)
        grad_time = grad_end_time - grad_start_time
        do i=1,cint_in_md%natom
            this%acc_new(i,:) = grad_in_md%force(i,:)/(cint_in_md%mass(i)*dalton2au)
            if (mole_in_md%md_job=='Velocity_Verlet') then
                this%vel_new(i,:) = this%vel_old(i,:) &
                                    + 0.5_f8*(this%acc_new(i,:)+this%acc_old(i,:))*this%dt 
            else if(mole_in_md%md_job == 'leapfrog') then
                this%vel_new(i,:) = this%vel_half(i,:) &
                                    + 0.5_f8*this%acc_new(i,:)*this%dt    
            end if             
        end do
        mv = 0.0_f8
        vel_cent = 0.0_f8 
        mtot = 0.0_f8
        do i=1,cint_in_md%natom
            mv(:) = mv(:) + this%vel_new(i,:)*(cint_in_md%mass(i)*dalton2au)
            mtot = mtot + cint_in_md%mass(i)*dalton2au
        end do 
        vel_cent(:) = mv(:)/mtot
        do i=1,cint_in_md%natom
            this%vel_new(i,:) = this%vel_new(i,:) - vel_cent(:)
        end do
        call this%cal_kin_pot(cint_in_md,scf_in_md,this%vel_new)
        call this%berendsen(cint_in_md)
        write(*,"(/,1x,a,f8.2,a)") 'Temperature: ',this%temp_cal,' K.'
        write(*,"(/,1x,a)") 'Position in Born: '
        write(*,'(/,2x,a)') ' Atom           X                 Y                 Z       '
        write(*,'(2x,a)')   '------  ----------------  ----------------  ----------------'
            do i=1,cint_in_md%natom 
            write(*,"(4x,a2,5x,F15.10,3x,F15.10,3x,F15.10)") mole_in_md%elem_list(i),this%pos_new(i,1)&
                                                            ,this%pos_new(i,2),this%pos_new(i,3)
        end do 
        write(*,"(/,1x,a)") 'Velocity: '
        write(*,'(/,2x,a)') ' Atom           X                 Y                 Z       '
        write(*,'(2x,a)')   '------  ----------------  ----------------  ----------------'
        do i=1,cint_in_md%natom 
            write(*,"(4x,a2,5x,F15.10,3x,F15.10,3x,F15.10)") mole_in_md%elem_list(i),this%vel_new(i,1)&
                                                            ,this%vel_new(i,2),this%vel_new(i,3)
        end do 
        if (trim(mole_in_md%QC_method)=='hf') then
            scf_in_md%Ent = Etemp1
        else if (trim(mole_in_md%QC_method)=='mp2(fc)' .or. trim(mole_in_md%QC_method)=='mp2(full)') then 
            scf_in_md%Ent = Etemp1
            scf_in_md%Emp2 = Etemp2
        end if 
        call this%E_out(scf_in_md,mole_in_md,n)
        this%pos_old = this%pos_new
        this%vel_old = this%vel_new
        if (mole_in_md%nbond_angle /= 0) then 
            call this%calc_angle(mole_in_md,scf_in_md)
        end if 
        if (mole_in_md%nbond_angle /= 0 .or. mole_in_md%nbond_length /= 0) then
            call this%pop_out(mole_in_md,scf_in_md,n)    
        end if 
    end do
    write(99,'(/,1x,a)') 'MD step has done!'
    write(*,'(/,1x,a)') 'MD step has done!'
end subroutine

subroutine xyz_out_sub(this,mole_in_md,nstep)
implicit none
    class(md) :: this
    type(mole_input_info),intent(in) :: mole_in_md
    integer,intent(in) :: nstep
    integer :: i

    write(13,"(I2)") mole_in_md%natom 
    if (trim(this%job)=='Velocity_Verlet') then 
        write(13,"(A9,I4,A9,F7.2,1x,A26)") '#Nstep = ',nstep,', time = ',this%dt*nstep/fs2au,'fs, using Velocity Verlet.'
    else if (trim(this%job)=='leapfrog') then
        write(13,"(A9,I4,A9,F7.2,1x,A26)") '#Nstep = ',nstep,', time = ',this%dt*nstep/fs2au,'fs, using leapfrog.'
    end if
    do i=1,mole_in_md%natom
        write(13,"(A2,3F15.10)") mole_in_md%elem_list(i),mole_in_md%x_list(i),mole_in_md%y_list(i),mole_in_md%z_list(i)
    end do
end subroutine

subroutine E_out_sub(this,scf_in_md,mole_in_md,nstep)
implicit none
    class(md) :: this
    type(scf_info),intent(in) :: scf_in_md
    type(mole_input_info),intent(in) :: mole_in_md
    integer,intent(in) :: nstep
    integer :: i
    if (trim(mole_in_md%QC_method)=='hf') then
    write(99,"(F8.2,3x,F16.10,3x,F16.10,3x,F16.10,3x,F8.2,3x,F10.3,3x,F10.3)") this%dt*nstep/fs2au,scf_in_md%Ent,&
                                                                               this%kin,this%pot,this%temp_cal,scf_time,grad_time
    else if (trim(mole_in_md%QC_method)=='mp2(fc)' .or. trim(mole_in_md%QC_method)=='mp2(full)') then
    write(99,"(F8.2,3x,F16.10,3x,F16.10,3x,F16.10,3x,F8.2,3x,F10.3,3x,F10.3)") this%dt*nstep/fs2au,scf_in_md%Ent+scf_in_md%Emp2,&
                                                                               this%kin,this%pot,this%temp_cal,scf_time,grad_time
    end if
end subroutine

subroutine calc_angle_sub(this,mole_in_md,scf_in_md)
implicit none
    class(md) :: this
    type(scf_info),intent(in) :: scf_in_md
    type(mole_input_info),intent(in) :: mole_in_md
    real(kind=f8),allocatable,dimension(:,:) :: ex
    real(kind=f8),allocatable,dimension(:,:) :: ey
    real(kind=f8),allocatable,dimension(:,:) :: ez
    real(kind=f8),allocatable,dimension(:,:,:) :: costheta
    integer :: i,j,k 
    integer :: status

    allocate (costheta(mole_in_md%natom,mole_in_md%natom,mole_in_md%natom),stat=status)
    allocate (ex(mole_in_md%natom,mole_in_md%natom),stat=status)
    allocate (ey(mole_in_md%natom,mole_in_md%natom),stat=status)
    allocate (ez(mole_in_md%natom,mole_in_md%natom),stat=status)
    
    this%bond_length = scf_in_md%Bond_R
    do j=1,mole_in_md%natom
        do k=1,mole_in_md%natom
            ex(j,k)=-((mole_in_md%x_list(j)-mole_in_md%x_list(k))/this%bond_length(j,k))
            ey(j,k)=-((mole_in_md%y_list(j)-mole_in_md%y_list(k))/this%bond_length(j,k)) 
            ez(j,k)=-((mole_in_md%z_list(j)-mole_in_md%z_list(k))/this%bond_length(j,k)) 
        end do
    end do 
    do i=1,mole_in_md%natom
        do j=1,mole_in_md%natom
            do k=1,mole_in_md%natom
                costheta(i,j,k)=(ex(j,i)*ex(j,k))+(ey(j,i)*ey(j,k))+(ez(j,i)*ez(j,k))
                this%bond_angle(i,j,k)=(acos(costheta(i,j,k))/pi)*180.0_f8
            end do
        end do
    end do
end subroutine

subroutine pop_out_sub(this,mole_in_md,scf_in_md,nstep)
implicit none
    class(md) :: this
    type(scf_info),intent(in) :: scf_in_md
    type(mole_input_info),intent(in) :: mole_in_md
    integer,intent(in) :: nstep
    integer :: i,j

    write(15,'(A)') '========================================================================'
    if (trim(this%job)=='Velocity_Verlet') then 
        write(15,"(A9,I4,A9,F7.2,1x,A26)") '#Nstep = ',nstep,', time = ',this%dt*nstep/fs2au,'fs, using Velocity Verlet.'
    else if (trim(this%job)=='leapfrog') then
        write(15,"(A9,I4,A9,F7.2,1x,A26)") '#Nstep = ',nstep,', time = ',this%dt*nstep/fs2au,'fs, using leapfrog.'
    end if
    write(15,'(A)') '========================================================================'
    if (mole_in_md%nbond_length /= 0) then 
        write(15,"(A,I2,2x,A)") 'You have chosen ',mole_in_md%nbond_length,'bond length groups, here they are(in Angstroms): '
    
        j = 1
        do i =1, mole_in_md%nbond_length
            write(15,'(3x,A,I2,A,I2,A,I2,2x,F8.5)') 'Bond length group',i,' :' &
                                                    ,mole_in_md%bond_length_index(j)&
                                                    ,'-',mole_in_md%bond_length_index(j+1) &
                                                    ,scf_in_md%Bond_R(mole_in_md%bond_length_index(j)&
                                                    ,mole_in_md%bond_length_index(j+1))
            j = j+2
        end do 
    else 
        write(15,"(A)") 'There is no bond length.'
    end if 
    if (mole_in_md%nbond_angle /= 0) then 
        write(15,"(A,I2,2x,A)") 'You have chosen ',mole_in_md%nbond_angle,'bond angle groups, here they are(in degree): '
        j = 1
        do i =1, mole_in_md%nbond_angle
            write(15,'(3x,A,I2,A,I2,A,I2,A,I2,2x,F10.4)') 'Bond angle group',i,' :' &
                                                        ,mole_in_md%bond_angle_index(j)&
                                                        ,'-',mole_in_md%bond_angle_index(j+1),'-'&
                                                        ,mole_in_md%bond_angle_index(j+2) &
                                                        ,this%bond_angle(mole_in_md%bond_angle_index(j)&
                                                        ,mole_in_md%bond_angle_index(j+1),mole_in_md%bond_angle_index(j+2))
            j = j+3
        end do 
    else 
        write(15,"(A)") 'There is no bond angle.'
    end if 
    write(15,*) ' '
end subroutine

end module







