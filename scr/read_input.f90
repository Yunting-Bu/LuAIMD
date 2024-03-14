module read_input
use machina_basic, only: i4, f8
implicit none

type,public :: mole_input_info
    character(len=200) :: input_name
    character(len=200) :: xyz_name
    integer :: mole_charge
    integer :: mole_mul
    character(len=2),allocatable,dimension(:) :: elem_list
    integer :: natom
    real(kind=f8),allocatable,dimension(:) :: x_list
    real(kind=f8),allocatable,dimension(:) :: y_list
    real(kind=f8),allocatable,dimension(:) :: z_list
    ! for QC
    real(kind=f8) :: damp
    character(len=200) :: QC_method
    character(len=200) :: basis 
    character(len=200) :: guess
    character(len=200) :: grad_job
    logical :: DIIS  
    ! for md
    character(len=200) :: md_job
    character(len=200) :: temp_job
    real(kind=8) :: init_temp
    real(kind=8) :: bath_temp
    real(kind=8) :: con_time
    real(kind=8) :: dt
    integer :: nstep
    ! for mole_analy
    integer :: nbond_length
    integer :: nbond_angle
    integer,allocatable,dimension(:) :: bond_length_index
    integer,allocatable,dimension(:) :: bond_angle_index
    
    contains
        procedure,public :: set_input => set_input_sub

end type

contains

subroutine set_input_sub(this,input_name)
implicit none 
    class(mole_input_info) :: this
    character(len=*),intent(in) :: input_name
    character(len=200) ctmp
	integer :: status,i

    this%input_name = input_name
    open (11,file=this%input_name,status='old')
!get the QC_ctrl
!$QC_ctrl
!    method:     hf   <== hf/mp2
!    basis:      sto-3g
!    grad:       analy  <== analy(only for hf)/num
!    guess:      Huckel <== core/GWH/Huckel/read, read isn't a key word, it will be used in AIMD
!    damp:       0.5
!    DIIS:       .true.
!$end 
    call loclabel(11,'$QC_ctrl')
    read (11,*) ctmp 
    read (11,*) ctmp,this%QC_method
    read (11,*) ctmp,this%basis 
    read (11,*) ctmp,this%grad_job
    read (11,*) ctmp,this%guess 
    read (11,*) ctmp,this%damp
    read (11,*) ctmp,this%DIIS
!get the MD_ctrl
!$MD_ctrl
!    method:     Velocity_Verlet <== Velocity_Verlet/leapfrog
!    init_temp:  298.15  <== K
!    dt:         0.5  <== fs
!    Nstep:      1000
!    thermostat: Berendsen
!    bath_temp:  298.15
!    con_time:   30.0  <== fs
!$end 
    call loclabel(11,'$MD_ctrl')
    read (11,*) ctmp 
    read (11,*) ctmp,this%md_job
    read (11,*) ctmp,this%init_temp
    read (11,*) ctmp,this%dt
    read (11,*) ctmp,this%nstep
    read (11,*) ctmp,this%temp_job
    read (11,*) ctmp,this%bath_temp
    read (11,*) ctmp,this%con_time

!get the MOLE_analy job, for bond_length, bond_angle and bond_mayer
!$MOLE_geom
!    number_of_bond_length: nbond_lenght
!    list: 1,4 4,5 <= such as
!    number_of_bond_angle: nbond_angle
!    list: 4,5,6
!$end 
! if nbond_length/angle = 0, don't calculate them
    call loclabel(11,'$MOLE_analy')
    read (11,*) ctmp 
    ! read bond_length
    read (11,*) ctmp,this%nbond_length 
    if (this%nbond_length == 0) then 
        read(11,*) ctmp 
    else 
        allocate(this%bond_length_index(2*this%nbond_length))
        read (11,*) ctmp, this%bond_length_index
    end if 
    ! read bond_angle
    read (11,*) ctmp,this%nbond_angle 
    if (this%nbond_angle == 0) then 
        read(11,*) ctmp 
    else 
        allocate(this%bond_angle_index(3*this%nbond_angle))
    end if 
    read (11,*) ctmp, this%bond_angle_index
!get the charge, multiplicity and geom
!$geom
!    0 1 <= charge and multiplicity
!    H2O.xyz
!$end
    call loclabel(11,'$geom')
    read (11,*) ctmp
    read (11,*) this%mole_charge,this%mole_mul
    read (11,*) this%xyz_name
    close (11)
!get the elements, atom number and xyz
    open (12,file=trim(this%xyz_name),status='old')
    read (12,*) this%natom
    allocate(this%elem_list(this%natom),stat=status)
    allocate(this%x_list(this%natom),stat=status)
    allocate(this%y_list(this%natom),stat=status)
    allocate(this%z_list(this%natom),stat=status)
    read (12,*) ctmp
    backspace(12)
    do i=1,this%natom
    read (12,*) this%elem_list(i),this%x_list(i),this%y_list(i),this%z_list(i)
    this%x_list(i) = this%x_list(i)
    this%y_list(i) = this%y_list(i)
    this%z_list(i) = this%z_list(i)
    end do
    close(12)
end subroutine



! From Multiwfn
!!-------- Locate the line where the label first appears in fileid
!Return ifound=1 if found the label, else return 0
!Default is rewind, if irewind=0 then will not rewind
!If the current line just has the label, calling this subroutine will do nothing
!maxline define the maximum number of lines that will be searched, default is search the whole file
subroutine loclabel(fileid,label,ifound,irewind,maxline)
    integer fileid,ierror,iline
    integer,optional :: ifound,irewind,maxline
    character c200tmp*200
    character(len=*) label
    if (.not.present(irewind)) then
	    rewind(fileid)
    else
	    if (irewind==1) rewind(fileid)
    end if
    if (.not.present(maxline)) then
	    do while(.true.)
		    read(fileid,"(a)",iostat=ierror) c200tmp
		    if (ierror/=0) exit
		    if (index(c200tmp,label)/=0) then
			    backspace(fileid)
			    if (present(ifound)) ifound=1 !Found result
			    return
		    end if
	    end do
    else
	    do iline=1,maxline
		    read(fileid,"(a)",iostat=ierror) c200tmp
		    if (ierror/=0) exit
		    if (index(c200tmp,label)/=0) then
			    backspace(fileid)
			    if (present(ifound)) ifound=1 !Found result
			    return
		    end if
	    end do
    end if
    if (present(ifound)) ifound=0
end subroutine

end module