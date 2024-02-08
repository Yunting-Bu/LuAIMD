module basis
use machina_basic, only: i4, f8
use read_input,only : loclabel
implicit none

type,public :: atom_basis_info
    integer :: atom_charge
    integer :: atom_nbas
    integer :: atom_nshl
    integer :: atom_nprm
    integer,allocatable,dimension(:) :: bas_angl
	integer,allocatable,dimension(:) :: prm_num
    real(kind=f8),allocatable,dimension(:) :: expnt
    real(kind=f8),allocatable,dimension(:) :: coeff
	real(kind=f8) :: atom_mass
	real(kind=f8) :: atom_mass_iso
	integer :: atom_fz_core

    contains
        procedure,public :: set_basis => set_basis_sub
    
end type atom_basis_info

contains


!!-------- Get the infomation of atom_basis
    subroutine set_basis_sub(this,basis_name,elem_label)
    implicit none 
        class(atom_basis_info) :: this
        character(len=*),intent(in) :: basis_name
		character(len=*),intent(in) :: elem_label
		character(len=200) label
		character(len=2) c2tmp
		character(len=2) elem
		character(len=200) ctmp,ctmpint,ctmpreal
		character(len=3) anlg_label
		real(kind=f8) :: creal,creal2int
		integer :: n,status,i,j,k,l,cint,count,nbas_c,nprm_c

    open(unit=10,file='basis/'//trim(basis_name),status='old')
!locate the atom
	rewind(10)
	call loclabel(10,trim(elem_label)//'     0')
!get the elem_name
	read (10,*) elem,cint
	i = 0
	j = 0
	count = 0
	nbas_c = 0
	nprm_c = 0
	do 
		cint = 0
		creal = 0.0
		read (10,*) c2tmp
		if (index(trim(c2tmp),'*') /= 0) exit
		if (index(trim(c2tmp),'S') /= 0 .or. index(trim(c2tmp),'P') /= 0 .or. &
			index(trim(c2tmp),'D') /= 0 .or. index(trim(c2tmp),'F') /= 0) then
			backspace(10)
			read (10,*) ctmp,ctmpint,ctmpreal
			read(ctmpint,'(i2)') cint 
			read(ctmpreal,'(f3.2)') creal
			anlg_label = trim(ctmp)
			select case (trim(anlg_label))
        			case ('S'); nbas_c = nbas_c + 1
        			case ('P'); nbas_c = nbas_c + 3
        			case ('D'); nbas_c = nbas_c + 5
					case ('F'); nbas_c = nbas_c + 7
					case ('SP'); nbas_c = nbas_c + 4
        	end select
			if (trim(anlg_label)=='SP') then 
				count = count + 2
				nprm_c = nprm_c + cint*2
			else 
				count = count + 1
				nprm_c = nprm_c + cint		
			end if 	
		end if 
	end do 
	this%atom_nshl = count
	this%atom_nbas = nbas_c
	this%atom_nprm = nprm_c
	allocate(this%bas_angl(this%atom_nshl),stat=status)
	allocate(this%prm_num(this%atom_nshl),stat=status)
	allocate(this%expnt(this%atom_nprm),stat=status)
	allocate(this%coeff(this%atom_nprm),stat=status)

	rewind(10)
	call loclabel(10,trim(elem_label)//'     0')
	read (10,*) elem,cint
	do 
		cint = 0
		creal = 0.0
		read (10,*) c2tmp
		if (index(trim(c2tmp),'*') /= 0) exit
		if (index(trim(c2tmp),'S') /= 0 .or. index(trim(c2tmp),'P') /= 0 .or. &
			index(trim(c2tmp),'D') /= 0 .or. index(trim(c2tmp),'F') /= 0) then
			i = i+1
			backspace(10)
			read (10,*) anlg_label,ctmpint,ctmpreal
			read(ctmpint,'(i2)') this%prm_num(i)
				select case (trim(anlg_label))
        			case ('S'); this%bas_angl(i) = 0
        			case ('P'); this%bas_angl(i) = 1
        			case ('D'); this%bas_angl(i) = 2
					case ('F'); this%bas_angl(i) = 3
					case ('SP')
						this%bas_angl(i) = 0
						this%bas_angl(i+1) = 1
						this%prm_num(i+1) = this%prm_num(i)
						i = i+1
        		end select
		else 
			j = j+1
			if (trim(anlg_label) == 'SP') then 
				backspace(10)
				do k = j,j+this%prm_num(i)-1
					read(10,*) this%expnt(k),this%coeff(k),creal
				end do 
				do l = 1, this%prm_num(i)
					backspace(10)
				end do 
				j = j + this%prm_num(i)
				do k = j,j+this%prm_num(i)-1
					read(10,*) this%expnt(k),creal,this%coeff(k)
				end do 
				j = j+ this%prm_num(i)-1
			else 
				backspace(10)
				read (10,*) this%expnt(j),this%coeff(j)
			end if
		end if 
	end do 
	close(10)

	open(unit=9,file='atom.dat',status='old')
	rewind(9)
	call loclabel(9,trim(elem_label))
	read(9,*) ctmp
	read(9,*) this%atom_charge
	read(9,*) this%atom_mass
	read(9,*) this%atom_mass_iso
	read(9,*) this%atom_fz_core
	close(9)

    end subroutine

end module basis