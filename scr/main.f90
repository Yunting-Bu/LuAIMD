program main 
use basis
use read_input
use cint
use scf
use gradient
use md_pro
use output
implicit none 

character(len=7) :: basis_name
character(len=200) :: input_name
type(atom_basis_info),dimension(:),allocatable :: atom_bas0 
type(mole_input_info) :: mole0
type(cint_info) :: cint0
type(scf_info) :: scf0
type(grad_info) :: grad0
type(md) :: md0
real(kind=f8) :: cpu_start_time
real(kind=f8) :: cpu_end_time
real(kind=f8) :: cputime
integer :: hour
integer :: min 
integer :: sec 
integer :: i,j,k,l,status
character(len=100) :: namepre
integer :: str_len


write (*,'(1x,a)') 'Please enter the input file: '
read(*,*) input_name
str_len = len_trim(input_name)
namepre = adjustl(input_name(1:str_len-4))
call cpu_time(cpu_start_time)
write(*,'(/,1x,a,/)')'======================================================================================='
write(*,'(1x,20x,a,20x,/)') 'LuAIMD: An ab initio molecular dynamics program'
write(*,'(1x,32x,a,26x,/)') 'Author: Yun-ting Bu'
write(*,'(1x,a,/)')'======================================================================================='
open(99,file=trim(namepre)//'.out',status='replace')
call out_title()
write(*,"(1x,a)") 'Read the input...'
call mole0%set_input(trim(input_name))
allocate(atom_bas0(mole0%natom),stat=status)

call out_mole_info(mole0%natom,mole0%basis,mole0%mole_charge,&
                   mole0%mole_mul,mole0%elem_list,mole0%x_list,mole0%y_list,mole0%z_list)
write(*,"(1x,a)") 'Setup the basis information...'
do i=1,mole0%natom 
    call atom_bas0(i)%set_basis(mole0%basis,mole0%elem_list(i))
end do 
write(*,"(1x,a)") 'Start the molecular dynamics!'
call md0%md_int(mole0,cint0,atom_bas0,scf0,grad0)
call cpu_time(cpu_end_time)
cputime = cpu_end_time - cpu_start_time
if (cputime >= 3600) then 
    hour = int(cputime)/3600 
    sec = mod(int(cputime),3600)
    if (sec >= 60) then 
        min = sec/60
        sec = mod(sec,60)
    else 
        min = 0
    end if 
elseif (cputime<3600 .and. cputime>=60) then 
    hour = 0
    min = int(cputime)/60
    sec = mod(int(cputime),60)
else 
    hour = 0
    min = 0
    sec = 0
end if 
call out_cputime(hour,min,sec,cputime)
    deallocate(md0%pos_old)
    deallocate(md0%pos_new)
    deallocate(md0%vel_old)
    deallocate(md0%vel_new)
    deallocate(md0%vel_half)
    deallocate(md0%acc_old)
    deallocate(md0%acc_new)
    deallocate(md0%bond_length)
    deallocate(md0%bond_angle)
    
    if(trim(mole0%grad_job)=='analy') then 
        deallocate(grad0%dS)
        deallocate(grad0%dH)
        deallocate(grad0%deri)
        deallocate(grad0%dVnn)
    end if 
    deallocate(grad0%force)
    deallocate(grad0%grad)

    deallocate(scf0%H_GWH)
    deallocate(scf0%Hcore)
    deallocate(scf0%P)
    deallocate(scf0%P_read)
    deallocate(scf0%F)
    deallocate(scf0%C)
    deallocate(scf0%G)
    deallocate(scf0%X)
    deallocate(scf0%E)
    deallocate(scf0%Bond_R)

    deallocate(cint0%cntr_odr)
    deallocate(cint0%angl)
    deallocate(cint0%shl_belong_to_atom)
    deallocate(cint0%charge)
    deallocate(cint0%shl_index)
    deallocate(cint0%atom_bas_index)
    deallocate(cint0%mass)
    deallocate(cint0%coor)
    deallocate(cint0%expnt)
    deallocate(cint0%coeff)
    deallocate(cint0%env)
    deallocate(cint0%atm)
    deallocate(cint0%bas)
    deallocate(cint0%zeta)
    deallocate(cint0%S)
    deallocate(cint0%T)
    deallocate(cint0%V)
    deallocate(cint0%eri)
    if(trim(mole0%grad_job)=='analy') then
        deallocate(cint0%ipovlp)
        deallocate(cint0%ipkin)
        deallocate(cint0%ipnuc)
        deallocate(cint0%iperi)
        deallocate(cint0%iprinv)
    end if 

    do i=1,mole0%natom
        deallocate(atom_bas0(i)%bas_angl)
        deallocate(atom_bas0(i)%prm_num)
        deallocate(atom_bas0(i)%expnt)
        deallocate(atom_bas0(i)%coeff)
    end do 

    deallocate(mole0%elem_list)
    deallocate(mole0%x_list)
    deallocate(mole0%y_list)
    deallocate(mole0%z_list)
    if (mole0%nbond_length /= 0) then 
        deallocate(mole0%bond_length_index)
    end if 
    if (mole0%nbond_angle /= 0) then 
        deallocate(mole0%bond_angle_index)
    end if
    if (mole0%nbond_angle /= 0 .and. mole0%nbond_length /=0) then 
        close(15)
    end if
    close(13)
write(*,"(1x,a)") 'Detailed information will be given in output file, pos.xyz and mole_geom.dat!'
write(*,'(1x,a)')'Normal termination, have a nice day! ( ^ ω ^ )'
close(99)
end program
