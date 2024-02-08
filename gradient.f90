module gradient
use machina_basic, only: i4, f8
use basis 
use read_input
use cint
use scf
implicit none 

type(cint_info) :: cint_in_grad
type(atom_basis_info),dimension(:),allocatable :: bas_in_grad
type(scf_info) :: scf_in_grad
type(mole_input_info) :: mole_in_grad

type,public :: grad_info
    real(kind=f8),allocatable,dimension(:,:,:,:) :: dS 
    real(kind=f8),allocatable,dimension(:,:,:,:) :: dH
    real(kind=f8),allocatable,dimension(:,:,:,:,:,:) :: deri
    real(kind=f8),allocatable,dimension(:,:) :: dVnn
    real(kind=f8),allocatable,dimension(:,:) :: grad
    real(kind=f8),allocatable,dimension(:,:) :: force

    contains
        procedure,public :: cal_analy_grad => cal_analy_grad_sub
        procedure,public :: cal_num_grad => cal_num_grad_sub

end type 

contains

subroutine cal_analy_grad_sub(this,cint_in_grad,bas_in_grad,scf_in_grad,mole_in_grad)
implicit none 
    class(grad_info) :: this
    type(atom_basis_info),dimension(:),allocatable,intent(in) :: bas_in_grad
    type(cint_info),intent(inout) :: cint_in_grad 
    type(scf_info),intent(in) :: scf_in_grad 
    type(mole_input_info),intent(in) :: mole_in_grad
    real(kind=f8),allocatable,dimension(:,:,:,:) :: temp
    real(kind=f8),allocatable,dimension(:,:,:,:) :: dS_mo
    real(kind=f8),allocatable,dimension(:,:) :: F_mo
    real(kind=f8),allocatable,dimension(:,:,:,:,:,:) :: temp1
    integer :: status,i,j,k,l,ibas,iatm,ebas
    real(kind=f8) :: tempVnn
    real(kind=f8),allocatable,dimension(:,:) :: part1
    real(kind=f8),allocatable,dimension(:,:) :: part2
    real(kind=f8),allocatable,dimension(:,:) :: part3

    allocate(this%dS(cint_in_grad%natom,cint_in_grad%nbas,cint_in_grad%nbas,3),stat=status)
    allocate(this%dH(cint_in_grad%natom,cint_in_grad%nbas,cint_in_grad%nbas,3),stat=status)
    allocate(this%dVnn(cint_in_grad%natom,3),stat=status)
    allocate(this%force(cint_in_grad%natom,3),stat=status)
    allocate(this%grad(cint_in_grad%natom,3),stat=status)
    allocate(part1(cint_in_grad%natom,3),stat=status)
    allocate(part2(cint_in_grad%natom,3),stat=status)
    allocate(part3(cint_in_grad%natom,3),stat=status)
    allocate(F_mo(cint_in_grad%nbas,cint_in_grad%nbas),stat=status)
    allocate(this%deri(cint_in_grad%natom,cint_in_grad%nbas,cint_in_grad%nbas,cint_in_grad%nbas,cint_in_grad%nbas,3),stat=status)
    allocate(temp(cint_in_grad%natom,cint_in_grad%nbas,cint_in_grad%nbas,3),stat=status)
    allocate(dS_mo(cint_in_grad%natom,cint_in_grad%nbas,cint_in_grad%nbas,3),stat=status)
    allocate(temp1(cint_in_grad%natom,cint_in_grad%nbas,cint_in_grad%nbas,cint_in_grad%nbas,cint_in_grad%nbas,3),stat=status)
    this%dS = 0.0_f8
    this%dH = 0.0_f8
    this%deri = 0.0_f8
    temp = 0.0_f8
    temp1 = 0.0_f8
    ! dS    
    do iatm =1,cint_in_grad%natom
        ibas = cint_in_grad%atom_bas_index(iatm)
        ebas = ibas + bas_in_grad(iatm)%atom_nbas - 1
        temp(iatm,ibas:ebas,:,:) = - cint_in_grad%ipovlp(ibas:ebas,:,:)
    end do 
    do i=1,cint_in_grad%nbas 
        do j=1,cint_in_grad%nbas
            this%dS(:,i,j,:) = this%dS(:,i,j,:) + temp(:,i,j,:) + temp(:,j,i,:)
        end do 
    end do 
    !deri
    do iatm =1,cint_in_grad%natom
        ibas = cint_in_grad%atom_bas_index(iatm)
        ebas = ibas + bas_in_grad(iatm)%atom_nbas - 1
        temp1(iatm,ibas:ebas,:,:,:,:) = - cint_in_grad%iperi(ibas:ebas,:,:,:,:)   
    end do
    do i=1,cint_in_grad%nbas 
        do j=1,cint_in_grad%nbas
            do k=1,cint_in_grad%nbas
                do l=1,cint_in_grad%nbas                                       
                    this%deri(:,i,j,k,l,:) = this%deri(:,i,j,k,l,:) &
                                                + temp1(:,i,j,k,l,:) + temp1(:,j,i,k,l,:)&
                                                + temp1(:,k,l,i,j,:) + temp1(:,l,k,i,j,:)
                end do
            end do
        end do 
    end do 
!            write(*,*) '...'
!do i=1,cint_in_grad%nbas
!    write(*,"(*(F15.10))") cint_in_grad%ipovlp(i,:,1)
!end do 
!write(*,*) '...'
!do i=1,cint_in_grad%nbas
!    write(*,"(*(F15.10))") cint_in_grad%ipovlp(i,:,2)
!end do 
!write(*,*) '...'
!do i=1,cint_in_grad%nbas
!    write(*,"(*(F15.10))") cint_in_grad%ipovlp(i,:,3)
!end do 
    !dH
    temp = 0.0_f8
    call cint_in_grad%set_cint_info(mole_in_grad,2)
    do iatm =1,cint_in_grad%natom
        ibas = cint_in_grad%atom_bas_index(iatm)
        ebas = ibas + bas_in_grad(iatm)%atom_nbas - 1
        temp(iatm,ibas:ebas,:,:) = temp(iatm,ibas:ebas,:,:) - cint_in_grad%ipkin(ibas:ebas,:,:)
        temp(iatm,ibas:ebas,:,:) = temp(iatm,ibas:ebas,:,:) - cint_in_grad%ipnuc(ibas:ebas,:,:)
        call cint_in_grad%rinv_at_nucleus(iatm)
        call cint_in_grad%cal_iprinv()
        cint_in_grad%env(5:7) = 0.0_f8
        temp(iatm,:,:,:) = temp(iatm,:,:,:) - cint_in_grad%iprinv(:,:,:) * cint_in_grad%charge(iatm)
    end do 
    do i=1,cint_in_grad%nbas 
        do j=1,cint_in_grad%nbas
            this%dH(:,i,j,:) = this%dH(:,i,j,:) + temp(:,i,j,:) + temp(:,j,i,:)
        end do 
    end do
    !dVnn
    this%dVnn = 0.0_f8
    do i =1,cint_in_grad%natom
        do j =1,cint_in_grad%natom
            if (i == j) then
            tempVnn = 0.0_f8
            else
            tempVnn = cint_in_grad%charge(i)*cint_in_grad%charge(j)*(mole_in_grad%x_list(j)-mole_in_grad%x_list(i))/0.529177249_f8 &
                      /((scf_in_grad%Bond_R(i,j)/0.529177249_f8)**3.0_f8)
            this%dVnn(i,1) = tempVnn + this%dVnn(i,1)
            tempVnn = cint_in_grad%charge(i)*cint_in_grad%charge(j)*(mole_in_grad%y_list(j)-mole_in_grad%y_list(i))/0.529177249_f8 &
                      /((scf_in_grad%Bond_R(i,j)/0.529177249_f8)**3.0_f8)
            this%dVnn(i,2) = tempVnn + this%dVnn(i,2)
            tempVnn = cint_in_grad%charge(i)*cint_in_grad%charge(j)*(mole_in_grad%z_list(j)-mole_in_grad%z_list(i))/0.529177249_f8 &
                      /((scf_in_grad%Bond_R(i,j)/0.529177249_f8)**3.0_f8)
            this%dVnn(i,3) = tempVnn + this%dVnn(i,3)
            end if
        end do 
    end do 
    part1 = 0.0_f8
        do i=1,cint_in_grad%nbas 
            do j=1,cint_in_grad%nbas
                part1(:,1) = part1(:,1) + scf_in_grad%P(i,j)*this%dH(:,i,j,1)
                part1(:,2) = part1(:,2) + scf_in_grad%P(i,j)*this%dH(:,i,j,2)
                part1(:,3) = part1(:,3) + scf_in_grad%P(i,j)*this%dH(:,i,j,3)
            end do 
        end do 
    part2 = 0.0_f8
        do i=1,cint_in_grad%nbas 
            do j=1,cint_in_grad%nbas
                do k=1,cint_in_grad%nbas
                    do l=1,cint_in_grad%nbas
                        part2(:,1) = part2(:,1) &
                                        + (0.5_f8*scf_in_grad%P(i,j)*scf_in_grad%P(k,l) &
                                        -0.25_f8*scf_in_grad%P(i,k)*scf_in_grad%P(j,l)) &
                                        * this%deri(:,i,j,k,l,1)              
                        part2(:,2) = part2(:,2) &
                                        + (0.5_f8*scf_in_grad%P(i,j)*scf_in_grad%P(k,l) &
                                        -0.25_f8*scf_in_grad%P(i,k)*scf_in_grad%P(j,l)) &
                                        * this%deri(:,i,j,k,l,2)   
                        part2(:,3) = part2(:,3) &
                                        + (0.5_f8*scf_in_grad%P(i,j)*scf_in_grad%P(k,l) &
                                        -0.25_f8*scf_in_grad%P(i,k)*scf_in_grad%P(j,l)) &
                                        * this%deri(:,i,j,k,l,3)   
                    end do
                end do 
            end do 
        end do 
    part3 = 0.0_f8
    do i=1,cint_in_grad%nbas
        do j=1,cint_in_grad%nbas
            F_mo(i,j) = 0.0_f8
            do k=1,cint_in_grad%nbas
                do l=1,cint_in_grad%nbas
                    F_mo(i,j)=F_mo(i,j)+scf_in_grad%C(k,i)*scf_in_grad%F(k,l)*scf_in_grad%C(l,j)
                end do
            end do
        end do
    end do
    do i=1,cint_in_grad%nbas
        do j=1,cint_in_grad%nbas
            dS_mo(:,i,j,:) = 0.0_f8
            do k=1,cint_in_grad%nbas
                do l=1,cint_in_grad%nbas
                    dS_mo(:,i,j,:)=dS_mo(:,i,j,:)+scf_in_grad%C(k,i)*this%dS(:,k,l,:)*scf_in_grad%C(l,j)
                end do
            end do
        end do
    end do
    do i=1,scf_in_grad%nocc
        do j=1,scf_in_grad%nocc
            part3(:,1) = 2.0_f8 * F_mo(i,j) * dS_mo(:,i,j,1) + part3(:,1)
            part3(:,2) = 2.0_f8 * F_mo(i,j) * dS_mo(:,i,j,2) + part3(:,2)
            part3(:,3) = 2.0_f8 * F_mo(i,j) * dS_mo(:,i,j,3) + part3(:,3)
        end do
    end do
    this%force = 0.0_f8
    do iatm =1,cint_in_grad%natom
        this%grad(iatm,1) = part1(iatm,1) + part2(iatm,1) - part3(iatm,1) + this%dVnn(iatm,1)
        this%grad(iatm,2) = part1(iatm,2) + part2(iatm,2) - part3(iatm,2) + this%dVnn(iatm,2)
        this%grad(iatm,3) = part1(iatm,3) + part2(iatm,3) - part3(iatm,3) + this%dVnn(iatm,3)
    end do
write(*,'(/,1x,a,/)') '=== Gradient ==='
    write(*,'(/,2x,a)') ' Atom           X                   Y                   Z       '
    write(*,'(2x,a)')   '------  ----------------    ----------------    ----------------'
do i=1,cint_in_grad%natom
   write(*,"(4x,a2,*(5x,F15.10))") mole_in_grad%elem_list(i),this%grad(i,:)
end do
this%force = - this%grad
deallocate(temp)
deallocate(dS_mo)
deallocate(F_mo)
deallocate(temp1)
deallocate(part1)
deallocate(part2)
deallocate(part3)
end subroutine

subroutine cal_num_grad_sub(this,cint_in_grad,scf_in_grad,mole_in_grad)
implicit none 
    class(grad_info) :: this
    type(cint_info),intent(inout) :: cint_in_grad 
    type(scf_info),intent(inout) :: scf_in_grad 
    type(mole_input_info),intent(inout) :: mole_in_grad
    real(kind=f8) :: h = 1.0D-4
    real(kind=f8) :: e_negh
    real(kind=f8) :: e_posh
    real(kind=f8),allocatable,dimension(:,:) :: coor_temp
    real(kind=f8),allocatable,dimension(:,:) :: coor_h
    integer :: i,j,status

    allocate(this%force(cint_in_grad%natom,3),stat=status)
    allocate(this%grad(cint_in_grad%natom,3),stat=status)
    allocate(coor_h(cint_in_grad%natom,3),stat=status)
    allocate(coor_temp(cint_in_grad%natom,3),stat=status)

    coor_temp(:,1) = mole_in_grad%x_list
    coor_temp(:,2) = mole_in_grad%y_list
    coor_temp(:,3) = mole_in_grad%z_list
    do i = 1,mole_in_grad%natom 
        do j = 1,3
            coor_h(:,1) = mole_in_grad%x_list/0.529177249_f8
            coor_h(:,2) = mole_in_grad%y_list/0.529177249_f8
            coor_h(:,3) = mole_in_grad%z_list/0.529177249_f8
            select case (j)
                case(1); mole_in_grad%x_list(i) = (coor_h(i,1)-h)*0.529177249_f8
                case(2); mole_in_grad%y_list(i) = (coor_h(i,2)-h)*0.529177249_f8   
                case(3); mole_in_grad%z_list(i) = (coor_h(i,3)-h)*0.529177249_f8 
            end select
            call cint_in_grad%set_cint_info(mole_in_grad,1)
            if (trim(mole_in_grad%QC_method)=='hf') then
                call cint_in_grad%cal_cint()
                call scf_in_grad%scf_pro(cint_in_grad,mole_in_grad)
                e_negh = scf_in_grad%Ent 
            else if (trim(mole_in_grad%QC_method)=='mp2(fc)' .or. trim(mole_in_grad%QC_method)=='mp2(full)') then 
                call cint_in_grad%cal_cint()
                call scf_in_grad%scf_pro(cint_in_grad,mole_in_grad)
                call scf_in_grad%init_mp2(cint_in_grad)
                call scf_in_grad%cal_mp2(cint_in_grad)
                e_negh = scf_in_grad%Ent + scf_in_grad%Emp2
            end if 
            mole_in_grad%x_list = coor_temp(:,1)
            mole_in_grad%y_list = coor_temp(:,2)
            mole_in_grad%z_list = coor_temp(:,3)
            coor_h = coor_temp/0.529177249_f8 
            select case (j)
                case(1); mole_in_grad%x_list(i) = (coor_h(i,1)+h)*0.529177249_f8
                case(2); mole_in_grad%y_list(i) = (coor_h(i,2)+h)*0.529177249_f8   
                case(3); mole_in_grad%z_list(i) = (coor_h(i,3)+h)*0.529177249_f8 
            end select
            call cint_in_grad%set_cint_info(mole_in_grad,1)
            if (trim(mole_in_grad%QC_method)=='hf') then
                call cint_in_grad%cal_cint()
                call scf_in_grad%scf_pro(cint_in_grad,mole_in_grad)
                e_posh = scf_in_grad%Ent 
            else if (trim(mole_in_grad%QC_method)=='mp2(fc)' .or. trim(mole_in_grad%QC_method)=='mp2(full)') then 
                call cint_in_grad%cal_cint()
                call scf_in_grad%scf_pro(cint_in_grad,mole_in_grad)
                call scf_in_grad%init_mp2(cint_in_grad)
                call scf_in_grad%cal_mp2(cint_in_grad)
                e_posh = scf_in_grad%Ent + scf_in_grad%Emp2
            end if 
            this%grad(i,j) = (e_posh-e_negh)/(2.0*h)
            mole_in_grad%x_list = coor_temp(:,1)
            mole_in_grad%y_list = coor_temp(:,2)
            mole_in_grad%z_list = coor_temp(:,3)
        end do 
    end do 
    write(*,'(/,1x,a,/)') '=== Gradient ==='
    write(*,'(/,2x,a)') ' Atom           X                   Y                   Z       '
    write(*,'(2x,a)')   '------  ----------------    ----------------    ----------------'
do i=1,mole_in_grad%natom
   write(*,"(4x,a2,*(5x,F15.10))") mole_in_grad%elem_list(i),this%grad(i,:)
end do
    this%force = -this%grad 
    deallocate(coor_temp)
    deallocate(coor_h)

end subroutine

end module