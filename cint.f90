module cint
use machina_basic
use basis 
use read_input
implicit none 
!the parameters that libcint needs
    integer, parameter :: CHARGE_OF = 1
    integer, parameter :: PTR_COORD = 2
    integer, parameter :: NUC_MOD_OF = 3
    integer, parameter :: PTR_ZETA = 4
    
    integer, parameter :: ATM_SLOTS = 6

    integer, parameter :: ATOM_OF = 1
    integer, parameter :: ANG_OF = 2
    integer, parameter :: NPRIM_OF = 3
    integer, parameter :: NCTR_OF = 4
    integer, parameter :: KAPPA_OF = 5
    integer, parameter :: PTR_EXP = 6
    integer, parameter :: PTR_COEFF = 7
    integer, parameter :: BAS_SLOTS = 8

    integer, parameter :: PTR_RINV_ORIG = 5
    integer, parameter :: PTR_RINV_ZETA = 8
    integer, parameter :: AS_RINV_ORIG_ATOM = 18
    integer, parameter :: PTR_ENV_START = 21

    type(atom_basis_info),dimension(:),allocatable :: atom_bas
    type(mole_input_info) :: mole

type,public :: cint_info
!for set_mole_bas_info_sub
    integer :: natom
    integer :: nshl
    integer :: nbas 
    integer :: nprm 
    integer,allocatable,dimension(:) :: cntr_odr
    integer,allocatable,dimension(:) :: angl
    integer,allocatable,dimension(:) :: shl_belong_to_atom
    integer,allocatable,dimension(:) :: charge
    
    integer,allocatable,dimension(:) :: shl_index
    integer,allocatable,dimension(:) :: atom_bas_index
    real(kind=f8),allocatable,dimension(:) :: mass
    real(kind=f8),allocatable,dimension(:) :: mass_iso
    real(kind=f8),allocatable,dimension(:,:) :: coor
    real(kind=f8),allocatable,dimension(:) :: expnt
    real(kind=f8),allocatable,dimension(:) :: coeff
    integer,allocatable,dimension(:) :: fz_core
!for set cint_info_sub
    integer, dimension(:,:), allocatable :: atm
    integer, dimension(:,:), allocatable :: bas
    real(kind=f8), dimension(:), allocatable :: env
    real(kind=f8),allocatable,dimension(:) :: zeta
!for cal_cint_sub
    real(kind=f8),allocatable,dimension(:,:) :: S
    real(kind=f8),allocatable,dimension(:,:) :: T
    real(kind=f8),allocatable,dimension(:,:) :: V
    real(kind=f8),allocatable,dimension(:,:,:,:) :: eri
!for cal_gradint_sub
    real(kind=f8),allocatable,dimension(:,:,:) :: ipovlp
    real(kind=f8),allocatable,dimension(:,:,:) :: ipkin
    real(kind=f8),allocatable,dimension(:,:,:) :: ipnuc
    real(kind=f8),allocatable,dimension(:,:,:) :: iprinv
    real(kind=f8),allocatable,dimension(:,:,:,:,:) :: iperi

    contains
        procedure,public :: set_mole_bas_info => set_mole_bas_info_sub
        procedure,public :: set_cint_info => set_cint_info_sub
        procedure,public :: dyall_nuc_mod => dyall_nuc_mod_sub
        procedure,public :: normalize => normalize_sub
        procedure,public :: cal_cint => cal_cint_sub
        procedure,public :: cal_gradint => cal_gradint_sub
        procedure,public :: cal_iprinv => cal_iprinv_sub
        procedure,public :: rinv_at_nucleus => rinv_at_nucleus_sub

end type 

contains
    subroutine set_mole_bas_info_sub(this,atom_bas,mole)
    implicit none 
    class(cint_info) :: this
    type(atom_basis_info),allocatable,dimension(:),intent(in) :: atom_bas
    type(mole_input_info),intent(in) :: mole
    integer :: status,i,j,off,ioff,nbas_per_shl,ibas

    this%natom = mole%natom
    this%nshl = 0
    this%nbas = 0
    this%nprm = 0
    do i=1,this%natom
        this%nshl = this%nshl + atom_bas(i)%atom_nshl
        this%nbas = this%nbas + atom_bas(i)%atom_nbas
        this%nprm = this%nprm + atom_bas(i)%atom_nprm
    end do
    allocate(this%cntr_odr(this%nshl),stat=status)
    allocate(this%shl_index(this%nshl),stat=status)
    allocate(this%charge(this%natom),stat=status)
    allocate(this%fz_core(this%natom),stat=status)
    allocate(this%angl(this%nshl),stat=status)
    allocate(this%shl_belong_to_atom(this%nshl),stat=status)
    allocate(this%expnt(this%nprm),stat=status)
    allocate(this%coeff(this%nprm),stat=status)
    allocate(this%mass(this%natom),stat=status)
    allocate(this%mass_iso(this%natom),stat=status)
    allocate(this%zeta(this%natom),stat=status)
    allocate(this%coor(this%natom,3),stat=status)
    allocate(this%atom_bas_index(this%natom),stat=status)
    
    ibas = 1
    do i=1,this%natom 
        this%atom_bas_index(i) = ibas
        ibas = ibas+atom_bas(i)%atom_nbas
    end do
    off = 1
    ioff = 1
    do i = 1,this%natom
        this%shl_belong_to_atom(off:) = i
        this%mass(i) = atom_bas(i)%atom_mass
        this%mass_iso(i) = atom_bas(i)%atom_mass_iso
        this%fz_core(i) = atom_bas(i)%atom_fz_core
        this%cntr_odr(off:(off+atom_bas(i)%atom_nshl-1)) = atom_bas(i)%prm_num
        this%angl(off:(off+atom_bas(i)%atom_nshl-1)) = atom_bas(i)%bas_angl
        off = off + atom_bas(i)%atom_nshl
        this%charge(i) = atom_bas(i)%atom_charge
        this%expnt(ioff:(ioff+atom_bas(i)%atom_nprm-1)) = atom_bas(i)%expnt
        this%coeff(ioff:(ioff+atom_bas(i)%atom_nprm-1)) = atom_bas(i)%coeff
        this%coor(i,1) = mole%x_list(i)
        this%coor(i,2) = mole%y_list(i)
        this%coor(i,3) = mole%z_list(i)
        ioff = ioff + atom_bas(i)%atom_nprm
    end do 
    do i =1,this%nshl
        if (i == 1) then
            this%shl_index(i) = 1
        else 
		    select case (this%angl(i-1))
        	case (0); nbas_per_shl = 1
            case (1); nbas_per_shl = 3
            case (2); nbas_per_shl = 5   
            end select     
            this%shl_index(i) = this%shl_index(i-1) + nbas_per_shl
        end if 
    end do 
    end subroutine

    subroutine set_cint_info_sub(this,mole,nucl_mod)
    implicit none
    class(cint_info) :: this
    type(mole_input_info),intent(in) :: mole
    integer,intent(in) :: nucl_mod
    integer :: off, prim_off,status
    integer :: iatom, ishl, iprim, ioff
    real(kind=f8), external :: CINTgto_norm
    
    allocate(this%atm(ATM_SLOTS,this%natom),stat=status)
    allocate(this%bas(BAS_SLOTS,this%nshl),stat=status)
    allocate(this%env(PTR_ENV_START+4*this%natom+2*this%nprm),stat=status)
    this%env(1:21) = 0.0_f8
    off = PTR_ENV_START
    do iatom = 1, this%natom
        this%atm(CHARGE_OF,iatom) = this%charge(iatom)
        this%atm(PTR_COORD,iatom) = off
        this%atm(NUC_MOD_OF,iatom) = nucl_mod
        if (nucl_mod == 1) then
            this%zeta(iatom) = 0.0_f8
        else
            call this%dyall_nuc_mod(this%mass_iso(iatom),iatom)
        end if
        this%atm(PTR_ZETA,iatom) = off + 4
        this%atm(5,iatom) = 0
        this%atm(ATM_SLOTS,iatom) = 0 
        this%env(off+1) = mole%x_list(iatom)/0.529177249_f8
        this%env(off+2) = mole%y_list(iatom)/0.529177249_f8
        this%env(off+3) = mole%z_list(iatom)/0.529177249_f8
        this%env(off+4) = this%zeta(iatom)
        off = off + 4
    end do
    prim_off = 1
    do ishl = 1, this%nshl
        this%bas(ATOM_OF,ishl) = this%shl_belong_to_atom(ishl) - 1
        this%bas(ANG_OF,ishl) = this%angl(ishl)
        this%bas(NPRIM_OF,ishl) = this%cntr_odr(ishl)
        this%bas(NCTR_OF,ishl) = 1
        this%bas(KAPPA_OF,ishl) = 0
        this%bas(PTR_EXP,ishl) = off
        this%bas(BAS_SLOTS,ishl) = 0
        ioff = 0
        do iprim = prim_off, prim_off + this%cntr_odr(ishl)-1
            this%env(off+1+ioff) = this%expnt(iprim)
            ioff = ioff + 1
        end do
        off = off + this%cntr_odr(ishl)
        this%bas(PTR_COEFF,ishl) = off
        ioff = 0
        do iprim = prim_off, prim_off + this%cntr_odr(ishl)-1
            this%env(off+1+ioff) = this%coeff(iprim) * CINTgto_norm(this%angl(ishl), this%expnt(iprim))
            ioff = ioff + 1
        end do
        off = off + this%cntr_odr(ishl)
        prim_off = prim_off + this%cntr_odr(ishl)
    end do
    end subroutine 

    subroutine dyall_nuc_mod_sub(this,mass,iatom)
        implicit none
        class(cint_info) :: this
        real(kind=f8),intent(in) :: mass
        integer,intent(in) :: iatom
        real(kind=f8) :: r

        r = (0.836_f8*(mass**(1.0/3.0))+0.570_f8)/52917.7249_f8
        this%zeta(iatom) = 1.5_f8/(r**2.0)
        end subroutine

    subroutine normalize_sub(this)
        implicit none
        class(cint_info) :: this
        integer :: ishl, iprim
        integer :: di, dj
        real(kind=f8), dimension(:,:), allocatable :: buf1e
        integer, dimension(2) :: shls
        integer, external :: CINTcgto_spheric
  
        do ishl = 1, this%nshl
            shls(1) = ishl - 1
            shls(2) = ishl - 1
            di = CINTcgto_spheric(ishl-1, this%bas)
            dj = CINTcgto_spheric(ishl-1, this%bas)
            allocate(buf1e(di,dj))
            call cint1e_ovlp_sph(buf1e,shls,this%atm,this%natom,this%bas,this%nshl,this%env)
            do iprim = 1, this%cntr_odr(ishl)
                this%env(this%bas(PTR_COEFF,ishl)+iprim) = this%env(this%bas(PTR_COEFF,ishl)+iprim)/sqrt(buf1e(1,1))
            end do
            deallocate(buf1e)
        end do
    end subroutine 

    subroutine cal_cint_sub(this)
        implicit none 
        class(cint_info) :: this
        integer :: i, j, k, l, di, dj, dk, dl,status,x,y,z,w
        real(kind=f8), dimension(:,:), allocatable :: buf1e
        real(kind=f8), dimension(:,:,:,:), allocatable :: buf2e
        integer, dimension(4) :: shls
        integer(kind=i8) :: opt
        integer, external :: CINTcgto_spheric

        allocate(this%S(this%nbas,this%nbas),stat=status)
        allocate(this%T(this%nbas,this%nbas),stat=status)
        allocate(this%V(this%nbas,this%nbas),stat=status)
        allocate(this%eri(this%nbas,this%nbas,this%nbas,this%nbas),stat=status)
        call this%normalize()
        call cint2e_sph_optimizer(opt,this%atm,this%natom,this%bas,this%nshl,this%env)
        do i = 1, this%nshl
            shls(1) = i - 1
            di = CINTcgto_spheric(i-1, this%bas)
            do j = 1, this%nshl
                shls(2) = j - 1
                dj = CINTcgto_spheric(j-1, this%bas)
                allocate(buf1e(di,dj))
                x = this%shl_index(i)
                y = this%shl_index(j)
                call cint1e_ovlp_sph(buf1e,shls,this%atm,this%natom,this%bas,this%nshl,this%env)
                this%S(x:(x+di-1),y:(y+dj-1)) = buf1e(:,:)
                call cint1e_kin_sph(buf1e,shls,this%atm,this%natom,this%bas,this%nshl,this%env)
                this%T(x:(x+di-1),y:(y+dj-1)) = buf1e(:,:)
                call cint1e_nuc_sph(buf1e,shls,this%atm,this%natom,this%bas,this%nshl,this%env)
                this%V(x:(x+di-1),y:(y+dj-1)) = buf1e(:,:)
                deallocate(buf1e)
            end do
        end do
        do i = 1, this%nshl
            shls(1) = i - 1
            di = CINTcgto_spheric(i-1, this%bas)
            do j = 1, this%nshl
                shls(2) = j - 1
                dj = CINTcgto_spheric(j-1, this%bas)
                do k = 1, this%nshl
                    shls(3) = k - 1
                    dk = CINTcgto_spheric(k-1, this%bas)
                    do l = 1, this%nshl
                        shls(4) = l - 1
                        dl = CINTcgto_spheric(l-1, this%bas)
                        allocate(buf2e(di,dj,dk,dl))
                        x = this%shl_index(i)
                        y = this%shl_index(j)
                        z = this%shl_index(k)
                        w = this%shl_index(l)
                        call cint2e_sph(buf2e,shls,this%atm,this%natom,this%bas,this%nshl,this%env,opt)
                        this%eri(x:(x+di-1),y:(y+dj-1),z:z+dk-1,w:w+dl-1) = buf2e(:,:,:,:)
                        call CINTdel_optimizer(opt)
                        deallocate(buf2e)
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine cal_gradint_sub(this)
    implicit none
    class(cint_info) :: this
    integer :: i, j, k, l, di, dj, dk, dl,status,x,y,z,w
    real(kind=f8), dimension(:,:,:), allocatable :: buf1e
    real(kind=f8), dimension(:,:,:,:,:), allocatable :: buf2e
    integer, dimension(4) :: shls
    integer(kind=i8) :: opt
    integer, external :: CINTcgto_spheric

    allocate(this%ipovlp(this%nbas,this%nbas,3),stat=status)
    allocate(this%ipkin(this%nbas,this%nbas,3),stat=status)
    allocate(this%ipnuc(this%nbas,this%nbas,3),stat=status)
    allocate(this%iperi(this%nbas,this%nbas,this%nbas,this%nbas,3),stat=status)
    !call this%normalize()
        do i = 1, this%nshl
            shls(1) = i - 1
            di = CINTcgto_spheric(i-1, this%bas)
            do j = 1, this%nshl
                shls(2) = j - 1
                dj = CINTcgto_spheric(j-1, this%bas)
                allocate(buf1e(di,dj,3))
                x = this%shl_index(i)
                y = this%shl_index(j)
                call cint1e_ipovlp_sph(buf1e,shls,this%atm,this%natom,this%bas,this%nshl,this%env)
                this%ipovlp(x:(x+di-1),y:(y+dj-1),:) = buf1e(:,:,:)
                call cint1e_ipkin_sph(buf1e,shls,this%atm,this%natom,this%bas,this%nshl,this%env)
                this%ipkin(x:(x+di-1),y:(y+dj-1),:) = buf1e(:,:,:)
                call cint1e_ipnuc_sph(buf1e,shls,this%atm,this%natom,this%bas,this%nshl,this%env)
                this%ipnuc(x:(x+di-1),y:(y+dj-1),:) = buf1e(:,:,:)
                deallocate(buf1e)
            end do
        end do
        call cint2e_ip1_sph_optimizer(opt,this%atm,this%natom,this%bas,this%nshl,this%env)
        do i = 1, this%nshl
            shls(1) = i - 1
            di = CINTcgto_spheric(i-1, this%bas)
            do j = 1, this%nshl
                shls(2) = j - 1
                dj = CINTcgto_spheric(j-1, this%bas)
                do k = 1, this%nshl
                    shls(3) = k - 1
                    dk = CINTcgto_spheric(k-1, this%bas)
                    do l = 1, this%nshl
                        shls(4) = l - 1
                        dl = CINTcgto_spheric(l-1, this%bas)
                        allocate(buf2e(di,dj,dk,dl,3))
                        x = this%shl_index(i)
                        y = this%shl_index(j)
                        z = this%shl_index(k)
                        w = this%shl_index(l)
                        call cint2e_ip1_sph(buf2e,shls,this%atm,this%natom,this%bas,this%nshl,this%env,opt)
                        this%iperi(x:(x+di-1),y:(y+dj-1),z:z+dk-1,w:w+dl-1,:) = buf2e(:,:,:,:,:)
                        call CINTdel_optimizer(opt)
                        deallocate(buf2e)
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine cal_iprinv_sub(this)
    implicit none
    class(cint_info) :: this
    integer :: i,j,di,dj,status,x,y
    real(kind=f8), dimension(:,:,:), allocatable :: buf1e
    integer, dimension(4) :: shls
    integer, external :: CINTcgto_spheric

    allocate(this%iprinv(this%nbas,this%nbas,3),stat=status)
        do i = 1, this%nshl
            shls(1) = i - 1
            di = CINTcgto_spheric(i-1, this%bas)
            do j = 1, this%nshl
                shls(2) = j - 1
                dj = CINTcgto_spheric(j-1, this%bas)
                allocate(buf1e(di,dj,3))
                x = this%shl_index(i)
                y = this%shl_index(j)
                call cint1e_iprinv_sph(buf1e,shls,this%atm,this%natom,this%bas,this%nshl,this%env)
                this%iprinv(x:(x+di-1),y:(y+dj-1),:) = buf1e(:,:,:)
                deallocate(buf1e)
            end do
        end do
    end subroutine

    subroutine rinv_at_nucleus_sub(this,iatom)
    implicit none
    class(cint_info) :: this
    integer,intent(in) :: iatom

    this%env(AS_RINV_ORIG_ATOM) = iatom
    this%env(PTR_RINV_ZETA) = 0.0_f8
    this%env(PTR_RINV_ORIG) = this%coor(iatom,1)/0.529177249_f8
    this%env(PTR_RINV_ORIG+1) = this%coor(iatom,2)/0.529177249_f8
    this%env(PTR_RINV_ORIG+2) = this%coor(iatom,3)/0.529177249_f8
    end subroutine



end module