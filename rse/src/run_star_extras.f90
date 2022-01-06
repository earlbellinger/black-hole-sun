! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none
      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         real(dp) :: M_BH, L_BH, new_core_mass, core_avg_eps, &
            core_avg_rho, R_B, M_dot_Bondi, L_Bondi, L_Eddington, kap_face, rad_eff, con_eff
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% x_logical_ctrl(1)) then
             write(*, *) '** extras_startup'
             write(*, *) 'put the black hole in the center'
             write(*, *) 'calculate the Bondi radius and the luminosities'
             write(*, *) 'assume that there is a cavity of equal mass to the black hole'
             
             rad_eff = s% x_ctrl(1) ! radiative efficiency 
             con_eff = 0.1 !s% x_ctrl(2) ! convective efficiency 
             !inner_bound = s% x_ctrl(2) ! r_0 = 1/b * Bondi radius
             
             new_core_mass = s% job% new_core_mass
             M_BH = Msun * new_core_mass / 2 ! black hole seed mass 
             
             !kap_face = interp_val_to_pt(s% opacity, s% nz, s% nz, s% dq, 'kap_face')
             kap_face = s% opacity(s% nz)
             L_BH = pi4 * clight * s% cgrav(s% nz) * M_BH / kap_face ! Eddington luminosity 
             L_Eddington = L_BH
             !L_BH = 0.13 / kap_face * M_BH / (1d-5 * Msun) * Lsun ! Clayton 1975
             
             M_dot_Bondi = 16*pi*con_eff/(rad_eff / (1-rad_eff))/s% gamma1(s% nz) * pow(s% cgrav(s% nz) * M_BH, 2) / s% csound(s% nz) * s% rho(s% nz) / pow(clight, 2)
             L_Bondi = rad_eff * M_dot_Bondi * pow(clight, 2)
             
             s% xtra(8) = 0
             if (L_Bondi < L_BH) then
                write(*, *) 'Bondi accretion'
                L_BH = L_Bondi
                s% xtra(8) = 1
                else
                write(*, *) 'Eddington accretion'
             end if 
             
             write(*, *) 'initial M_center', s% M_center
             write(*, *) 'initial csound(s% nz)', s% csound(s% nz)
             write(*, *) 'initial kap_face', kap_face
             
             core_avg_eps = L_BH / (new_core_mass * Msun)
             
             ! Bondi radius, Warrick Ball's Ph.D. thesis Eqn. (3.2) 
             R_B = 2 * s% cgrav(s% nz) * M_BH / pow(s% csound(s% nz), 2)
             core_avg_rho = 1/(4/3*pi) * new_core_mass * Msun / pow(R_B, 3)
             
             s% xtra(1) = M_BH
             s% xtra(2) = L_BH
             s% xtra(3) = R_B
             s% xtra(6) = s% X(1) ! initial hydrogen abundance of the photosphere 
             s% xtra(7) = kap_face ! opacity of the centerpoint
             s% xtra(9) = L_Bondi
             s% xtra(10) = L_Eddington
             
             write(*, *) 'star mass', s% Mstar / Msun
             write(*, *) 'M_BH', M_BH / Msun
             write(*, *) 'L_BH', L_BH / Lsun
             write(*, *) 'R_B',  R_B  / Rsun
             
             write(*, *) 'new_core_mass', new_core_mass
             write(*, *) 'core_avg_eps',  core_avg_eps
             write(*, *) 'core_avg_rho',  core_avg_rho
             
             call star_relax_core( &
                   id, new_core_mass, s% job% dlg_core_mass_per_step, &
                   s% job% relax_core_years_for_dt, &
                   core_avg_rho, core_avg_eps, ierr)
         end if
         
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         
         real(dp) :: M_BH, L_BH, new_core_mass, core_avg_eps, &
            core_avg_rho, M_dot, M_dot_Bondi, L_Bondi, L_Eddington, rad_eff, con_eff, R_B, kap_face, M_BH_new, dm_dt
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
         
         ! star_jobs: new_core_mass,
         !            dlg_core_mass_per_step, 
         !            relax_core_years_for_dt
         !
         ! x_ctrl: (1) radiative efficiency 
         !         (2) convective efficiency (NOT USED)
         !         (3) Bondi radius factor (NOT USED)
         
         if (s% x_logical_ctrl(1)) then
             write(*, *) '** extras_start_step'
             write(*, *) 'accrete onto the black hole'
             write(*, *) 'use the Eddington luminosity to find the accretion rate M_dot'
             write(*, *) 'update the mass and then calculate the new Bondi radius'
             
             rad_eff = s% x_ctrl(1) ! radiative efficiency 
             con_eff = 0.1 !s% x_ctrl(2) ! convective efficiency 
             !inner_bound = s% x_ctrl(2) ! r_0 = 1/b * Bondi radius
             
             !! Calculate core_avg_eps, new_core_mass, and core_avg_rho
             
             M_BH = s% xtra(1) ! s% M_center / 2 ! g
             L_BH = s% xtra(2)
             
             ! get new black hole mass from the accretion rate 
             M_dot = L_BH / (rad_eff * pow(clight, 2)) ! g/s 
             dm_dt = (1 - rad_eff) * M_dot * s% dt ! g 
             M_BH_new = M_BH + dm_dt ! g 
             
             !new_core_mass = s% M_center / Msun + dm_dt / Msun
             new_core_mass = 2 * M_BH_new / Msun ! assume cavity is of equal mass to BH 
             
             ! recalculate luminosity based off the new BH mass 
             kap_face = s% opacity(s% nz)
             L_BH = pi4 * clight * s% cgrav(s% nz) * M_BH_new / kap_face ! erg/s
             L_Eddington = L_BH
             
             M_dot_Bondi = 16*pi*con_eff/(rad_eff / (1-rad_eff))/s% gamma1(s% nz) * pow(s% cgrav(s% nz) * M_BH_new, 2) / s% csound(s% nz) * s% rho(s% nz) / pow(clight, 2)
             L_Bondi = rad_eff * M_dot_Bondi * pow(clight, 2)
             
             s% xtra(8) = 0
             write(*, *) "L_BH", L_BH
             write(*, *) "L_B" , L_Bondi
             if (L_Bondi < L_BH) then
                write(*, *) 'Bondi accretion'
                L_BH = L_Bondi
                s% xtra(8) = 1
                else
                write(*, *) 'Eddington accretion'
             end if 
             
             core_avg_eps = L_BH / (new_core_mass * Msun)
             
             ! calculate Bondi radius, Warrick Ball's Ph.D. thesis Eqn. (3.2) 
             R_B = 2 * s% cgrav(s% nz) * M_BH_new / pow(s% csound(s% nz), 2) ! cm 
             core_avg_rho = 1/(4/3*pi) * new_core_mass * Msun / pow(R_B, 3) ! g/cm^3
             
             s% xtra(1) = M_BH_new
             s% xtra(2) = L_BH
             s% xtra(3) = R_B
             s% xtra(4) = M_dot
             s% xtra(5) = dm_dt
             !s% xtra(6) ! initial X
             s% xtra(7) = kap_face ! opacity of the centerpoint
             !s% xtra(8) ! Bondi or Eddington accretion 
             s% xtra(9) = L_Bondi
             s% xtra(10) = L_Eddington
             
             write(*, *) 'star mass', s% Mstar / Msun
             write(*, *) 'M_BH',     M_BH      / Msun
             write(*, *) 'M_BH_new', M_BH_new  / Msun
             write(*, *) 'L_BH',     L_BH      / Lsun
             write(*, *) 'R_B',      R_B       / Rsun
             write(*, *) 'M_dot',    M_dot     / Msun
             
             write(*, *) 'opacity', s% opacity(s% nz)
             write(*, *) 'M_center',      s% M_center
             write(*, *) 'new_core_mass', new_core_mass
             write(*, *) 'core_avg_eps',  core_avg_eps
             write(*, *) 'core_avg_rho',  core_avg_rho
             
             call star_relax_core( &
                   id, new_core_mass, s% job% dlg_core_mass_per_step, &
                   s% job% relax_core_years_for_dt, &
                   core_avg_rho, core_avg_eps, ierr)
             
         end if
         
      end function extras_start_step

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 11
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         integer :: i
         real(dp) :: X0, mX0, rX0
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         X0 = s% xtra(6) 
         do i = s% nz, 1, -1
            if (X0 - s% X(i) < 0.001) then
                mX0 = s% m(i) / Msun
                rX0 = s% r(i) / Rsun
                exit 
            end if 
         end do
         
         if (s% x_logical_ctrl(1)) then 
             names(1) = "M_BH"
             vals(1) = s% xtra(1) / Msun  ! M_BH / Msun
             names(2) = "L_BH"
             vals(2)  = s% xtra(2) / Lsun ! L_BH / Lsun
             names(3) = "R_B"
             vals(3)  = s% xtra(3) / Rsun ! R_B / Rsun
             names(4) = "M_dot"
             vals(4)  = s% xtra(4) / Msun ! M_dot / Msun
             names(5) = "dm/dt"
             vals(5)  = s% xtra(5) / Msun ! dm_dt
             names(6) = "mX0"
             vals(6)  = mX0
             names(7) = "rX0"
             vals(7)  = rX0
             names(8) = "kap_center" ! central opacity
             vals(8)  = s% xtra(7) 
             names(9) = "Bondi"
             vals(9)  = s% xtra(8) ! Bondi accretion 
             names(10) = "L_B"
             vals(10) = s% xtra(9) / Lsun ! Bondi luminosity 
             names(11) = "L_E"
             vals(11) = s% xtra(10) / Lsun ! Eddington luminosity 
         end if
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      real(dp) function interp_val_to_pt(v,k,sz,dq,str)
         use interp_1d_lib, only: interp_4_to_1
         integer, intent(in) :: k, sz
         real(dp), intent(in) :: v(:), dq(:)
         character (len=*), intent(in) :: str
         integer :: ierr
         include 'formats'
         if (k == 1) then
            interp_val_to_pt = v(k)
            return
         end if
         if (k > 2 .and. k < sz) then
            ierr = 0
            call interp_4_to_1( &
               0.5d0*(dq(k-2)+dq(k-1)), &
               0.5d0*(dq(k-1)+dq(k)), &
               0.5d0*(dq(k)+dq(k+1)), &
               0.5d0*dq(k-2)+dq(k-1), &
               v(k-2), v(k-1), v(k), v(k+1), &
               interp_val_to_pt, str, ierr)
            if (ierr == 0) return
            write(*,1) '0.5d0*(dq(k-2)+dq(k-1))', 0.5d0*(dq(k-2)+dq(k-1))
            write(*,1) '0.5d0*(dq(k-1)+dq(k))', 0.5d0*(dq(k-1)+dq(k))
            write(*,1) '0.5d0*(dq(k)+dq(k+1))', 0.5d0*(dq(k)+dq(k+1))
            write(*,2) 'dq(k-2)', k-2, dq(k-2)
            write(*,2) 'dq(k-1)', k-1, dq(k-1)
            write(*,2) 'dq(k)', k, dq(k)
            write(*,2) 'dq(k+1)', k+1, dq(k+1)

            call mesa_error(__FILE__,__LINE__,'interp_val_to_pt')
         endif
         interp_val_to_pt = (v(k)*dq(k-1) + v(k-1)*dq(k))/(dq(k-1) + dq(k))
      end function interp_val_to_pt


      end module run_star_extras
