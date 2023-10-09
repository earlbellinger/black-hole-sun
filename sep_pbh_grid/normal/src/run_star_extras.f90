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
      !use auto_diff
      
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
            core_avg_rho, R_B, M_dot_Bondi, L_Bondi, L_Eddington, rad_eff, con_eff, M_cav
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% x_logical_ctrl(1)) then
             rad_eff = s% x_ctrl(1) ! radiative efficiency 
             con_eff = 0.1 !s% x_ctrl(2) ! convective efficiency 
             
             M_BH = Msun * s% job% new_core_mass ! black hole seed mass 
             R_B = 2 * s% cgrav(s% nz) * M_BH / pow(s% csound(s% nz), 2)
             M_cav = 8*pi/3 * s% rho(s% nz) * pow(R_B, 3)
             new_core_mass = (M_BH + M_cav) / Msun
             
             L_Eddington = pi4 * clight * s% cgrav(s% nz) * M_BH / s% opacity(s% nz)
             M_dot_Bondi = 16*pi*con_eff/(rad_eff / (1-rad_eff))/s% gamma1(s% nz) * pow(s% cgrav(s% nz) * M_BH, 2) / s% csound(s% nz) * s% rho(s% nz) / pow(clight, 2)
             L_Bondi = rad_eff * M_dot_Bondi * pow(clight, 2)
             
             L_BH = L_Eddington
             s% xtra(8) = 0
             if (L_Bondi < L_Eddington) then
                L_BH = L_Bondi
                s% xtra(8) = 1
             end if 
             
             core_avg_eps = L_BH / (new_core_mass * Msun)
             core_avg_rho = 1/(4/3*pi) * (new_core_mass * Msun) / pow(R_B, 3) ! g/cm^3
             
             s% xtra(1) = M_BH
             s% xtra(2) = L_BH
             s% xtra(3) = R_B
             !s% xtra(4) = M_dot
             !s% xtra(5) = dm/dt
             s% xtra(6) = s% X(1) ! initial hydrogen abundance of the photosphere 
             s% xtra(7) = s% opacity(s% nz) ! opacity of the centerpoint
             !s% xtra(8) = Bondi or Eddington accretion
             s% xtra(9) = L_Bondi
             s% xtra(10) = L_Eddington
             s% xtra(11) = M_cav
             
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
            core_avg_rho, M_dot, M_dot_Bondi, L_Bondi, L_Eddington, rad_eff, con_eff, R_B, M_BH_new, dm, M_cav
         
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
             rad_eff = s% x_ctrl(1) ! radiative efficiency 
             con_eff = 0.1 !s% x_ctrl(2) ! convective efficiency 
             
             !! Calculate core_avg_eps, new_core_mass, and core_avg_rho
             
             M_BH = s% xtra(1) 
             L_BH = s% xtra(2) 
             
             ! get new black hole mass from the accretion rate 
             M_dot = L_BH / (rad_eff * pow(clight, 2)) ! g/s 
             dm = (1 - rad_eff) * M_dot * s% dt ! g 
             M_BH_new = M_BH + dm ! g 
             
             R_B = 2 * s% cgrav(s% nz) * M_BH_new / pow(s% csound(s% nz), 2) ! cm 
             M_cav = 8*pi/3 * s% rho(s% nz) * pow(R_B, 3)
             new_core_mass = (M_BH_new + M_cav) / Msun
             
             ! recalculate luminosity based off the new BH mass 
             L_Eddington = pi4 * clight * s% cgrav(s% nz) * M_BH_new / s% opacity(s% nz) ! erg/s
             M_dot_Bondi = 16*pi*con_eff/(rad_eff / (1-rad_eff))/s% gamma1(s% nz) * pow(s% cgrav(s% nz) * M_BH_new, 2) / s% csound(s% nz) * s% rho(s% nz) / pow(clight, 2)
             L_Bondi = rad_eff * M_dot_Bondi * pow(clight, 2)
             
             L_BH = L_Eddington
             s% xtra(8) = 0
             if (L_Bondi < L_BH) then
                L_BH = L_Bondi
                s% xtra(8) = 1
                else
             end if 
             
             s% max_timestep = M_BH / ((1-rad_eff) * M_dot)
             
             core_avg_eps = L_BH / (new_core_mass * Msun)
             core_avg_rho = 1/(4/3*pi) * (new_core_mass * Msun) / pow(R_B, 3) ! g/cm^3
             
             s% xtra(1) = M_BH_new
             s% xtra(2) = L_BH
             s% xtra(3) = R_B
             s% xtra(4) = M_dot
             s% xtra(5) = dm/s%dt
             !s% xtra(6) ! initial X
             s% xtra(7) = s% opacity(s% nz) 
             !s% xtra(8) ! Bondi or Eddington accretion 
             s% xtra(9) = L_Bondi
             s% xtra(10) = L_Eddington
             s% xtra(11) = M_cav
             
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
         !how_many_extra_history_columns = 0
         how_many_extra_history_columns = 13
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
         
         names(1:13) = 'empty'
         vals(1:13) = -1d99

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
             names(12) = "M_cav"
             vals(12) = s% xtra(11) / Msun ! M_cav / Msun
             names(13) = "cs_center"
             vals(13) = s% csound(s% nz)
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
         how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
                                              !(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, op_err, net_lwork
         logical :: okay
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
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

      end module run_star_extras
