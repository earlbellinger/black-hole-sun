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
      
      subroutine do1_relax_R_center(s, new_Rcenter, ierr)
         ! adjust all lnR's to keep same density for each cell as 1st guess for next model
         type (star_info), pointer :: s
         real(dp), intent(in) :: new_Rcenter ! cm
         integer, intent(out) :: ierr
         real(dp) :: dm, rho, dr3, rp13
         integer :: k
         ierr = 0
         s% R_center = new_Rcenter
         ! adjust lnR's
         rp13 = s% R_center*s% R_center*s% R_center
         do k = s% nz, 1, -1
            dm = s% dm(k)
            rho = s% rho(k)
            dr3 = dm/(rho*four_thirds_pi) ! dm/rho is cell volume
            s% xh(s% i_lnR,k) = log(rp13 + dr3)*one_third
            rp13 = rp13 + dr3
         end do
      end subroutine do1_relax_R_center
      
      subroutine black_hole_accretion(id, s, startup, ierr)
          integer, intent(in) :: id
          logical, intent(in) :: startup
          type (star_info), pointer :: s
          integer, intent(out) :: ierr
          
          real(dp) :: G, c2, c_s, rho, gamma1, opacity, dt
          real(dp) :: nabla_ad, P_rad, P_gas
          real(dp) :: M_BH, M_BH_new, M_cav
          real(dp) :: M_dot_BH, M_dot, dm, R_B
          real(dp) :: L_Bondi, L_Edd, L_BH
          real(dp) :: rad_eff, con_eff, timestep_factor
          real(dp) :: core_avg_rho, core_avg_eps, new_core_mass
          
          rad_eff = s% x_ctrl(1) ! epsilon 
          con_eff = s% x_ctrl(2) ! eta 
          timestep_factor = s% x_ctrl(3)
          
          dt      = s% dt             ! time step              (s)
          G       = s% cgrav(s% nz)   ! gravitational constant (cm^3 / g s^2)
          c_s     = s% csound(s% nz)  ! speed of sound         (cm / s)
          rho     = s% rho(s% nz)     ! density                (g / cm^3)
          opacity = s% opacity(s% nz) ! opacity                (cm^2 / g)
          P_rad   = s% prad(s% nz)    ! radiation pressure     (dyn / cm^2)
          P_gas   = s% pgas(s% nz)    ! gas pressure           (dyn / cm^2)
          gamma1  = s% gamma1(s% nz)  ! adiabatic index
          nabla_ad = 1 - 1 / gamma1
          c2 = pow(clight, 2)
          
          M_BH = s% xtra(1) ! black hole mass (g)
          
          M_dot_BH = 16*pi / (rad_eff / (1 - rad_eff)) * con_eff/gamma1/c_s*rho * pow(G*M_BH, 2) / c2 ! g/s
          
          L_Edd = 4*pi * clight * G * M_BH / opacity  ! erg/s
          L_Bondi = (rad_eff / (1 - rad_eff)) * M_dot_BH * c2 ! erg/s
          L_BH = L_Bondi
          if (s% x_logical_ctrl(2)) L_BH = min(L_Bondi, L_Edd)
          
          M_dot = L_BH / (rad_eff * c2) ! g/s
          
          dm = (1 - rad_eff) * M_dot * dt  ! g 
          M_BH_new = M_BH + dm             ! g
          
          R_B = 2 * G * M_BH_new / pow(c_s, 2)         ! Bondi radius (cm)
          M_cav = 8 * pi / 3 * rho * pow(R_B, 3)       ! mass of cavity (g)
          new_core_mass = (M_BH_new + M_cav) / Msun    ! new core mass (Msun)
          core_avg_eps = L_BH / (new_core_mass * Msun) ! average energy generation rate (erg / g s)
          core_avg_rho = 1 / (4 / 3 * pi) * (new_core_mass * Msun) / pow(R_B, 3) ! average core density (g / cm^3)
          
          s% max_timestep = timestep_factor * M_BH / ((1 - rad_eff) * M_dot) ! maximum timestep (s)
          
          s% xtra(1)  = M_BH_new
          s% xtra(2)  = L_BH
          s% xtra(3)  = R_B
          s% xtra(4)  = M_dot
          s% xtra(5)  = safe_log10(dm) - safe_log10(dt)
          s% xtra(6)  = rad_eff
          s% xtra(7)  = opacity
          s% xtra(8)  = L_Bondi
          s% xtra(9)  = L_Edd
          s% xtra(10) = M_cav
          s% xtra(11) = P_rad
          s% xtra(12) = P_gas
          s% xtra(13) = nabla_ad
          
          if (startup) then
              call star_relax_core( &
                  id, new_core_mass, s% job% dlg_core_mass_per_step, &
                  s% job% relax_core_years_for_dt, &
                  core_avg_rho, core_avg_eps, ierr)
          else
              s% M_center = new_core_mass
              s% mstar = s% mstar - rad_eff * M_dot * dt ! remove mass converted into photons
              s% xmstar = s% mstar - s% M_center
              s% L_center = L_BH
              call do1_relax_R_center(s, R_B, ierr)
          end if 
          
          print*, '--- Black Hole Properties ---'
          print*, 'M/M_sun: ',    s% mstar / Msun
          print*, 'M_BH/M_sun: ', M_BH_new / Msun
          print*, 'L_BH/L_sun: ', L_BH / Lsun
          print*, 'radiative efficiency: ', rad_eff
          print*, '-----------------------------'
          
      end subroutine black_hole_accretion
      
      
      subroutine extras_startup(id, restart, ierr)
          integer, intent(in) :: id
          logical, intent(in) :: restart
          integer, intent(out) :: ierr
          type (star_info), pointer :: s

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          if (s% x_logical_ctrl(1)) then
              s% xtra(1) = s% job% new_core_mass * Msun
              call black_hole_accretion(id, s, .true., ierr) 
          end if
      end subroutine extras_startup
      
      
      integer function extras_start_step(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_start_step = 0

          if (s% x_logical_ctrl(1)) then
              call black_hole_accretion(id, s, .false., ierr) 
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
         
         ! terminate if the BH mass exceeds x_ctrl(4)
         if (s% x_ctrl(4) > 0 .and. s% xtra(1) / Msun > s% x_ctrl(4)) then
            extras_check_model = terminate
            termination_code_str(t_xtra1) = 'black hole'
         end if
         
         ! terminate if the BH mass exceeds the stellar mass 
         if (s% xtra(1) + s% xtra(11) >= s% mstar) then
            extras_check_model = terminate
            write(*, *) 'M_BH >= Mstar'
            termination_code_str(t_xtra1) = 'black hole'
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
         how_many_extra_history_columns = 14
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
         
         names(1:14) = 'empty'
         vals(1:14) = -1d99
         
         if (s% x_logical_ctrl(1)) then 
             names(1) = "M_BH"
             vals(1) = s% xtra(1) / Msun   ! M_BH / Msun
             names(2) = "L_BH"
             vals(2)  = s% xtra(2) / Lsun  ! L_BH / Lsun
             names(3) = "R_B"
             vals(3)  = s% xtra(3) / Rsun  ! R_B  / Rsun
             names(4) = "M_dot"
             vals(4)  = s% xtra(4) / Msun  ! M_dot / Msun
             names(5) = "log10(dm/dt)"
             vals(5)  = s% xtra(5)         ! g/s
             names(6) = "rad_eff"
             vals(6)  = s% xtra(6)         ! epsilon  
             names(7) = "kap_center"
             vals(7)  = s% xtra(7)         ! cm^2/g
             names(8) = "L_B"
             vals(8) = s% xtra(9) / Lsun   ! Bondi luminosity 
             names(9) = "L_E"
             vals(9) = s% xtra(10) / Lsun ! Eddington luminosity 
             names(10) = "M_cav"
             vals(10) = s% xtra(10) / Msun ! M_cav / Msun
             names(11) = "prad_center"
             vals(11)  = s% xtra(11)
             names(12) = "pgas_center"
             vals(12)  = s% xtra(12)
             names(13) = "nabla_ad_center"
             vals(13)  = s% xtra(13)
             names(14) = "cs_center"
             vals(14) = s% csound(s% nz)   ! cm/s
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
