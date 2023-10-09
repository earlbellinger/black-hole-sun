# Black Hole Sun

**[Bellinger, E. P.](https://earlbellinger.com), Caplan, M. E., Ryu, T., Bollimpalli, D., Ball, W. H., Kühnel, F., Farmer, R., de Mink, S. E., Christensen-Dalsgaard, J.** (2023). *[Solar evolution models with a central black hole]()*. ApJ accepted.

```
In my eyes
Indisposed
In disguises no one knows
Hides the face
Lies the snake
And the sun in my disgrace
Boiling heat
Summer stench
'Neath the black, the sky looks dead
Call my name
Through the cream
And I'll hear you scream again

Black hole sun
Won't you come
And wash away the rain?
Black hole sun
Won't you come
Won't you come
Won't you come?
```

The supplied notebooks generate all the figures in the paper. 

```
Stuttering
Cold and damp
Steal the warm wind, tired friend
Times are gone
For honest men
Sometimes, far too long for snakes
In my shoes
Walking sleep
In my youth, I pray to keep
Heaven send
Hell away
No one sings like you anymore

Black hole sun
Won't you come

And wash away the rain?
Black hole sun
Won't you come
Won't you come?
```

The crux of the code is a modification to the [MESA](https://mesa.sourceforge.net/) `run_star_extras.f90` file, the most important contents of which is the following FORTRAN subroutine to model the growth of the black hole: 
```fortran
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
          core_avg_eps = L_BH / (new_core_mass * Msun) ! average specific energy generation rate (erg/g s)
          core_avg_rho = 1/(4/3*pi) * (new_core_mass * Msun) / pow(R_B, 3) ! avg core density (g/cm^3)
          
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
```

```
Hang my head
Drown my fear
Till you all just disappear
Black hole sun
Won't you come
```