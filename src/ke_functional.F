MODULE ke_functional
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_retain,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type,&
                                              pw_release
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE qs_environment_types,            ONLY: get_qs_env,&               
                                              qs_environment_type  
   USE dbcsr_api,                       ONLY: dbcsr_add,&                                            
                                              dbcsr_dot,&                                            
                                              dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE cp_control_types,                ONLY: dft_control_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&                                    
                                              pw_pool_give_back_pw,&                                 
                                              pw_pool_type

   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type

   USE qs_ks_types,                     ONLY: qs_ks_env_type
   use dbcsr_api,                       only: dbcsr_clear

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ke_functional'

   PUBLIC :: add_ke_func_contribution

CONTAINS

    subroutine add_ke_func_contribution(qs_env, ks_matrix, ekin)
    
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(KIND=dp), intent(inout)                          :: ekin                         

      TYPE(qs_rho_type), POINTER                         :: rho 
      TYPE(section_vals_type), POINTER                   :: kin_func_section, input
      integer                                            :: kinetic_op, ispin, nspins
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vkin_rho, &
                                                            vkin_tau
      REAL(KIND=dp)                                      :: alpha
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: density_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, matrix_s, matrix_t
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      LOGICAL                             :: calc_force
      
      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      rho=rho, &
                      input=input,&
                      dft_control=dft_control, &
                      kinetic_kp=matrix_t, &
                      pw_env=pw_env &
                      )

       ekin=0
       
       CALL section_vals_val_get(input, "DFT%KINETIC_OPERATOR", i_val=kinetic_op)
         
       print *, "Nimag ", dft_control % nimages

   
       if (kinetic_op /= 0 .and. dft_control % nimages == 1) then
           print *,"KINETIC FUNC" 
          
          call dbcsr_clear(matrix_t(1,1)%matrix)

          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

          kin_func_section => section_vals_get_subs_vals(input, "DFT%KG_METHOD%XC")


          CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=kin_func_section, &
                       vxc_rho=vkin_rho, vxc_tau=vkin_tau, exc=ekin)
          
          print *, "EKIN ", ekin
          CALL qs_rho_get(rho, rho_ao=density_matrix)

          alpha = 1.0
          if (dft_control%nspins==1) alpha = 2.0
          
          CPASSERT(dft_control%nspins==1)
          calc_force = .false.   
          DO ispin = 1, dft_control%nspins
      
            

             ! in KG my_alpha is added , why??
             vkin_rho(ispin)%pw%cr3d = alpha * vkin_rho(ispin)%pw%cr3d * vkin_rho(ispin)%pw%pw_grid%dvol
             
             ! Integrate kin-potential with external density for outer response forces
             CALL integrate_v_rspace(v_rspace=vkin_rho(ispin), &
                            pmat=density_matrix(ispin), hmat=ks_matrix(ispin), &
                            qs_env=qs_env,calculate_forces=calc_force)

             !CALL pw_release(vkin_rho(ispin)%pw)
             CALL pw_pool_give_back_pw(auxbas_pw_pool, vkin_rho(ispin)%pw)
          END DO
       endif 
 
    end subroutine
end module
