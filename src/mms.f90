! Manufactured solutions and source terms for operator and Navier-Stokes verification.
# include <petsc/finclude/petscsys.h>

module mms
   use types
   implicit none

   private

   public f1d
   public fone
   public fcos2
   public fsin2
   public fcos4
   public fcos4_2
   public rcos4_2
   public fsin4
   public fcos6
   public fsin6
   public fsin8
   public fcos8
   public fsin22
   public fcos22
   public fcos23
   public fcos42
   public fsin42
   public fcos43
   public fcos2cos4
   public fcos4cos4
   public fsin2sin4
   public fcos2cos2
   public fcos2sin2
   public fcos2sin4
   public fcos4sin2
   public fcos4sin4
   public fsin2cos4
   public fcos2sin8
   public fcos22cos4
   public fcos22sin4
   public fsin22cos4
   public fcos2sin22
   public fcos22sin2
   public fcos2cos42
   public fcos2sin42
   public fsin2cos42
   public fsin2cos23
   public fsin4cos42
   public fcos4sin42
   public fcos2sin2cos4
   public fcos2cos4sin4
   public fsin2sin4cos4
   public gsin4sin4sin4
   public gsin4cos4sin4
   public gsin4cos4sin4_2
   public rsin4cos4sin4_2
   public gcos4sin4cos4
   public gcos4sin4cos4_2
   public gcos2sin2cos2
   public gcos2sin2cos2_2
   public gcos4cos4cos4
   public gcos2cos4cos2
   public gcos2cos4cos2_2
   public rcos2cos4cos2_2
   public gcos4cos2cos4_2
   public gcos2cos4sin2sin4cos2cos4
   public gcos2cos4cos4cos4cos2cos4
   public gcos2cos4cos2cos4sin4cos4
   public grho_u
   public gcos2cos4cos2sin4cos4sin4
   public grho_v
   public gcos2cos4cos2cos2sin2cos2
   public grho_w
   public gcos2cos4sin4sin2sin4cos4cos2cos4sin4
   public dx_cos2cos4cos4cos4cos2cos4
   public dy_cos2cos4cos4cos4cos2cos4
   public dz_cos2cos4cos4cos4cos2cos4
   public div_cos2cos4cos4cos4cos2cos4
   public div_cos4cos4cos4
   public div_rho_vel
   public expt
   public INTexpt
   public u_convection
   public v_convection
   public w_convection
   public u_viscous
   public v_viscous
   public w_viscous
   public u_momentum
   public v_momentum
   public w_momentum
   public g_energy
   public fu
   public ftemp
contains

   function fu(v, t) result (integrand)
      PetscReal :: v, t, integrand

      PetscReal :: Umax

      Umax = 1.5
      integrand = Umax * (1.0 - (v - 1.0)**2) * (1.0 + t)
   end function fu

   function ftemp(v, t) result (integrand)
      PetscReal :: v, t, integrand

      PetscReal :: Umax

      Umax = 1.5
      integrand = Umax * (1.0 - (v - 1.0)**4)
   end function ftemp

   function f1d(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(4.0 * v) * EXP(-t)
   end function f1d

   function fone(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = 1.0
   end function fone

   function fcos2(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * EXP(-t)
   end function fcos2

   function fcos4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(4.0 * v) * EXP(-t)
   end function fcos4

   function fcos4_2(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = (COS(4.0 * v) + 2.0) * EXP(-t)
   end function fcos4_2

   function rcos4_2(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = 1.0 / ((COS(4.0 * v) + 2.0) * EXP(-t))
   end function rcos4_2

   function fsin2(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(2.0 * v) * EXP(-t)
   end function fsin2

   function fsin4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(4.0 * v) * EXP(-t)
   end function fsin4

   function fcos6(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(6.0 * v) * EXP(-t)
   end function fcos6

   function fsin6(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(6.0 * v) * EXP(-t)
   end function fsin6

   function fcos43(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(4.0 * v)**3 * EXP(-t)
   end function fcos43

   function fsin8(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(8.0 * v) * EXP(-t)
   end function fsin8

   function fcos8(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(8.0 * v) * EXP(-t)
   end function fcos8

   function fcos42(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(4.0 * v)**2 * EXP(-t)
   end function fcos42

   function fsin42(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(4.0 * v)**2 * EXP(-t)
   end function fsin42

   function fsin2cos4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(2.0 * v) * COS(4.0 * v) * EXP(-t)
   end function fsin2cos4

   function fcos2sin8(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * SIN(8.0 * v) * EXP(-t)
   end function fcos2sin8

   function fcos22cos4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v)**2 * COS(4.0 * v) * EXP(-t)
   end function fcos22cos4

   function fcos22sin4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v)**2 * SIN(4.0 * v) * EXP(-t)
   end function fcos22sin4

   function fcos2cos42(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * COS(4.0 * v)**2 * EXP(-t)
   end function fcos2cos42

   function fcos2sin42(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * SIN(4.0 * v)**2 * EXP(-t)
   end function fcos2sin42

   function fsin2cos42(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(2.0 * v) * COS(4.0 * v)**2 * EXP(-t)
   end function fsin2cos42

   function fsin2cos23(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(2.0 * v) * COS(2.0 * v)**3 * EXP(-t)
   end function fsin2cos23

   function fsin4cos42(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(4.0 * v) * COS(4.0 * v)**2 * EXP(-t)
   end function fsin4cos42

   function fcos4sin42(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(4.0 * v) * SIN(4.0 * v)**2 * EXP(-t)
   end function fcos4sin42

   function fcos2cos4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * COS(4.0 * v) * EXP(-t)
   end function fcos2cos4

   function fcos4cos4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(4.0 * v) * COS(4.0 * v) * EXP(-t)
   end function fcos4cos4

   function fsin2sin4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(2.0 * v) * SIN(4.0 * v) * EXP(-t)
   end function fsin2sin4

   function fcos2cos2(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * COS(2.0 * v) * EXP(-t)
   end function fcos2cos2

   function fcos2sin2(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * SIN(2.0 * v) * EXP(-t)
   end function fcos2sin2

   function fcos2sin4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * SIN(4.0 * v) * EXP(-t)
   end function fcos2sin4

   function fcos4sin2(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(4.0 * v) * SIN(2.0 * v) * EXP(-t)
   end function fcos4sin2

   function fcos4sin4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(4.0 * v) * SIN(4.0 * v) * EXP(-t)
   end function fcos4sin4

   function fsin22(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(2.0 * v)**2 * EXP(-t)
   end function fsin22

   function fcos22(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v)**2 * EXP(-t)
   end function fcos22

   function fcos23(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v)**3 * EXP(-t)
   end function fcos23

   function fcos2sin2cos4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * SIN(2.0 * v) * COS(4.0 * v) * EXP(-t)
   end function fcos2sin2cos4

   function fsin22cos4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(2.0 * v)**2 * COS(4.0 * v) * EXP(-t)
   end function fsin22cos4

   function fcos2sin22(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * SIN(2.0 * v)**2 * EXP(-t)
   end function fcos2sin22

   function fcos22sin2(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v)**2 * SIN(2.0 * v) * EXP(-t)
   end function fcos22sin2

   function fcos2cos4sin4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = COS(2.0 * v) * COS(4.0 * v) * SIN(4.0 * v) * EXP(-t)
   end function fcos2cos4sin4

   function fsin2sin4cos4(v, t) result (integrand)
      PetscReal :: v, t, integrand

      integrand = SIN(2.0 * v) * SIN(4.0 * v) * COS(4.0 * v) * EXP(-t)
   end function fsin2sin4cos4

   function gsin4sin4sin4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = SIN(4.0 * x) * SIN(4.0 * y) * SIN(4.0 * z) * EXP(-t)
   end function gsin4sin4sin4

   function gsin4cos4sin4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = SIN(4.0 * x) * COS(4.0 * y) * SIN(4.0 * z) * EXP(-t)
   end function gsin4cos4sin4

   function gsin4cos4sin4_2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = (SIN(4.0 * x) * COS(4.0 * y) * SIN(4.0 * z) + 2.0) * EXP(-t)
   end function gsin4cos4sin4_2

   function rsin4cos4sin4_2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = 1.0 / gsin4cos4sin4_2(x, y, z, t)
   end function rsin4cos4sin4_2

   function gcos4sin4cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(4.0 * x) * SIN(4.0 * y) * COS(4.0 * z) * EXP(-t)
   end function gcos4sin4cos4

   function gcos4sin4cos4_2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = (COS(4.0 * x) * SIN(4.0 * y) * COS(4.0 * z) + 2.0) * EXP(-t)
   end function gcos4sin4cos4_2

   function gcos2sin2cos2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(2.0 * x) * SIN(2.0 * y) * COS(2.0 * z) * EXP(-t)
   end function gcos2sin2cos2

   function gcos2sin2cos2_2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = (COS(2.0 * x) * SIN(2.0 * y) * COS(2.0 * z) + 2.0) * EXP(-t)
   end function gcos2sin2cos2_2

   function gcos4cos4cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(4.0 * x) * COS(4.0 * y) * COS(4.0 * z) * EXP(-t)
   end function gcos4cos4cos4

   function gcos2cos4cos2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z) * EXP(-t)
   end function gcos2cos4cos2

   function gcos2cos4cos2_2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = (COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z) + 2.0) * EXP(-t)
   end function gcos2cos4cos2_2

   function rcos2cos4cos2_2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = 1.0 / gcos2cos4cos2_2(x, y, z, t)
   end function rcos2cos4cos2_2

   function gcos4cos2cos4_2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = (COS(4.0 * x) * COS(2.0 * y) * COS(4.0 * z) + 2.0) * EXP(-t)
   end function gcos4cos2cos4_2

   function gcos2cos4sin2sin4cos2cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(2.0 * x) * SIN(2.0 * y) * COS(2.0 * z) * &
                  COS(4.0 * x) * SIN(4.0 * y) * COS(4.0 * z) * EXP(-t)
   end function gcos2cos4sin2sin4cos2cos4

   function gcos2cos4cos4cos4cos2cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z) * &
                  COS(4.0 * x) * COS(4.0 * y) * COS(4.0 * z) * EXP(-t)
   end function gcos2cos4cos4cos4cos2cos4

   function gcos2cos4cos2cos4sin4cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z) * &
                  COS(4.0 * x) * SIN(4.0 * y) * COS(4.0 * z) * EXP(-t)
   end function gcos2cos4cos2cos4sin4cos4

   function grho_u(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = gcos2cos4cos2_2(x, y, z, t) * gcos4sin4cos4(x, y, z, t)
   end function grho_u

   function gcos2cos4cos2sin4cos4sin4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z) * &
                  SIN(4.0 * x) * COS(4.0 * y) * SIN(4.0 * z) * EXP(-t)
   end function gcos2cos4cos2sin4cos4sin4

   function grho_v(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = gcos2cos4cos2_2(x, y, z, t) * gsin4cos4sin4(x, y, z, t)
   end function grho_v

   function gcos2cos4cos2cos2sin2cos2(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z) * &
                  COS(2.0 * x) * SIN(2.0 * y) * COS(2.0 * z) * EXP(-t)
   end function gcos2cos4cos2cos2sin2cos2

   function grho_w(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = gcos2cos4cos2_2(x, y, z, t) * gcos2sin2cos2(x, y, z, t)
   end function grho_w

   function gcos2cos4sin4sin2sin4cos4cos2cos4sin4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = COS(2.0 * x) * SIN(2.0 * y) * COS(2.0 * z) * &
                  COS(4.0 * x) * SIN(4.0 * y) * COS(4.0 * z) * &
                  SIN(4.0 * x) * COS(4.0 * y) * SIN(4.0 * z) * EXP(-t)
   end function gcos2cos4sin4sin2sin4cos4cos2cos4sin4

   function dx_cos2cos4cos4cos4cos2cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = - 0.5 * COS(4.0 * y)**2 * &
         (COS(2.0 * z) + COS(6.0 * z)) * &
         (SIN(2.0 * x) + 3.0 * SIN(6.0 * x)) * EXP(-t)
   end function dx_cos2cos4cos4cos4cos2cos4

   function dy_cos2cos4cos4cos4cos2cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = - (COS(2.0 * x) + COS(6.0 * x)) * &
                    (COS(2.0 * z) + COS(6.0 * z)) * &
                     SIN(8.0 * y) * EXP(-t)
   end function dy_cos2cos4cos4cos4cos2cos4

   function dz_cos2cos4cos4cos4cos2cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = - 0.5 * COS(4.0 * y)**2 * &
         (COS(2.0 * x) + COS(6.0 * x)) * &
         (SIN(2.0 * z) + 3.0 * SIN(6.0 * z)) * EXP(-t)
   end function dz_cos2cos4cos4cos4cos2cos4

   function div_cos2cos4cos4cos4cos2cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = dx_cos2cos4cos4cos4cos2cos4(x, y, z, t) + &
                  dy_cos2cos4cos4cos4cos2cos4(x, y, z, t) + &
                  dz_cos2cos4cos4cos4cos2cos4(x, y, z, t)
   end function div_cos2cos4cos4cos4cos2cos4

   function div_cos4cos4cos4(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = -4.0 * (COS(4.0 * x) * SIN(4.0 * y) * COS(4.0 * z) + &
         COS(4.0 * y) * SIN(4.0 * (x + z))) * EXP(-t)
   end function div_cos4cos4cos4

   function expt(x, t) result (integrand)
      PetscReal :: x, t, integrand

      integrand = EXP(-t**2)
   end function expt

   function INTexpt(x, t) result (integrand)
      PetscReal :: x, t, integrand

      integrand = 0.5 * SQRT(PETSC_PI) * ERF(t) + 1.0
   end function INTexpt

   function div_rho_vel(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = &
         -4.0 * EXP(-2.0 * t) * (2.0 + COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z)) * &
         COS(4.0 * z) * SIN(4.0 * x) * SIN(4.0 * y) - 2.0 * EXP(-2.0 * t) * &
         COS(4.0 * y) * COS(2.0 * z) * SIN(2.0 * x) * &
         (2.0 + COS(4.0 * x) * COS(4.0 * z) * SIN(4.0 * y)) &
         -4.0 * EXP(-2.0 * t) * (2.0 + COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z)) * &
         SIN(4.0 * x) * SIN(4.0 * y) * SIN(4.0 * z) - 4.0 * EXP(-2.0 * t) * &
         COS(2.0 * x) * SIN(4.0 * y) * COS(2.0 * z) * &
         (2.0 + SIN(4.0 * x) * COS(4.0 * y) * SIN(4.0 * z)) &
         -2.0 * EXP(-2.0 * t) * (2.0 + COS(2.0 * x) * COS(4.0 * y) * COS(2.0 * z)) * &
         COS(2.0 * x) * SIN(2.0 * y) * SIN(2.0 * z) - 2.0 * EXP(-2.0 * t) * &
         COS(2.0 * x) * COS(4.0 * y) * SIN(2.0 * z) * &
         (2.0 + COS(2.0 * x) * SIN(2.0 * y) * COS(2.0 * z))
   end function div_rho_vel

   function u_convection(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = Cos(2*x)*Cos(4*x)*Cos(4*y)*Cos(2*z)*&
      (4.*Cos(4*z)**2*Sin(4*x)*Sin(4*y)**2 - &
      4.*Cos(4*y)**2*Cos(4*z)*Sin(4*x)*Sin(4*z) + &
      4.*Cos(2*x)*Cos(2*z)*Sin(2*y)*Sin(4*y)*Sin(4*z))
   end function u_convection

   function v_convection(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = Cos(2*x)*Cos(4*y)**2*Cos(2*z)*&
         (-4.*Cos(2*x)*Cos(2*z)*Cos(4*z)*Sin(4*x)*Sin(2*y) + &
         Sin(4*y)*Sin(4*z)*(-4.*Cos(4*x)**2*Cos(4*z) + &
         4.*Sin(4*x)**2*Sin(4*z)))
   end function v_convection

   function w_convection(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = Cos(2*x)*Cos(4*y)*Cos(2*z)**2*&
         (2.*Cos(4*x)*Cos(4*z)*Sin(2*x)*Sin(2*y)*Sin(4*y) + &
         Cos(2*x)*(2.*Cos(2*x)*Sin(2*y)**2*Sin(2*z) -&
         2.*Cos(2*y)*Cos(4*y)*Sin(4*x)*Sin(4*z)))
   end function w_convection

   function u_viscous(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = (-160*Cos(4*x)*Cos(4*z)*Sin(4*y))/3. +  &
      (4*Sin(2*x)*Sin(2*y)*Sin(2*z) - 16*Cos(4*x)*Sin(4*y)*Sin(4*z))/3.
   end function u_viscous

   function v_viscous(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = (-16*Cos(4*y)*Cos(4*z)*Sin(4*x) - &
      4*Cos(2*x)*Cos(2*y)*Sin(2*z))/3. - (160*Cos(4*y)*Sin(4*x)*Sin(4*z))/3.
   end function v_viscous

   function w_viscous(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = (-40*Cos(2*x)*Cos(2*z)*Sin(2*y))/3. + &
      (-16*Cos(4*z)*Sin(4*x)*Sin(4*y) + 16*Sin(4*x)*Sin(4*y)*Sin(4*z))/3.
   end function w_viscous

   function u_momentum(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = &
   (120*Cos(2*y)*Cos(4*z)*Sin(4*x) - &
     84*Cos(2*x)*Cos(4*x)*Cos(4*y)**3*Cos(2*z)*Cos(4*z)*Sin(4*x)*Sin(4*z) + &
     28*(3*Cos(4*z)**2*Sin(8*x)*Sin(4*y)**2 + Sin(2*x)*Sin(2*y)*Sin(2*z) - &
     2*Cos(4*x)*Sin(4*y)*(20*Cos(4*z) +&
     (2 - 3*Cos(2*x)*Cos(2*z)*Sin(2*y))*Sin(4*z))) + &
     6*Cos(4*y)*Cos(2*z)*(5*(2 + Cos(4*x)*Cos(2*y)*Cos(4*z))*Sin(2*x) + &
     Cos(2*x)*(10*Cos(2*y)*Cos(4*z)*Sin(4*x) + &
     7*Sin(4*y)*(Cos(4*z)**2*Sin(8*x)*Sin(4*y) + &
     2*Cos(2*x)*Cos(4*x)*Cos(2*z)*Sin(2*y)*Sin(4*z)))) - &
     42*Cos(4*y)**2*Sin(8*x)*Sin(8*z))/(21.*(2 + Cos(2*x)*Cos(4*y)*Cos(2*z)))
   end function u_momentum

   function v_momentum(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = &
    (-2*(-60*Cos(2*x)*Cos(2*z)*Sin(4*y) - &
     30*Cos(4*x)*Cos(4*z)*(Sin(2*y) + Cos(2*x)*Cos(2*y)*Cos(2*z)*Sin(4*y)) + &
     14*Cos(2*x)*Cos(2*y)*Sin(2*z) - 42*Sin(4*x)**2*Sin(8*y)*Sin(4*z)**2 + &
     Cos(4*y)*(Cos(4*z)*(-15*Cos(2*x)*Cos(4*x)*Cos(2*z)*Sin(2*y) + &
           28*Sin(4*x)*(2 + 3*Cos(2*x)*Cos(2*z)*Sin(2*y))) + &
     560*Sin(4*x)*Sin(4*z)) + &
     21*Cos(4*x)**2*Sin(8*y)*Sin(8*z) + &
     21*Cos(2*x)*Cos(4*y)**2*Cos(2*z)*&
     (2*Cos(2*x)**2*(Cos(2*z) + Cos(6*z))*Sin(2*x)*Sin(2*y) + &
     Sin(4*y)*(-2*Sin(4*x)**2*Sin(4*z)**2 + Cos(4*x)**2*Sin(8*z)))))/ &
     (21.*(2 + Cos(2*x)*Cos(4*y)*Cos(2*z)))
   end function v_momentum

   function w_momentum(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = &
     (28*Cos(4*z)*(-4*Sin(4*x) + 3*Cos(4*x)*Cos(2*z)*Sin(2*x)*Sin(2*y))*Sin(4*y) + &
      42*Cos(2*x)**3*Cos(4*y)*Cos(2*z)**2* &
      (-4*Cos(2*y)*Cos(4*y)*Cos(2*z)*Sin(2*x) + Sin(2*y)**2)*Sin(2*z) + &
      42*Cos(2*x)**2*Sin(2*y)**2*Sin(4*z) + &
      2*Cos(2*x)*(15*Cos(4*y)*(2 + Cos(4*x)*Cos(2*y)*Cos(4*z))*Sin(2*z) - &
      2*Cos(2*z)*(70*Sin(2*y) + 3*Cos(2*y)*Cos(4*y)*(-5*Cos(4*x) + &
      7*Sin(4*x))*Sin(4*z))) + &
      (Sin(4*z)*(3840*Cos(4*x)*Cos(2*y) + 3584*Sin(4*x)*Sin(4*y) + &
      21/(Sin(2*z)**2) *Sin(8*x)*Sin(2*y)*Sin(8*y)*Sin(8*z)))/32.)/ &
      (21.*(2 + Cos(2*x)*Cos(4*y)*Cos(2*z)))
   end function w_momentum

   function g_energy(x, y, z, t) result (integrand)
      PetscReal :: x, y, z, t, integrand

      integrand = &
        (2*(336*Cos(4*y)**2*Cos(4*z)**2*Sin(4*x)**2 + &
      420*Cos(2*y)**2*Cos(4*z)**2*Sin(8*x)*Sin(2*y) + &
      84*Cos(2*z)**2*Sin(2*x)**2*Sin(2*y)**2 + &
      240*Cos(4*z)*Sin(4*x)*Sin(4*y) + &
      448*Cos(4*z)**2*Sin(4*x)**2*Sin(4*y)**2 + &
      150*Cos(6*x)*Cos(2*z)**2*Sin(4*y)*Sin(2*z) + &
      240*Sin(4*x)*Sin(4*y)*Sin(4*z) + &
      448*Sin(4*x)**2*Sin(4*y)**2*Sin(4*z)**2 + &
      15/Sin(2*x)/Sin(2*z)*Sin(4*x)**2*Sin(8*y)*Sin(4*z)**2 - &
      42*Cos(4*x)*(45*Cos(2*y)*Cos(4*z) - &
      8*Cos(2*z)*Sin(2*x)*Sin(2*y)*Sin(4*y)*Sin(4*z)) + &
      336*Cos(4*x)**2*(Sin(4*y)**2*Sin(4*z)**2 + &
      Cos(4*y)**2*(Cos(4*z) + Sin(4*z))**2) + &
      60*Cos(2*y)**2*Sin(8*x)*Sin(2*y)*Sin(8*z) + &
      (75*Cos(4*y)*Sin(8*x)*Sin(2*y)*Sin(8*z))/2. - &
      224*Sin(4*x)**2*Sin(4*y)**2*Sin(8*z) + &
      Cos(2*x)**2*(84*Cos(2*y)**2*Cos(2*z)**2 + &
      60*Cos(2*z)*Sin(2*x)*Sin(8*y) + &
         60*Cos(6*z)*Sin(2*x)*Sin(8*y) + &
         75*Cos(4*x)*Cos(2*z)**3*Sin(8*y)*Sin(2*z) + &
         112*Sin(2*y)**2*Sin(2*z)**2 - &
         15*Sin(2*y)*Sin(4*z) + 15*Sin(6*y)*Sin(4*z) + &
         (15*Cos(4*x)*Sin(8*y)*Sin(8*z))/8.) + &
      Cos(2*x)*Sin(2*y)*Sin(2*z)*(120 + &
      60*Cos(2*y)*(5*Cos(2*z)**2 + Cos(4*x)*Cos(4*z)) + &
         75*Cos(4*y)**2*Cos(2*z)**2*Cos(4*z)*Sin(8*x) + &
         30*Cos(2*y)**2*Cos(4*y)*Cos(4*z)*(4*Cos(2*z)**2 + &
         7*Cos(4*z)/Tan(2*z))*Sin(8*x) - &
         224*Cos(4*z)*Sin(4*x)*Sin(4*y) - 224*Sin(4*x)*Sin(4*y)*Sin(4*z) + &
         21/(Sin(2*y)**2)/(Sin(2*z)**2)*Sin(4*x)*Sin(8*y)*Sin(8*z))))/ &
      (75.*(2 + Cos(2*x)*Cos(4*y)*Cos(2*z)))

   end function g_energy

end module mms

! vim: set ts=3 sw=3 sts=3:
