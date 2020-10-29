//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//=========================================================================================
//! \file hlle.cpp
//  \brief HLLE Riemann solver for hydrodynamics
//
//  Computes 1D fluxes using the Harten-Lax-van Leer (HLL) Riemann solver.  This flux is
//  very diffusive, especially for contacts, and so it is not recommended for use in
//  applications.  However, as shown by Einfeldt et al.(1991), it is positively
//  conservative (cannot return negative densities or pressure), so it is a useful
//  option when other approximate solvers fail and/or when extra dissipation is needed.
//
// REFERENCES:
// - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//   Springer-Verlag, Berlin, (1999) chpt. 10.
// - Einfeldt et al., "On Godunov-type methods near low densities", JCP, 92, 273 (1991)
// - A. Harten, P. D. Lax and B. van Leer, "On upstream differencing and Godunov-type
//   schemes for hyperbolic conservation laws", SIAM Review 25, 35-61 (1983).

// C/C++ headers
#include <algorithm> // max(), min()
#include <cmath>     // sqrt()
#include <memory>    // shared_ptr

// Athena headers
#include "../../main.hpp"
#include "config.hpp"
#include "riemann.hpp"

using parthenon::MeshBlockVarFluxPack;
using parthenon::MeshBlockVarPack;
using parthenon::Real;
//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//  \brief The HLLE Riemann solver for hydrodynamics (both adiabatic and isothermal)

void RiemannSolver(const int kl, const int ku, const int jl, const int ju, const int il,
                   const int iu, const int ivx, const MeshBlockVarPack<Real> &wl_pack,
                   const MeshBlockVarPack<Real> &wr_pack,
                   MeshBlockVarFluxPack<Real> &cons_pack, const AdiabaticHydroEOS &peos) {

  int ivy = IVX + ((ivx - IVX) + 1) % 3;
  int ivz = IVX + ((ivx - IVX) + 2) % 3;
  Real gamma = peos.GetGamma();
  Real gm1 = gamma - 1.0;
  Real iso_cs = 1.0;
  bool NON_BAROTROPIC_EOS = true;

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Riemann hydro HLLE", parthenon::DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, kl, ku, jl, ju, il, iu,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &wl = wl_pack(b);
        const auto &wr = wr_pack(b);
        // auto &flx = cons_pack(b).flux(ivx);
        auto &cons = cons_pack(b);
        Real wli[(NHYDRO)], wri[(NHYDRO)], wroe[(NHYDRO)];
        Real fl[(NHYDRO)], fr[(NHYDRO)], flxi[(NHYDRO)];

        //--- Step 1.  Load L/R states into local variables

        wli[IDN] = wl(IDN, k, j, i);
        wli[IVX] = wl(ivx, k, j, i);
        wli[IVY] = wl(ivy, k, j, i);
        wli[IVZ] = wl(ivz, k, j, i);
        if (NON_BAROTROPIC_EOS) wli[IPR] = wl(IPR, k, j, i);

        wri[IDN] = wr(IDN, k, j, i);
        wri[IVX] = wr(ivx, k, j, i);
        wri[IVY] = wr(ivy, k, j, i);
        wri[IVZ] = wr(ivz, k, j, i);
        if (NON_BAROTROPIC_EOS) wri[IPR] = wr(IPR, k, j, i);

        //--- Step2.  Compute Roe-averaged state

        Real sqrtdl = std::sqrt(wli[IDN]);
        Real sqrtdr = std::sqrt(wri[IDN]);
        Real isdlpdr = 1.0 / (sqrtdl + sqrtdr);

        wroe[IDN] = sqrtdl * sqrtdr;
        wroe[IVX] = (sqrtdl * wli[IVX] + sqrtdr * wri[IVX]) * isdlpdr;
        wroe[IVY] = (sqrtdl * wli[IVY] + sqrtdr * wri[IVY]) * isdlpdr;
        wroe[IVZ] = (sqrtdl * wli[IVZ] + sqrtdr * wri[IVZ]) * isdlpdr;

        // Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
        // rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
        Real el, er, hroe;
        if (NON_BAROTROPIC_EOS) {
          el = wli[IPR] / gm1 +
               0.5 * wli[IDN] * (SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
          er = wri[IPR] / gm1 +
               0.5 * wri[IDN] * (SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
          hroe = ((el + wli[IPR]) / sqrtdl + (er + wri[IPR]) / sqrtdr) * isdlpdr;
        }

        //--- Step 3.  Compute sound speed in L,R, and Roe-averaged states

        Real cl = peos.SoundSpeed(wli);
        Real cr = peos.SoundSpeed(wri);
        Real a = iso_cs;
        if (NON_BAROTROPIC_EOS) {
          Real q = hroe - 0.5 * (SQR(wroe[IVX]) + SQR(wroe[IVY]) + SQR(wroe[IVZ]));
          a = (q < 0.0) ? 0.0 : std::sqrt(gm1 * q);
        }

        //--- Step 4. Compute the max/min wave speeds based on L/R and Roe-averaged values

        Real al = fmin((wroe[IVX] - a), (wli[IVX] - cl));
        Real ar = fmax((wroe[IVX] + a), (wri[IVX] + cr));

        Real bp = ar > 0.0 ? ar : 0.0;
        Real bm = al < 0.0 ? al : 0.0;

        //-- Step 5. Compute L/R fluxes along the lines bm/bp: F_L - (S_L)U_L; F_R -
        //(S_R)U_R

        Real vxl = wli[IVX] - bm;
        Real vxr = wri[IVX] - bp;

        fl[IDN] = wli[IDN] * vxl;
        fr[IDN] = wri[IDN] * vxr;

        fl[IVX] = wli[IDN] * wli[IVX] * vxl;
        fr[IVX] = wri[IDN] * wri[IVX] * vxr;

        fl[IVY] = wli[IDN] * wli[IVY] * vxl;
        fr[IVY] = wri[IDN] * wri[IVY] * vxr;

        fl[IVZ] = wli[IDN] * wli[IVZ] * vxl;
        fr[IVZ] = wri[IDN] * wri[IVZ] * vxr;

        if (NON_BAROTROPIC_EOS) {
          fl[IVX] += wli[IPR];
          fr[IVX] += wri[IPR];
          fl[IEN] = el * vxl + wli[IPR] * wli[IVX];
          fr[IEN] = er * vxr + wri[IPR] * wri[IVX];
        } else {
          fl[IVX] += (iso_cs * iso_cs) * wli[IDN];
          fr[IVX] += (iso_cs * iso_cs) * wri[IDN];
        }

        //--- Step 6. Compute the HLLE flux at interface.

        Real tmp = 0.0;
        if (bp != bm) tmp = 0.5 * (bp + bm) / (bp - bm);

        flxi[IDN] = 0.5 * (fl[IDN] + fr[IDN]) + (fl[IDN] - fr[IDN]) * tmp;
        flxi[IVX] = 0.5 * (fl[IVX] + fr[IVX]) + (fl[IVX] - fr[IVX]) * tmp;
        flxi[IVY] = 0.5 * (fl[IVY] + fr[IVY]) + (fl[IVY] - fr[IVY]) * tmp;
        flxi[IVZ] = 0.5 * (fl[IVZ] + fr[IVZ]) + (fl[IVZ] - fr[IVZ]) * tmp;
        if (NON_BAROTROPIC_EOS)
          flxi[IEN] = 0.5 * (fl[IEN] + fr[IEN]) + (fl[IEN] - fr[IEN]) * tmp;

        // flx(IDN, k, j, i) = flxi[IDN];
        // flx(ivx, k, j, i) = flxi[IVX];
        // flx(ivy, k, j, i) = flxi[IVY];
        // flx(ivz, k, j, i) = flxi[IVZ];
        // if (NON_BAROTROPIC_EOS) flx(IEN, k, j, i) = flxi[IEN];
        cons.flux(ivx, IDN, k, j, i) = flxi[IDN];
        cons.flux(ivx, ivx, k, j, i) = flxi[IVX];
        cons.flux(ivx, ivy, k, j, i) = flxi[IVY];
        cons.flux(ivx, ivz, k, j, i) = flxi[IVZ];
        if (NON_BAROTROPIC_EOS) cons.flux(ivx, IEN, k, j, i) = flxi[IEN];
      });
}