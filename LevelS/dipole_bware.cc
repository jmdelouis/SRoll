#include <assert.h>

#include "dipole.h"
#include "pointing.h"
#include "vec3.h"
#include "sat_info.h"
#include "arr.h"
#include "focalplane_db.h"
#include "lsconstants.h"
#include "paramfile.h"


using namespace std;

namespace {

// DiporbLib constants
const double diplib_speedOfLight = 2.99792458e8;
const double diplib_tcmb         = 2.72548; /* As stated in D. J. Fixsen 2009, ApJ 709 (2009) (see arXiv:0911.1955) */

// SOLAR DIPOLE VALUES
// see http://pmwiki.sylvainmottet.fr/index.php?n=Main.DipoleValues

#if 0
// original (2013?) LevelS solar dipole value
// https://lambda.gsfc.nasa.gov/product/map/dr3/pub_papers/fiveyear/basic_results/wmap5basic_reprint.pdf
const char *soldipname = "LevelS";
//! Colatitude of the solar system motion relative to CMB
//! (ecliptical coordinates) (Hinshaw et al. 2009, ApJS, 180, 225).
const double solsysdir_ecl_theta = 1.765248346;

//! Longitude of the solar system motion relative to CMB
//! (ecliptical coordinates) (Hinshaw et al. 2009, ApJS, 180, 225).
const double solsysdir_ecl_phi = 2.995840906;

//! Speed of the solar system motion relative to CMB in m/s
//! (Hinshaw et al. 2009, ApJS, 180, 225).
const double solsysspeed = 369000.0;
#endif

#if 0
// DiporbLib v00-00-26 aka Planck 2015 values
// used for Planck2018/FFP10 simulations
const char  *soldipname          = "Planck2015";
const double diplib_DTsol        = 0.0033645; /* updated with planck dipole parameters */
const double solsysdir_ecl_theta = 1.7656131194952;
const double solsysdir_ecl_phi   = 2.9958896005736;
const double solsysspeed         = diplib_speedOfLight * diplib_DTsol / diplib_tcmb;
#endif

#if 0
// HFI 2017: GLon = 264.021 deg., GLat = 48.253 deg., amplitude = 3362.71 microK
// used for JAN18 SRoll2 simulations
const char  *soldipname          = "HFI2017";
const double diplib_DTsol        = 0.00336271;
const double solsysdir_ecl_theta = 1.7654364921748555;
const double solsysdir_ecl_phi   = 2.9961777544330825;
const double solsysspeed         = diplib_speedOfLight * diplib_DTsol / diplib_tcmb;
#endif

#if 1
// HFI 2018: GLon = 264.021 deg., GLat = 48.253 deg., amplitude = 3362.08 microK
// https://arxiv.org/abs/1807.06207, chapter 4.2.2, eq. 10, page 25
// used for APR20 SRoll3 simulations
const char  *soldipname          = "HFI2018";
const double diplib_DTsol        = 0.00336208;
const double solsysdir_ecl_theta = 1.7654364921748555;
const double solsysdir_ecl_phi   = 2.9961777544330825;
const double solsysspeed         = diplib_speedOfLight * diplib_DTsol / diplib_tcmb;
#endif

inline double T_ant (double t_thermo, double hnydk)
  { return hnydk/(exp(hnydk/t_thermo)-1); }

} // unnamed namespace


/*----------------------------------------------------------------------------*/

Dipole::Dipole (Sat_Info &info, focalplane_db &focal, paramfile &params,
                bool source_dipole, bool source_fsldp)
  : satinfo(info), fpdb(focal), do_dipole(source_dipole),
    do_fsldp(source_fsldp)
  {
  string speedstr = params.find<string>("dipole_speed","TOTAL");
  speed = TOTAL;
  if (speedstr=="TOTAL") speed = TOTAL;
  else if (speedstr=="SOLSYS") speed = SOLSYS;
  else if (speedstr=="SATELLITE") speed = SATELLITE;
  else planck_fail ("Incorrect value '" + speedstr + "' for dipole_speed");
  thermotemp = params.find<bool>("dipole_thermotemp",false);
  outputtype = params.find<int>("dipole_type",2);
  planck_assert ((outputtype>=1) && (outputtype<=8), "dipole: wrong type");
  solsysdir_v = pointing(solsysdir_ecl_theta, solsysdir_ecl_phi).to_vec3();
  dip_norm = params.find<double>("dip_norm",1);
  cout << "cxxmod/dipole.cc: using " << soldipname << " solar dipole value (" << (solsysspeed * diplib_tcmb / diplib_speedOfLight * 1e6) << " microK)" << endl;
  }


/*----------------------------------------------------------------------------*/

void Dipole::Add_Intensities (const string &det_id, const arr<vec3> &vdetpt,
  const arr<vec3> &vsldppt, arr<double> &intensity) const
  {
  planck_assert (multiequal(vdetpt.size(),intensity.size(),vsldppt.size()),
    "array size mismatch");

  if (do_fsldp)
    {
    double detfreq = fpdb.getValue<double>(det_id,"nu_cen");
    double fact = 0.5 * (1 + fpdb.getValue<double>(det_id,"epsilon"))*dip_norm;
    double hnydk = hPlanck*detfreq/kBoltzmann;
    double x = hnydk/tcmb;
    double expx = exp(x);
    double planck = x/(expx-1);
    double thermo2ant = expx *planck*planck;
    vec3 vsat = satinfo.velocity();

    vec3 detspeed;
    if (speed==TOTAL)     detspeed = solsysdir_v*solsysspeed + vsat;
    if (speed==SOLSYS)    detspeed = solsysdir_v*solsysspeed;
    if (speed==SATELLITE) detspeed = vsat;
    vec3 speed_dir = detspeed; speed_dir.Normalize();
    double detspeed_value = detspeed.Length();

    double beta = detspeed_value/speedOfLight;

#pragma omp parallel
{
    int m, sz=vdetpt.size();
#pragma omp for schedule(static)
    for (m=0; m<sz; ++m)
      {
      /* Start by adding far sidelobe dipole since this is always a plain
         dipole irrespective of the outputtype switch */
      double fsldp_thermo_dp = tcmb * beta * dotprod(speed_dir,vsldppt[m]);
      intensity[m] += fact *(thermotemp ?
        fsldp_thermo_dp :
        thermo2ant * fsldp_thermo_dp);
      }
}
    }

  /* now add contribution through main beam */
  Add_Intensities (det_id, vdetpt, intensity);
  }


/*----------------------------------------------------------------------------*/

void Dipole::Add_Intensities (const string &det_id, const arr<vec3> &vdetpt,
  arr<double> &intensity) const
  {
  planck_assert (vdetpt.size()==intensity.size(), "array size mismatch");

  if (!do_dipole) return;

  double detfreq = fpdb.getValue<double>(det_id,"nu_cen");
  double fact = 0.5 * (1 + fpdb.getValue<double>(det_id,"epsilon"))*dip_norm;
  double hnydk = hPlanck*detfreq/kBoltzmann;
  double x = hnydk/tcmb;
  double expx = exp(x);
  double planck = x/(expx-1);
  double thermo2ant = expx *planck*planck;
  double tcmbant = T_ant(tcmb, hnydk);

  vec3 vsat = satinfo.velocity();

  vec3 detspeed;
  if (speed==TOTAL)     detspeed = solsysdir_v*solsysspeed + vsat;
  if (speed==SOLSYS)    detspeed = solsysdir_v*solsysspeed;
  if (speed==SATELLITE) detspeed = vsat;
  vec3 speed_dir = detspeed; speed_dir.Normalize();
  double detspeed_value = detspeed.Length();

  double beta = detspeed_value/speedOfLight;
  double gamma = 1/sqrt(1-beta*beta);

#pragma omp parallel
{
  int m, sz=vdetpt.size();
#pragma omp for schedule(static)
  for (m=0; m<sz; ++m)
    {
    double cosdir = dotprod(speed_dir,vdetpt[m]);

    switch (outputtype)
      {
      case 1:
        // multimod_par.txt:
        // "1: total relativistic Doppler effect (monopole plus relativistic Doppler effect / "dipole anisotropy")"
        {
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        intensity[m] += fact * (thermotemp ?
          t_thermo_all :
          T_ant(t_thermo_all,hnydk));
        }
        break;

      case 2:
        // multimod_par.txt:
        // "2: relativistic Doppler effect ("dipole anisotropy")"
        {
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        intensity[m] += fact * (thermotemp ?
          (t_thermo_all-tcmb) :
          (T_ant(t_thermo_all,hnydk)-tcmbant));
        }
        break;

      case 3:
        // multimod_par.txt:
        // "3: total non-relativistic Doppler effect (pure dipole plus monopole)"
        {
        double t_thermo_mp_dp = tcmb * (1+beta*cosdir);
        intensity[m] += fact * (thermotemp ?
          t_thermo_mp_dp :
          T_ant(t_thermo_mp_dp, hnydk));
        }
        break;

      case 4:
        // multimod_par.txt:
        // "4: non-relativistic Doppler effect (pure dipole)"
        {
        double t_thermo_mp_dp = tcmb * (1+beta*cosdir);
        intensity[m] += fact * (thermotemp ?
          (t_thermo_mp_dp-tcmb) :
          (T_ant(t_thermo_mp_dp, hnydk)-tcmbant));
        }
        break;

      case 5:
        // multimod_par.txt:
        // "5: quadrupole component of relativistic Doppler effect"
        {
        double t_thermo_mp_dp = tcmb * (1+beta*cosdir);
        double t_thermo_mp_dp_qp = tcmb *
          (1+beta*cosdir + beta*beta*(cosdir*cosdir-0.5));
        intensity[m] += fact * (thermotemp ?
          (t_thermo_mp_dp_qp-t_thermo_mp_dp) :
          (T_ant(t_thermo_mp_dp_qp,hnydk)-T_ant(t_thermo_mp_dp,hnydk)));
        }
        break;

      case 6:
        // multimod_par.txt:
        // "6: higher components of relativistic Doppler effect"
        {
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        double t_thermo_mp_dp_qp = tcmb *
          (1+beta*cosdir + beta*beta*(cosdir*cosdir-0.5));
        intensity[m] += fact * (thermotemp ?
          (t_thermo_all-t_thermo_mp_dp_qp) :
          (T_ant(t_thermo_all,hnydk)-T_ant(t_thermo_mp_dp_qp,hnydk)));
        }
        break;

      case 7:
        // multimod_par.txt:
        // "7: relativistic component of the Doppler effect (difference between relativistic and non-relativistic Doppler effect)"
        {
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        double t_thermo_mp_dp = tcmb * (1+beta*cosdir);
        intensity[m] += fact * (thermotemp ?
          (t_thermo_all-t_thermo_mp_dp) :
          (T_ant(t_thermo_all,hnydk)-T_ant(t_thermo_mp_dp,hnydk)));
        }
        break;

      case 8:
        // 'case 2' (relativistic Doppler effect ("dipole anisotropy")) + kinetic doppler quadrupole
/*

* Notari, A. & Quartin, M., On the proper kinetic quadrupole CMB removal and the quadrupole anomalies. 2015, https://arxiv.org/pdf/1504.02076.pdf
* Pagano, L., Calibration systematics in future CMB space Missions 2017, http://pmwiki.sylvainmottet.fr/uploads/Main/DipoleValues/paper_syst_v6.pdf

  The amplitude of the effect for the solar dipole only should be something like 
    beta^2*TCMB*q(x)
  with
    beta = dipole/monopole
  and 
    q(x) = x/2*(exp(x)+1)/(exp(x)-1)
  and
    x = h*nu/k/TCMB
  putting all in the same formula at 100ghz
    beta^2*TCMB*q(x) = (3362.08/2.72548e6)**2 * 2.72548e6 * 0.2460 = 1.0203 microK

*/
        {
        assert( thermotemp && "<dipole_thermotemp> must be 'T(rue)' to use <dipole_type> = 8\n");
        double t_thermo_all = tcmb/(gamma*(1-beta*cosdir));
        intensity[m] += fact * (T_ant(t_thermo_all,hnydk)-tcmbant)/thermo2ant;
        }
        break;
      }
    }
  }
}
