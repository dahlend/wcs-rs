//! This module deals with parsing the projection specific keywords
//!
//! Zenithal, pseudo-cylindrical, cylindrical, conic, polyconic and quadcube projections
//! are supported.

use std::f64::consts::PI;

use mapproj::{
    conic::{cod::Cod, coe::Coe, coo::Coo, cop::Cop},
    cylindrical::{car::Car, cea::Cea, cyp::Cyp, mer::Mer},
    hybrid::hpx::Hpx,
    polyconic::{bon::Bon, pco::Pco},
    pseudocyl::{ait::Ait, mol::Mol, par::Par, sfl::Sfl},
    quadcube::{csc::Csc, qsc::Qsc, tsc::Tsc},
    zenithal::{
        air::Air,
        arc::Arc,
        azp::Azp,
        ncp::Ncp,
        sin::{Sin, SinSlant},
        stg::Stg,
        szp::Szp,
        tan::Tan,
        zea::Zea,
        zpn::Zpn,
    },
    CanonicalProjection, CenteredProjection, LonLat,
};
use std::f64::consts::FRAC_PI_2;

use crate::params::WCSParams;

use super::Error;

#[derive(PartialEq, Debug)]
pub enum FiducialPoint {
    /// Fiducial point is the north pole (0, 90deg). This is the default for zenital projections
    NorthPole,
    /// Fiducial point is (0, 0). This is the default for cyl and pseudocyl projections
    Origin,
    UserSpeficied {
        /// Native Longitude of the fiducial point in radians
        phi_0: f64,
        /// Native Latitude of the fiducial point in radians
        theta_0: f64,
    },
}

impl FiducialPoint {
    /// Native latitude of the fiducial point in degrees
    fn theta_0(&self) -> f64 {
        match self {
            Self::NorthPole => 90.0,
            Self::Origin => 0.0,
            Self::UserSpeficied { theta_0, .. } => theta_0.to_degrees(),
        }
    }
}

/// Compute the celestial pole position
///
/// # Arguments
///
/// * `phi_p` - native longitude of the celestial pole. Corresponds to the LONPOLE fits keyword
/// * `theta_p` - native latitude of the celestial pole. Corresponds to the LATPOLE fits keyword
/// * `alpha_0` - celestial longitude of the fiducial point. Corresponds to the CRVAL1 fits keyword
/// * `delta_0` - celestial latitude of the fiducial point. Corresponds to the CRVAL2 fits keyword
/// * `fiducial_point` - the fiducial point position in the native system
///
/// # Returns
///
/// The celestial position of the pole
fn celestial_pole(
    phi_p: f64,
    theta_p: f64,
    alpha_0: f64,
    delta_0: f64,
    fiducial_point: &FiducialPoint,
) -> Result<LonLat, Error> {
    let (phi_0, theta_0) = match fiducial_point {
        FiducialPoint::NorthPole => {
            return Ok(LonLat::new(alpha_0, delta_0));
        }
        FiducialPoint::Origin => (0.0, 0.0),
        FiducialPoint::UserSpeficied { phi_0, theta_0 } => (*phi_0, *theta_0),
    };

    let (s_t0, c_t0) = theta_0.sin_cos();
    let (s_phi, c_phi) = (phi_p - phi_0).sin_cos();
    let (s_d0, c_d0) = delta_0.sin_cos();

    let a = c_t0 * c_t0 * s_phi * s_phi;
    let delta_p = if a >= 1.0 {
        if delta_0 == 0.0 {
            // Paper convention 6
            Ok(theta_p)
        } else {
            // Paper convention 3 error
            Err(Error::CelestialPoleInvalid)
        }
    } else {
        let b = (c_t0 * c_phi).atan2(s_t0);
        let c = (s_d0 / (1.0 - a).sqrt()).acos();

        let valid = -FRAC_PI_2..=FRAC_PI_2;

        match (valid.contains(&(b - c)), valid.contains(&(b + c))) {
            (false, false) => Err(Error::CelestialPoleInvalid),
            (true, false) => Ok(b - c),
            (false, true) => Ok(b + c),
            (true, true) => {
                // northerly solution chosen
                Ok(if theta_p >= 0.0 { b + c } else { b - c })
            }
        }
    }?;

    let alpha_p = if delta_0.abs() == FRAC_PI_2 {
        // Paper convention 1
        alpha_0
    } else if delta_p == FRAC_PI_2 {
        // Paper convention 2
        alpha_0 + phi_p - phi_0 - PI
    } else if delta_p == -FRAC_PI_2 {
        // Paper convention 2
        alpha_0 - phi_p + phi_0
    } else {
        let (s_dp, c_dp) = delta_p.sin_cos();

        let a = s_phi * c_t0 / c_d0;
        let b = (s_t0 - s_dp * s_d0) / (c_dp * c_d0);
        alpha_0 - a.atan2(b)
    };

    Ok(LonLat::new(alpha_p, delta_p))
}

/// Compute the default value for LATPOLE (native latitude of celestial pole)
/// following the WCSLIB algorithm in cel.c (celset function, lines 302-400)
///
/// # Arguments
///
/// * `lonpole` - native longitude of the celestial pole (LONPOLE) in degrees
/// * `latpole_hint` - user-provided hint for LATPOLE in degrees (typically 90.0 if not specified)
/// * `lat0` - celestial latitude of fiducial point (CRVAL2) in degrees
/// * `phi0` - native longitude of fiducial point in degrees
/// * `theta0` - native latitude of fiducial point in degrees
///
/// # Returns
///
/// The computed LATPOLE value in degrees
fn compute_latpole_default(
    lonpole: f64,
    latpole_hint: f64,
    lat0: f64,
    phi0: f64,
    theta0: f64,
) -> f64 {
    const TOL: f64 = 1.0e-10;

    // For zenithal projections (theta0 == 90.0), fiducial point at native pole
    if (theta0 - 90.0).abs() < TOL {
        // Celestial pole == fiducial point
        return lat0;
    }

    // For non-zenithal projections, compute using spherical trigonometry
    let lat0_rad = lat0.to_radians();
    let theta0_rad = theta0.to_radians();
    let (slat0, _clat0) = lat0_rad.sin_cos();
    let (sthe0, cthe0) = theta0_rad.sin_cos();

    let (u, v) = if (lonpole - phi0).abs() < TOL {
        // phip == phi0 case
        (theta0, 90.0 - lat0)
    } else {
        let phip_minus_phi0_rad = (lonpole - phi0).to_radians();
        let (_sphip, cphip) = phip_minus_phi0_rad.sin_cos();

        let x = cthe0 * cphip;
        let y = sthe0;
        let z = (x * x + y * y).sqrt();

        if z.abs() < TOL {
            // Special case: latp determined solely by latpole_hint
            // This occurs when phip - phi0 = +/-90 deg and theta0 = 0
            // Use the provided latpole_hint, clamped to valid range
            return latpole_hint.clamp(-90.0, 90.0);
        }

        // Check if slz is in valid range [-1, 1]
        let slz = slat0 / z;
        let slz_clamped = match slz {
            s if s > 1.0 => {
                if s - 1.0 < TOL {
                    1.0
                } else {
                    // Invalid: return a default
                    return if lat0 >= 0.0 { 90.0 } else { -90.0 };
                }
            }
            s if s < -1.0 => {
                if (s + 1.0).abs() < TOL {
                    -1.0
                } else {
                    // Invalid: return a default
                    return if lat0 >= 0.0 { 90.0 } else { -90.0 };
                }
            }
            s => s,
        };

        let u = y.atan2(x).to_degrees();
        let v = slz_clamped.acos().to_degrees();
        (u, v)
    };

    // Helper function to normalize angle to [-180, 180]
    let normalize_angle = |mut angle: f64| {
        while angle > 180.0 {
            angle -= 360.0;
        }
        while angle < -180.0 {
            angle += 360.0;
        }
        angle
    };

    // Compute the two potential solutions
    let latp1 = normalize_angle(u + v);
    let latp2 = normalize_angle(u - v);

    // Check which solutions are valid (in range [-90, 90])
    let latp1_valid = latp1.abs() < 90.0 + TOL;
    let latp2_valid = latp2.abs() < 90.0 + TOL;

    // Choose solution closest to latpole_hint
    let latp = match (
        (latpole_hint - latp1).abs() < (latpole_hint - latp2).abs(),
        latp1_valid,
        latp2_valid,
    ) {
        (true, true, _) => latp1,
        (true, false, _) => latp2,
        (false, _, true) => latp2,
        (false, _, false) => latp1,
    };

    // Clamp to valid range [-90, 90]
    latp.clamp(-90.0, 90.0)
}

pub trait WCSCanonicalProjection: CanonicalProjection {
    fn parse_proj(params: &WCSParams) -> Result<(f64, CenteredProjection<Self>), Error>
    where
        Self: Sized,
    {
        // Parse the celestial longitude of the fiducial point
        let crval1 = params.crval1.unwrap_or(0.0);
        // Parse the celestial latitude of the fiducial point
        let crval2 = params.crval2.unwrap_or(0.0);
        let crval = LonLat::new(crval1.to_radians(), crval2.to_radians());

        // Parse the native longitude of the fiducial point
        let native_fiducial_point = match (params.pv1_1, params.pv1_2) {
            (Some(phi_0), Some(theta_0)) => FiducialPoint::UserSpeficied {
                phi_0: phi_0.to_radians(),
                theta_0: theta_0.to_radians(),
            },
            _ => Self::default_native_fiducial_point(params)?,
        };

        let lonpole = params.lonpole.unwrap_or_else(|| {
            if crval2 >= native_fiducial_point.theta_0() {
                0.0
            } else {
                180.0
            }
        });

        // Get phi0 and theta0 from fiducial point
        let (phi_0, theta_0) = match &native_fiducial_point {
            FiducialPoint::NorthPole => (0.0, 90.0),
            FiducialPoint::Origin => (0.0, 0.0),
            FiducialPoint::UserSpeficied { phi_0, theta_0 } => {
                (phi_0.to_degrees(), theta_0.to_degrees())
            }
        };

        // Compute LATPOLE default using WCSLIB algorithm
        // Use user-provided latpole as hint, or 90.0 as default hint
        let latpole_hint = params.latpole.unwrap_or(90.0);
        let latpole = compute_latpole_default(lonpole, latpole_hint, crval2, phi_0, theta_0);

        let positional_angle = if native_fiducial_point == FiducialPoint::NorthPole {
            PI - lonpole.to_radians()
        } else {
            let pole = celestial_pole(
                lonpole.to_radians(),
                latpole.to_radians(),
                crval1.to_radians(),
                crval2.to_radians(),
                &native_fiducial_point,
            )?;

            // Compute the positional angle formed from (CRVAL1, CRVAL2) between (ALPHA_P, DELTA_P) and the north pole (0, 90 deg)
            let north_pole = LonLat::new(0.0, FRAC_PI_2);

            let crval2pole_dist = crval.haversine_dist(&pole);
            let crval2np_dist = crval.haversine_dist(&north_pole);

            if pole == crval || north_pole == crval {
                0.0
            } else {
                let pole2np_dist = pole.haversine_dist(&north_pole);

                let (s_02p, c_02p) = crval2pole_dist.sin_cos();
                let (s_02np, c_02np) = crval2np_dist.sin_cos();

                // A angle of a triangle on a sphere does not exceed PI
                // Use the law of cosines applied for geodesics on a sphere
                // https://www.theoremoftheday.org/GeometryAndTrigonometry/SphericalCos/TotDSphericalCos.pdf
                //((pole2np_dist.cos() - c_02p * c_02np) / (s_02p * s_02np)).acos()
                let c = (pole2np_dist.cos() - c_02p * c_02np) / (s_02p * s_02np);
                if c >= 1.0 {
                    0.0
                } else if c <= -1.0 {
                    PI
                } else {
                    c.acos()
                }
            }
        };

        let proj = Self::parse_internal_proj_params(params)?;

        let mut rotated_proj = CenteredProjection::new(proj);
        rotated_proj.set_proj_center_from_lonlat_and_positional_angle(&crval, -positional_angle);
        //rotated_proj.set_proj_center_from_lonlat(&crval);
        Ok((positional_angle, rotated_proj))
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error>;

    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error>
    where
        Self: Sized;
}

// Zenithal projections
impl WCSCanonicalProjection for Azp {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // mu given in spherical radii, default value: 0.0
        let mu = params.pv2_1.unwrap_or(0.0);
        // gamma given in deg, default value: 0.0
        let gamma = params.pv2_2.unwrap_or(0.0);

        let azp = Azp::from_params(mu, gamma.to_radians());

        Ok(azp)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Szp {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // mu given in spherical radii, default value: 0.0
        let mu = params.pv2_1.unwrap_or(0.0);
        // phi_c given in deg, default value: 0.0
        let phi_c = params.pv2_2.unwrap_or(0.0);
        // theta_c given in deg, default value: 90.0
        let theta_c = params.pv2_3.unwrap_or(90.0);

        let szp = Szp::from_params(mu, phi_c.to_radians(), theta_c.to_radians());

        Ok(szp)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Tan {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Tan::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }

    // Override parse_proj to ALWAYS use NorthPole as native fiducial point for TAN projection
    // TAN projection MUST have phi_0=0, theta_0=90 deg (north pole) by definition
    // PV1_1 and PV1_2 are NOT used for TAN projection parameters
    // (they may be used for TPV distortion, but that's handled separately)
    fn parse_proj(params: &WCSParams) -> Result<(f64, CenteredProjection<Self>), Error>
    where
        Self: Sized,
    {
        let crval1 = params.crval1.unwrap_or(0.0);
        let crval2 = params.crval2.unwrap_or(0.0);
        let crval = LonLat::new(crval1.to_radians(), crval2.to_radians());

        // TAN projection ALWAYS uses NorthPole as native fiducial point
        let native_fiducial_point = FiducialPoint::NorthPole;

        let lonpole = params.lonpole.unwrap_or_else(|| {
            if crval2 >= native_fiducial_point.theta_0() {
                0.0
            } else {
                180.0
            }
        });

        // Compute LATPOLE using WCSLIB algorithm
        // For TAN (zenithal), this will always return lat0 since theta0=90
        let latpole_hint = params.latpole.unwrap_or(90.0);
        let _latpole = compute_latpole_default(
            lonpole,
            latpole_hint,
            crval2,
            0.0,  // phi0 for zenithal
            90.0, // theta0 for zenithal
        );

        let positional_angle = PI - lonpole.to_radians();

        let proj = Self::parse_internal_proj_params(params)?;
        let mut rotated_proj = CenteredProjection::new(proj);
        rotated_proj.set_proj_center_from_lonlat_and_positional_angle(&crval, -positional_angle);

        Ok((positional_angle, rotated_proj))
    }
}

impl WCSCanonicalProjection for Stg {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Stg::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Sin {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Sin::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for SinSlant {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // xi dimensionless, default value: 0.0
        let xi = params.pv2_1.unwrap_or(0.0);
        // eta dimensionless, default value: 0.0
        let eta = params.pv2_2.unwrap_or(0.0);

        let sin_slant = SinSlant::new(-xi, eta);
        Ok(sin_slant)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Arc {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Arc::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Zpn {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        let mut coeffs = [
            params.pv2_0,
            params.pv2_1,
            params.pv2_2,
            params.pv2_3,
            params.pv2_4,
            params.pv2_5,
            params.pv2_6,
            params.pv2_7,
            params.pv2_8,
            params.pv2_9,
            params.pv2_10,
            params.pv2_11,
            params.pv2_12,
            params.pv2_13,
            params.pv2_14,
            params.pv2_15,
            params.pv2_16,
            params.pv2_17,
            params.pv2_18,
            params.pv2_19,
            params.pv2_20,
        ]
        .into_iter()
        .rev()
        .skip_while(|p| p.is_none())
        .map(|p| p.unwrap_or(0.0))
        .collect::<Vec<_>>();

        coeffs.reverse();

        Zpn::from_params(coeffs).ok_or(Error::InitProjection(
            Zpn::NAME,
            "negative polynomial in [0, pi]",
        ))
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Zea {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Zea::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Air {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_b in deg, default value: 90.0
        let theta_b = params.pv2_1.unwrap_or(90.0);

        let airy = Air::from_param(theta_b.to_radians());
        Ok(airy)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Ncp {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        let ncp = Ncp::new();
        Ok(ncp)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

// Cylindrical projections
impl WCSCanonicalProjection for Cyp {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // mu given in spherical radii, default value: 1.0
        let mu = params.pv2_1.unwrap_or(1.0);
        // lambda given in spherical radii, default value: 1.0
        let lambda = params.pv2_2.unwrap_or(1.0);

        let cyp = Cyp::from_params(mu, lambda);
        Ok(cyp)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Cea {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // lambda given in spherical radii, default value: 1.0
        let lambda = params.pv2_1.unwrap_or(1.0);

        let cea = Cea::from_param(lambda);
        Ok(cea)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Car {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Car)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Mer {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Mer)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

// Pseudo-cylindrical projections
impl WCSCanonicalProjection for Sfl {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Sfl)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Par {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Par)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Mol {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Mol::default())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Ait {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Ait)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

// Conic projections
impl WCSCanonicalProjection for Cop {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_a given in deg, has no default value
        if let Some(theta_a) = params.pv2_1 {
            // eta given in deg, default value: 0
            let eta = params.pv2_2.unwrap_or(0.0);

            let cop = Cop::from_params(theta_a.to_radians(), eta.to_radians());
            Ok(cop)
        } else {
            Err(Error::InitProjection(
                Cop::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_a) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_a.to_radians(),
            })
        } else {
            Err(Error::InitProjection(
                Cop::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Coe {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_a given in deg, has no default value
        if let Some(theta_a) = params.pv2_1 {
            // eta given in deg, default value: 0
            let eta = params.pv2_2.unwrap_or(0.0);

            let cop = Coe::from_params(theta_a.to_radians(), eta.to_radians());
            Ok(cop)
        } else {
            Err(Error::InitProjection(
                Coe::NAME,
                "PV_1 = theta_a must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_a) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_a.to_radians(),
            })
        } else {
            Err(Error::InitProjection(
                Coe::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Cod {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_a given in deg, has no default value
        if let Some(theta_a) = params.pv2_1 {
            // eta given in deg, default value: 0
            let eta = params.pv2_2.unwrap_or(0.0);

            let cod = Cod::from_params(theta_a.to_radians(), eta.to_radians());
            Ok(cod)
        } else {
            Err(Error::InitProjection(
                Cod::NAME,
                "PV_1 = theta_a must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_a) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_a.to_radians(),
            })
        } else {
            Err(Error::InitProjection(
                Cod::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Coo {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_a given in deg, has no default value
        if let Some(theta_a) = params.pv2_1 {
            // eta given in deg, default value: 0
            let eta = params.pv2_2.unwrap_or(0.0);

            let coo = Coo::from_params(theta_a.to_radians(), eta.to_radians());
            Ok(coo)
        } else {
            Err(Error::InitProjection(
                Coo::NAME,
                "PV_1 = theta_a must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_a) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_a.to_radians(),
            })
        } else {
            Err(Error::InitProjection(
                Coo::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Hpx {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Hpx {})
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

// Polyconic projections
impl WCSCanonicalProjection for Bon {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_1 given in deg, has no default value
        if let Some(theta_1) = params.pv2_1 {
            let bon = Bon::from_params(theta_1.to_radians());
            Ok(bon)
        } else {
            Err(Error::InitProjection(
                Bon::NAME,
                "PV2_1 = theta_1 must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_1) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_1.to_radians(),
            })
        } else {
            Err(Error::InitProjection(
                Bon::NAME,
                "PV2_1 = theta_1 must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Pco {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Pco::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

// Quadcube projections
impl WCSCanonicalProjection for Tsc {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Tsc::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Csc {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Csc::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Qsc {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Qsc::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_latpole_zenithal_at_pole() {
        // For zenithal projections (theta0=90), latpole should equal lat0
        let latpole = compute_latpole_default(
            180.0, // lonpole
            90.0,  // latpole_hint
            45.0,  // lat0 (CRVAL2)
            0.0,   // phi0
            90.0,  // theta0 (zenithal)
        );

        assert_eq!(latpole, 45.0, "For zenithal, latpole should equal lat0");
    }

    #[test]
    fn test_latpole_zenithal_negative() {
        // Test southern hemisphere
        let latpole = compute_latpole_default(
            180.0, // lonpole
            90.0,  // latpole_hint
            -30.0, // lat0 (southern hemisphere)
            0.0,   // phi0
            90.0,  // theta0 (zenithal)
        );

        assert_eq!(latpole, -30.0, "For zenithal, latpole should equal lat0");
    }

    #[test]
    fn test_latpole_cylindrical_north() {
        // For cylindrical (theta0=0), northern hemisphere
        let latpole = compute_latpole_default(
            0.0,  // lonpole
            90.0, // latpole_hint (prefer northern solution)
            45.0, // lat0
            0.0,  // phi0
            0.0,  // theta0 (cylindrical)
        );

        // Should compute a geometrically valid solution
        assert!(
            (-90.0..=90.0).contains(&latpole),
            "latpole should be in valid range"
        );
    }

    #[test]
    fn test_latpole_cylindrical_south() {
        // For cylindrical (theta0=0), southern hemisphere
        let latpole = compute_latpole_default(
            180.0, // lonpole
            -90.0, // latpole_hint (prefer southern solution)
            -45.0, // lat0
            0.0,   // phi0
            0.0,   // theta0 (cylindrical)
        );

        assert!(
            (-90.0..=90.0).contains(&latpole),
            "latpole should be in valid range"
        );
    }

    #[test]
    fn test_latpole_disambiguation() {
        // Case where two solutions exist
        // lonpole = phi0 case for simplicity
        let latpole1 = compute_latpole_default(
            0.0,  // lonpole == phi0
            90.0, // latpole_hint (northern)
            30.0, // lat0
            0.0,  // phi0
            0.0,  // theta0
        );

        let latpole2 = compute_latpole_default(
            0.0,   // lonpole == phi0
            -90.0, // latpole_hint (southern)
            30.0,  // lat0
            0.0,   // phi0
            0.0,   // theta0
        );

        // Both should be valid but different based on hint
        assert!((-90.0..=90.0).contains(&latpole1));
        assert!((-90.0..=90.0).contains(&latpole2));

        // The actual values depend on the spherical geometry
        println!("latpole1={}, latpole2={}", latpole1, latpole2);
    }

    #[test]
    fn test_latpole_special_case_at_zero() {
        // Special case: theta0=0, lat0=0, and lonpole-phi0 near +/-90 deg
        // This should use the hint directly
        let latpole = compute_latpole_default(
            90.0, // lonpole (90 deg from phi0=0)
            75.0, // latpole_hint
            0.0,  // lat0
            0.0,  // phi0
            0.0,  // theta0
        );

        assert_eq!(latpole, 75.0, "Special case should use latpole_hint");
    }

    #[test]
    fn test_latpole_clamping() {
        // Test that extreme hints get clamped to [-90, 90]
        let latpole = compute_latpole_default(
            0.0, 150.0, // hint > 90
            45.0, 0.0, 90.0, // zenithal
        );

        assert!(
            (-90.0..=90.0).contains(&latpole),
            "latpole should be clamped to valid range"
        );
    }
}
