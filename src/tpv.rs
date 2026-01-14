//! This module is an implementation of the TPV standard
//!
//! For the TPV convention, see
//! <https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html>
//! and related papers:
//! "Representations of distortions in FITS world coordinate systems" by Calabretta et al. (2004)
//! and "Representations of celestial coordinates in FITS" by Calabretta & Greisen (2002).
//!
//! TPV is a distortion convention used primarily with TAN projection that uses PV1_* and PV2_*
//! keywords to define polynomial distortion corrections applied AFTER the linear CD matrix
//! transformation. The evaluation steps are:
//!
//! 1. Apply the CD matrix to pixel coordinates to get intermediate world coordinates (xi, eta) in degrees
//! 2. Apply the polynomial distortion using PV coefficients to get corrected coordinates (xi', eta')
//! 3. Apply the tangent plane projection to (xi', eta') to get celestial coordinates

use mapproj::tpv::{Tpv, TpvCoeff, TpvPV};

use crate::error::Error;

use crate::params::WCSParams;

use paste::paste;

/// A macro that returns TPV coefficients for a given axis
macro_rules! combi_tpv_coeff {
    ($params:ident, $axis:tt, $( $idx:literal ),+ ) => {
        paste! {
            vec![
                $(
                    $params.[< pv $axis _ $idx >].unwrap_or(0.0)
                ),*
            ]
        }
    };
}

/// Parse TPV distortion parameters from WCS FITS header
///
/// # Arguments
///
/// * `params` - WCS parameters containing PV keywords
///
/// # Returns
///
/// A `Tpv` structure containing the polynomial distortion coefficients
///
/// # Notes
///
/// As of the latest mapproj update, range checks have been removed from the TPV
/// implementation, so image dimensions and reference pixels are no longer needed.
pub fn parse_tpv(params: &WCSParams) -> Result<Tpv, Error> {
    // According to TPV standard:
    // - PV1_1 and PV2_1 default to 1.0 (not 0.0)
    // - All other PV coefficients default to 0.0
    // - Coefficients go from PV*_0 to PV*_39 (max order 7)

    // Extract PV1 coefficients (for xi/longitude axis)
    let pv1_coeffs = combi_tpv_coeff!(
        params, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
        22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39
    );

    // Extract PV2 coefficients (for eta/latitude axis)
    let pv2_coeffs = combi_tpv_coeff!(
        params, 2, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
        22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39
    );

    // Apply the special defaults for PV1_1 and PV2_1
    // Note: The TpvCoeff::new_axis1/new_axis2 functions already handle the default
    // of PV1_1 and PV2_1 = 1.0, but we need to respect if they're explicitly set
    let mut pv1 = pv1_coeffs;
    let mut pv2 = pv2_coeffs;

    // PV1_1 defaults to 1.0 if not specified
    if params.pv1_1.is_none() && pv1.len() > 1 {
        pv1[1] = 1.0;
    }

    // PV2_1 defaults to 1.0 if not specified
    if params.pv2_1.is_none() && pv2.len() > 1 {
        pv2[1] = 1.0;
    }

    let pv1_coeff = TpvCoeff::new_axis1(&pv1);
    let pv2_coeff = TpvCoeff::new_axis2(&pv2);

    // Create TpvPV structure combining both axes
    let tpv_pv = TpvPV::new(pv1_coeff, pv2_coeff);

    Ok(Tpv::new(tpv_pv))
}
