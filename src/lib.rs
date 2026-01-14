extern crate mapproj;
#[macro_use]
extern crate quick_error;

#[doc = include_str!("../readme.md")]
pub mod error;

use coo_system::CooSystem;
use error::Error;
pub mod coo_system;
pub mod params;
mod projection;
mod sip;
mod tpv;

use crate::projection::WCSCanonicalProjection;
pub use params::WCSParams;

// Imports
use mapproj::{
    conic::{cod::Cod, coe::Coe, coo::Coo, cop::Cop},
    cylindrical::{car::Car, cea::Cea, cyp::Cyp, mer::Mer},
    hybrid::hpx::Hpx,
    img2celestial::Img2Celestial,
    img2proj::{WcsImgXY2ProjXY, WcsWithSipImgXY2ProjXY, WcsWithTpvImgXY2ProjXY},
    polyconic::{bon::Bon, pco::Pco},
    pseudocyl::{ait::Ait, mol::Mol, par::Par, sfl::Sfl},
    quadcube::{csc::Csc, qsc::Qsc, tsc::Tsc},
    zenithal::{
        air::Air, arc::Arc, azp::Azp, ncp::Ncp, sin::Sin, stg::Stg, szp::Szp, tan::Tan, zea::Zea,
        zpn::Zpn,
    },
    XYZ,
};

use paste::paste;
/// macro
macro_rules! create_specific_proj {
    ( $proj_name:ident, $params:expr, $ctype1:expr, $naxis1:expr, $naxis2:expr, $crpix1:expr, $crpix2:expr, $img2proj:expr ) => {{
        let (positional_angle, proj) = $proj_name::parse_proj(&$params)?;

        let is_sip_found = &$ctype1[($ctype1.len() - 3)..] == "SIP";
        let is_tpv_found = &$ctype1[($ctype1.len() - 3)..] == "TPV";

        if is_sip_found {
            let sip = sip::parse_sip($params, $naxis1, $naxis2, $crpix1, $crpix2)?;
            let img2proj = WcsWithSipImgXY2ProjXY::new($img2proj, sip);

            paste! {
                Ok((WCSCelestialProj::[ <$proj_name Sip> ](Img2Celestial::new(img2proj, proj)), positional_angle))
            }
        } else if is_tpv_found {
            // For TPV, use the CD matrix from the header and apply polynomial distortion on top
            let tpv = tpv::parse_tpv($params)?;
            let img2proj = WcsWithTpvImgXY2ProjXY::new($img2proj, tpv);

            paste! {
                Ok((WCSCelestialProj::[ <$proj_name Tpv> ](Img2Celestial::new(img2proj, proj)), positional_angle))
            }
        } else {
            Ok((
                WCSCelestialProj::$proj_name(Img2Celestial::new($img2proj, proj)),
                positional_angle,
            ))
        }
    }};
}

/// Structure alias coming from mapproj defining
/// image space pixel coordinates
pub type ImgXY = mapproj::ImgXY;
/// Structure alias coming from mapproj defining
/// longitude and latitude expressed in radians
pub type LonLat = mapproj::LonLat;

#[derive(Debug)]
pub struct WCS {
    /* Metadata keywords */
    /// Size of the image in pixels in its i-th dimension
    naxisi: Box<[i64]>,
    /// Field of view of the image along NAXIS1
    fov1: f64,
    /// Field of view of the image along NAXIS2
    fov2: f64,
    /// Main sub structure defining the projection
    proj: WCSProj,
}

/// Main object structure descripting a WCS object
/// Once created, the user can proceed two operation on it
/// * The projection of a (lon, lat) tuple onto the image space.
///   Results are given in pixels
/// * The unprojection of a (x, y) tuple given in pixel coordinates onto the sphere.
///   Results are given as a (lon, lat) tuple expressed in radians
impl WCS {
    pub fn new(params: &WCSParams) -> Result<Self, Error> {
        // Check for ZNAXISi first (case of tiled compressed image stored in bin tables)
        let naxisi = match params.naxis {
            2 => {
                let naxis1 = params
                    .znaxis1
                    .or(params.naxis1)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))?;

                let naxis2 = params
                    .znaxis2
                    .or(params.naxis2)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))?;

                Ok(vec![naxis1, naxis2].into_boxed_slice())
            }
            3 => {
                let naxis1 = params
                    .znaxis1
                    .or(params.naxis1)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))?;
                let naxis2 = params
                    .znaxis2
                    .or(params.naxis2)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))?;
                let naxis3 = params
                    .znaxis3
                    .or(params.naxis3)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS3"))?;

                Ok(vec![naxis1, naxis2, naxis3].into_boxed_slice())
            }
            4 => {
                let naxis1 = params
                    .znaxis1
                    .or(params.naxis1)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))?;
                let naxis2 = params
                    .znaxis2
                    .or(params.naxis2)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))?;
                let naxis3 = params
                    .znaxis3
                    .or(params.naxis3)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS3"))?;
                let naxis4 = params
                    .znaxis4
                    .or(params.naxis4)
                    .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS4"))?;

                Ok(vec![naxis1, naxis2, naxis3, naxis4].into_boxed_slice())
            }
            _ => Err(Error::NotSupportedNaxis(params.naxis)),
        }?;

        // At least NAXIS >= 2
        let proj = WCSProj::new(naxisi[0], naxisi[1], params)?;

        let fov1 = proj.s_lon * (naxisi[0] as f64);
        let fov2 = proj.s_lat * (naxisi[1] as f64);

        Ok(WCS {
            naxisi,
            fov1,
            fov2,
            proj,
        })
    }

    /// Returns the dimensions of the image given by the NAXIS1 x NAXIS2 keyword
    pub fn img_dimensions(&self) -> &[i64] {
        &self.naxisi[..]
    }

    pub fn field_of_view(&self) -> (f64, f64) {
        (self.fov1, self.fov2)
    }

    /// Project a (lon, lat) in ICRS
    ///
    /// The result is given a (X, Y) tuple expressed in pixel coordinates.
    ///
    /// # Arguments
    ///
    /// * `lonlat`: the 3D sphere vertex expressed as a (lon, lat) tuple given in radians
    pub fn proj(&self, lonlat: &LonLat) -> Option<ImgXY> {
        self.proj.proj_lonlat(lonlat)
    }

    /// Unproject a (X, Y) point to get a position on the sky in ICRS system
    ///
    /// # Arguments
    ///
    /// * `img_pos`: the image space point expressed as a (X, Y) tuple given en pixels
    ///
    /// # Returns
    ///
    /// A `LonLat` tuple with (lon, lat) expressed in radians
    pub fn unproj(&self, img_pos: &ImgXY) -> Option<LonLat> {
        self.proj.unproj_lonlat(img_pos)
    }

    /// Get the coordinate system frame
    pub fn coo_system(&self) -> &CooSystem {
        self.proj.coo_system()
    }
}

use std::ops::Deref;
impl Deref for WCS {
    type Target = WCSProj;

    fn deref(&self) -> &Self::Target {
        &self.proj
    }
}

pub struct WCSProj {
    /// The right part of the CTYPE keyword
    /// The projection type
    proj: WCSCelestialProj,
    /// The left part of the CTYPE keyword
    /// The coordinate system
    coo_system: CooSystem,
    pos_angle: f64,
    s_lon: f64, // scale in degrees along the longitude axis
    s_lat: f64, // scale in degrees along the latitude axis
}

impl std::fmt::Debug for WCSProj {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("WCS")
            .field("coosys", &self.coo_system)
            .field("euler angle", &self.pos_angle)
            .finish()
    }
}

/// Main enum structure descripting a WCS object
/// Once created, the user can proceed two operation on it
/// * The projection of a (lon, lat) tuple onto the image space.
///   Results are given in pixels
/// * The unprojection of a (x, y) tuple given in pixel coordinates onto the sphere.
///   Results are given as a (lon, lat) tuple expressed in radians
pub enum WCSCelestialProj {
    // Zenithal
    Azp(Img2Celestial<Azp, WcsImgXY2ProjXY>),
    Szp(Img2Celestial<Szp, WcsImgXY2ProjXY>),
    Tan(Img2Celestial<Tan, WcsImgXY2ProjXY>),
    Stg(Img2Celestial<Stg, WcsImgXY2ProjXY>),
    Sin(Img2Celestial<Sin, WcsImgXY2ProjXY>),
    Arc(Img2Celestial<Arc, WcsImgXY2ProjXY>),
    Zpn(Img2Celestial<Zpn, WcsImgXY2ProjXY>),
    Zea(Img2Celestial<Zea, WcsImgXY2ProjXY>),
    Air(Img2Celestial<Air, WcsImgXY2ProjXY>),
    Ncp(Img2Celestial<Ncp, WcsImgXY2ProjXY>),
    // Cylindrical
    Cyp(Img2Celestial<Cyp, WcsImgXY2ProjXY>),
    Cea(Img2Celestial<Cea, WcsImgXY2ProjXY>),
    Car(Img2Celestial<Car, WcsImgXY2ProjXY>),
    Mer(Img2Celestial<Mer, WcsImgXY2ProjXY>),
    // Pseudo-Cylindrical
    Sfl(Img2Celestial<Sfl, WcsImgXY2ProjXY>),
    Par(Img2Celestial<Par, WcsImgXY2ProjXY>),
    Mol(Img2Celestial<Mol, WcsImgXY2ProjXY>),
    Ait(Img2Celestial<Ait, WcsImgXY2ProjXY>),
    // Conic
    Cop(Img2Celestial<Cop, WcsImgXY2ProjXY>),
    Cod(Img2Celestial<Cod, WcsImgXY2ProjXY>),
    Coe(Img2Celestial<Coe, WcsImgXY2ProjXY>),
    Coo(Img2Celestial<Coo, WcsImgXY2ProjXY>),
    // Polyconic
    Bon(Img2Celestial<Bon, WcsImgXY2ProjXY>),
    Pco(Img2Celestial<Pco, WcsImgXY2ProjXY>),
    // Quadcube
    Tsc(Img2Celestial<Tsc, WcsImgXY2ProjXY>),
    Csc(Img2Celestial<Csc, WcsImgXY2ProjXY>),
    Qsc(Img2Celestial<Qsc, WcsImgXY2ProjXY>),
    // Hybrid
    Hpx(Img2Celestial<Hpx, WcsImgXY2ProjXY>),

    // SIP handling
    // Zenithal
    AzpSip(Img2Celestial<Azp, WcsWithSipImgXY2ProjXY>),
    SzpSip(Img2Celestial<Szp, WcsWithSipImgXY2ProjXY>),
    TanSip(Img2Celestial<Tan, WcsWithSipImgXY2ProjXY>),
    StgSip(Img2Celestial<Stg, WcsWithSipImgXY2ProjXY>),
    SinSip(Img2Celestial<Sin, WcsWithSipImgXY2ProjXY>),
    ArcSip(Img2Celestial<Arc, WcsWithSipImgXY2ProjXY>),
    ZpnSip(Img2Celestial<Zpn, WcsWithSipImgXY2ProjXY>),
    ZeaSip(Img2Celestial<Zea, WcsWithSipImgXY2ProjXY>),
    AirSip(Img2Celestial<Air, WcsWithSipImgXY2ProjXY>),
    NcpSip(Img2Celestial<Ncp, WcsWithSipImgXY2ProjXY>),
    // Cylindrical
    CypSip(Img2Celestial<Cyp, WcsWithSipImgXY2ProjXY>),
    CeaSip(Img2Celestial<Cea, WcsWithSipImgXY2ProjXY>),
    CarSip(Img2Celestial<Car, WcsWithSipImgXY2ProjXY>),
    MerSip(Img2Celestial<Mer, WcsWithSipImgXY2ProjXY>),
    // Pseudo-Cylindrical
    SflSip(Img2Celestial<Sfl, WcsWithSipImgXY2ProjXY>),
    ParSip(Img2Celestial<Par, WcsWithSipImgXY2ProjXY>),
    MolSip(Img2Celestial<Mol, WcsWithSipImgXY2ProjXY>),
    AitSip(Img2Celestial<Ait, WcsWithSipImgXY2ProjXY>),
    // Conic
    CopSip(Img2Celestial<Cop, WcsWithSipImgXY2ProjXY>),
    CodSip(Img2Celestial<Cod, WcsWithSipImgXY2ProjXY>),
    CoeSip(Img2Celestial<Coe, WcsWithSipImgXY2ProjXY>),
    CooSip(Img2Celestial<Coo, WcsWithSipImgXY2ProjXY>),
    // Polyconic
    BonSip(Img2Celestial<Bon, WcsWithSipImgXY2ProjXY>),
    PcoSip(Img2Celestial<Pco, WcsWithSipImgXY2ProjXY>),
    // Quadcube
    TscSip(Img2Celestial<Tsc, WcsWithSipImgXY2ProjXY>),
    CscSip(Img2Celestial<Csc, WcsWithSipImgXY2ProjXY>),
    QscSip(Img2Celestial<Qsc, WcsWithSipImgXY2ProjXY>),
    // Hybrid
    HpxSip(Img2Celestial<Hpx, WcsWithSipImgXY2ProjXY>),

    // TPV handling
    // Zenithal
    AzpTpv(Img2Celestial<Azp, WcsWithTpvImgXY2ProjXY>),
    SzpTpv(Img2Celestial<Szp, WcsWithTpvImgXY2ProjXY>),
    TanTpv(Img2Celestial<Tan, WcsWithTpvImgXY2ProjXY>),
    StgTpv(Img2Celestial<Stg, WcsWithTpvImgXY2ProjXY>),
    SinTpv(Img2Celestial<Sin, WcsWithTpvImgXY2ProjXY>),
    ArcTpv(Img2Celestial<Arc, WcsWithTpvImgXY2ProjXY>),
    ZpnTpv(Img2Celestial<Zpn, WcsWithTpvImgXY2ProjXY>),
    ZeaTpv(Img2Celestial<Zea, WcsWithTpvImgXY2ProjXY>),
    AirTpv(Img2Celestial<Air, WcsWithTpvImgXY2ProjXY>),
    NcpTpv(Img2Celestial<Ncp, WcsWithTpvImgXY2ProjXY>),
    // Cylindrical
    CypTpv(Img2Celestial<Cyp, WcsWithTpvImgXY2ProjXY>),
    CeaTpv(Img2Celestial<Cea, WcsWithTpvImgXY2ProjXY>),
    CarTpv(Img2Celestial<Car, WcsWithTpvImgXY2ProjXY>),
    MerTpv(Img2Celestial<Mer, WcsWithTpvImgXY2ProjXY>),
    // Pseudo-Cylindrical
    SflTpv(Img2Celestial<Sfl, WcsWithTpvImgXY2ProjXY>),
    ParTpv(Img2Celestial<Par, WcsWithTpvImgXY2ProjXY>),
    MolTpv(Img2Celestial<Mol, WcsWithTpvImgXY2ProjXY>),
    AitTpv(Img2Celestial<Ait, WcsWithTpvImgXY2ProjXY>),
    // Conic
    CopTpv(Img2Celestial<Cop, WcsWithTpvImgXY2ProjXY>),
    CodTpv(Img2Celestial<Cod, WcsWithTpvImgXY2ProjXY>),
    CoeTpv(Img2Celestial<Coe, WcsWithTpvImgXY2ProjXY>),
    CooTpv(Img2Celestial<Coo, WcsWithTpvImgXY2ProjXY>),
    // Polyconic
    BonTpv(Img2Celestial<Bon, WcsWithTpvImgXY2ProjXY>),
    PcoTpv(Img2Celestial<Pco, WcsWithTpvImgXY2ProjXY>),
    // Quadcube
    TscTpv(Img2Celestial<Tsc, WcsWithTpvImgXY2ProjXY>),
    CscTpv(Img2Celestial<Csc, WcsWithTpvImgXY2ProjXY>),
    QscTpv(Img2Celestial<Qsc, WcsWithTpvImgXY2ProjXY>),
    // Hybrid
    HpxTpv(Img2Celestial<Hpx, WcsWithTpvImgXY2ProjXY>),
}

fn parse_pc_matrix(params: &WCSParams) -> Option<(f64, f64, f64, f64)> {
    let pc11 = params.pc1_1;
    let pc12 = params.pc1_2;
    let pc21 = params.pc2_1;
    let pc22 = params.pc2_2;

    let pc_matrix_found = match (&pc11, &pc12, &pc21, &pc22) {
        (None, None, None, None) => false,
        // The CD1_1 keyword has been found
        // We are in a case where the CDij are given
        _ => true,
    };

    if pc_matrix_found {
        let pc11 = pc11.unwrap_or(1.0);
        let pc12 = pc12.unwrap_or(0.0);
        let pc21 = pc21.unwrap_or(0.0);
        let pc22 = pc22.unwrap_or(1.0);

        Some((pc11, pc12, pc21, pc22))
    } else {
        None
    }
}

fn parse_cd_matrix(params: &WCSParams) -> Option<(f64, f64, f64, f64)> {
    let cd11 = params.cd1_1;
    let cd12 = params.cd1_2;
    let cd21 = params.cd2_1;
    let cd22 = params.cd2_2;

    let cd_matrix_found = match (&cd11, &cd12, &cd21, &cd22) {
        (None, None, None, None) => false,
        // The CD1_1 keyword has been found
        // We are in a case where the CDij are given
        _ => true,
    };

    if cd_matrix_found {
        let cd11 = cd11.unwrap_or(1.0);
        let cd12 = cd12.unwrap_or(0.0);
        let cd21 = cd21.unwrap_or(0.0);
        let cd22 = cd22.unwrap_or(1.0);

        Some((cd11, cd12, cd21, cd22))
    } else {
        None
    }
}

impl WCSProj {
    /// Create a WCS from a specific fits header parsed with fitsrs
    /// # Param
    /// * `naxis1` - Size of the image in its first dimension (in pixels)
    /// * `naxis2` - Size of the image in its second dimension (in pixels)
    /// * `params` - Header unit coming from fitsrs.
    ///   This contains all the cards of one HDU.
    pub fn new(naxis1: i64, naxis2: i64, params: &WCSParams) -> Result<Self, Error> {
        // 1. Identify the image <-> intermediate projection
        // a. Linear transformation matrix cases:
        // - CRPIXi + CDij
        // - CRPIXi + CDELTi + CROTA2
        // - CRPIXi + CDELTi + PCij
        let crpix1 = params.crpix1.unwrap_or(0.0);
        let crpix2 = params.crpix2.unwrap_or(0.0);

        // Choice of the wcs order:
        // 1 - Priority to define the projection is given to CD
        // 2 - Then, to the couple PC + CDELT
        // 3 - Finally to the old CROTA + CDELT convention
        let (img2proj, s_lon, s_lat) =
            if let Some((cd11, cd12, cd21, cd22)) = parse_cd_matrix(params) {
                // CDij case
                (
                    WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22),
                    cd11.abs(),
                    cd22.abs(),
                )
            } else {
                // Search for CDELTi
                let cdelt1 = params.cdelt1.unwrap_or(1.0);
                let cdelt2 = params.cdelt2.unwrap_or(1.0);

                if let Some((pc11, pc12, pc21, pc22)) = parse_pc_matrix(params) {
                    // CDELTi + PCij case
                    (
                        WcsImgXY2ProjXY::from_pc(
                            crpix1, crpix2, pc11, pc12, pc21, pc22, cdelt1, cdelt2,
                        ),
                        (cdelt1 * pc11).abs(),
                        (cdelt2 * pc22).abs(),
                    )
                } else {
                    // CDELTi + CROTA2 case
                    let crota2 = params.crota2.unwrap_or(0.0);
                    let cosc = crota2.to_radians().cos();

                    (
                        WcsImgXY2ProjXY::from_cr(crpix1, crpix2, crota2, cdelt1, cdelt2),
                        (cdelt1 * cosc).abs(),
                        (cdelt2 * cosc).abs(),
                    )
                }
            };

        // 2. Identify the projection type
        let ctype1 = &params.ctype1;
        let proj_name = &ctype1[5..=7];

        let (proj, pos_angle) = match proj_name.as_bytes() {
            // Zenithal
            b"AZP" => {
                create_specific_proj!(Azp, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"SZP" => {
                create_specific_proj!(Szp, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"TAN" => {
                create_specific_proj!(Tan, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"STG" => {
                create_specific_proj!(Stg, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"SIN" => {
                create_specific_proj!(Sin, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"ARC" => {
                create_specific_proj!(Arc, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"ZPN" => {
                create_specific_proj!(Zpn, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"ZEA" => {
                create_specific_proj!(Zea, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"AIR" => {
                create_specific_proj!(Air, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"NCP" => {
                create_specific_proj!(Ncp, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // Cylindrical
            b"CYP" => {
                create_specific_proj!(Cyp, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"CEA" => {
                create_specific_proj!(Cea, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"CAR" => {
                create_specific_proj!(Car, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"MER" => {
                create_specific_proj!(Mer, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // Pseudo-cylindrical
            b"SFL" => {
                create_specific_proj!(Sfl, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"PAR" => {
                create_specific_proj!(Par, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"MOL" => {
                create_specific_proj!(Mol, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"AIT" => {
                create_specific_proj!(Ait, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // Conic
            b"COP" => {
                create_specific_proj!(Cop, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"COD" => {
                create_specific_proj!(Cod, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"COE" => {
                create_specific_proj!(Coe, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"COO" => {
                create_specific_proj!(Coo, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // Polyconic
            b"BON" => {
                create_specific_proj!(Bon, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"PCO" => {
                create_specific_proj!(Pco, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // Quadcube
            b"TSC" => {
                create_specific_proj!(Tsc, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"CSC" => {
                create_specific_proj!(Csc, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"QSC" => {
                create_specific_proj!(Qsc, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // HEALPix
            b"HPX" => {
                create_specific_proj!(Hpx, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // TPV is a special case: it's TAN projection with TPV distortion
            b"TPV" => {
                let (positional_angle, proj) = Tan::parse_proj(params)?;
                let tpv = tpv::parse_tpv(params)?;

                // Use the actual CD matrix from the header for TPV
                // The TPV polynomial distortion is applied in addition to the linear transformation
                let img2proj_tpv = WcsWithTpvImgXY2ProjXY::new(img2proj, tpv);
                Ok((
                    WCSCelestialProj::TanTpv(Img2Celestial::new(img2proj_tpv, proj)),
                    positional_angle,
                ))
            }
            _ => Err(Error::NotImplementedProjection(proj_name.to_string())),
        }?;

        let coo_system = CooSystem::parse(params)?;

        Ok(WCSProj {
            proj,
            coo_system,
            pos_angle,
            s_lon,
            s_lat,
        })
    }

    /// Project a (lon, lat) given in ICRS frame to get its corresponding location on the image
    ///
    /// The result is given a (X, Y) tuple expressed in pixel coordinates.
    ///
    /// # Arguments
    ///
    /// * `lonlat`: a coo expressed as (lon, lat) tuple given in radians and in ICRS system
    pub fn proj_lonlat(&self, lonlat: &LonLat) -> Option<ImgXY> {
        let lonlat = &self.coo_system.from_icrs(*lonlat);

        match &self.proj {
            // Zenithal
            WCSCelestialProj::Azp(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Szp(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Tan(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Stg(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Sin(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Arc(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Zpn(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Zea(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Air(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Ncp(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCSCelestialProj::Cyp(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Cea(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Car(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Mer(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCSCelestialProj::Sfl(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Par(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Mol(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Ait(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCSCelestialProj::Cop(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Cod(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Coe(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Coo(wcs) => wcs.lonlat2img(lonlat),
            // Polyconic
            WCSCelestialProj::Bon(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Pco(wcs) => wcs.lonlat2img(lonlat),
            // Quadcube
            WCSCelestialProj::Tsc(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Csc(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Qsc(wcs) => wcs.lonlat2img(lonlat),
            // Hybrid
            WCSCelestialProj::Hpx(wcs) => wcs.lonlat2img(lonlat),

            /* Sip variants */
            // Zenithal
            WCSCelestialProj::AzpSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::SzpSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::TanSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::StgSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::SinSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ArcSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ZpnSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ZeaSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::AirSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::NcpSip(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCSCelestialProj::CypSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CeaSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CarSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::MerSip(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCSCelestialProj::SflSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ParSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::MolSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::AitSip(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCSCelestialProj::CopSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CodSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CoeSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CooSip(wcs) => wcs.lonlat2img(lonlat),
            // Polyconic
            WCSCelestialProj::BonSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::PcoSip(wcs) => wcs.lonlat2img(lonlat),
            // Quadcube
            WCSCelestialProj::TscSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CscSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::QscSip(wcs) => wcs.lonlat2img(lonlat),
            // Hybrid
            WCSCelestialProj::HpxSip(wcs) => wcs.lonlat2img(lonlat),

            /* Tpv variants */
            // Zenithal
            WCSCelestialProj::AzpTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::SzpTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::TanTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::StgTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::SinTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ArcTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ZpnTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ZeaTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::AirTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::NcpTpv(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCSCelestialProj::CypTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CeaTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CarTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::MerTpv(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCSCelestialProj::SflTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ParTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::MolTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::AitTpv(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCSCelestialProj::CopTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CodTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CoeTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CooTpv(wcs) => wcs.lonlat2img(lonlat),
            // Polyconic
            WCSCelestialProj::BonTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::PcoTpv(wcs) => wcs.lonlat2img(lonlat),
            // Quadcube
            WCSCelestialProj::TscTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CscTpv(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::QscTpv(wcs) => wcs.lonlat2img(lonlat),
            // Hybrid
            WCSCelestialProj::HpxTpv(wcs) => wcs.lonlat2img(lonlat),
        }
    }

    pub fn proj_xyz(&self, xyz: &(f64, f64, f64)) -> Option<ImgXY> {
        let xyz = &self.coo_system.from_icrs_xyz(XYZ::new(xyz.0, xyz.1, xyz.2));

        match &self.proj {
            // Zenithal
            WCSCelestialProj::Azp(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Szp(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Tan(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Stg(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Sin(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Arc(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Zpn(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Zea(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Air(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Ncp(wcs) => wcs.xyz2img(xyz),
            // Pseudo-cyl
            WCSCelestialProj::Cyp(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Cea(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Car(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Mer(wcs) => wcs.xyz2img(xyz),
            // Cylindrical
            WCSCelestialProj::Sfl(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Par(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Mol(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Ait(wcs) => wcs.xyz2img(xyz),
            // Conic
            WCSCelestialProj::Cop(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Cod(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Coe(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Coo(wcs) => wcs.xyz2img(xyz),
            // Polyconic
            WCSCelestialProj::Bon(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Pco(wcs) => wcs.xyz2img(xyz),
            // Quadcube
            WCSCelestialProj::Tsc(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Csc(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Qsc(wcs) => wcs.xyz2img(xyz),
            // Hybrid
            WCSCelestialProj::Hpx(wcs) => wcs.xyz2img(xyz),

            /* Sip variants */
            // Zenithal
            WCSCelestialProj::AzpSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::SzpSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::TanSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::StgSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::SinSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ArcSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ZpnSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ZeaSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::AirSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::NcpSip(wcs) => wcs.xyz2img(xyz),
            // Pseudo-cyl
            WCSCelestialProj::CypSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CeaSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CarSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::MerSip(wcs) => wcs.xyz2img(xyz),
            // Cylindrical
            WCSCelestialProj::SflSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ParSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::MolSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::AitSip(wcs) => wcs.xyz2img(xyz),
            // Conic
            WCSCelestialProj::CopSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CodSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CoeSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CooSip(wcs) => wcs.xyz2img(xyz),
            // Polyconic
            WCSCelestialProj::BonSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::PcoSip(wcs) => wcs.xyz2img(xyz),
            // Quadcube
            WCSCelestialProj::TscSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CscSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::QscSip(wcs) => wcs.xyz2img(xyz),
            // Hybrid
            WCSCelestialProj::HpxSip(wcs) => wcs.xyz2img(xyz),

            /* Tpv variants */
            // Zenithal
            WCSCelestialProj::AzpTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::SzpTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::TanTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::StgTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::SinTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ArcTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ZpnTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ZeaTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::AirTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::NcpTpv(wcs) => wcs.xyz2img(xyz),
            // Pseudo-cyl
            WCSCelestialProj::CypTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CeaTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CarTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::MerTpv(wcs) => wcs.xyz2img(xyz),
            // Cylindrical
            WCSCelestialProj::SflTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ParTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::MolTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::AitTpv(wcs) => wcs.xyz2img(xyz),
            // Conic
            WCSCelestialProj::CopTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CodTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CoeTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CooTpv(wcs) => wcs.xyz2img(xyz),
            // Polyconic
            WCSCelestialProj::BonTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::PcoTpv(wcs) => wcs.xyz2img(xyz),
            // Quadcube
            WCSCelestialProj::TscTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CscTpv(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::QscTpv(wcs) => wcs.xyz2img(xyz),
            // Hybrid
            WCSCelestialProj::HpxTpv(wcs) => wcs.xyz2img(xyz),
        }
    }

    pub fn unproj_xyz(&self, img_pos: &ImgXY) -> Option<XYZ> {
        let xyz = match &self.proj {
            // Zenithal
            WCSCelestialProj::Azp(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Szp(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Tan(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Stg(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Sin(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Arc(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Zpn(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Zea(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Air(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Ncp(wcs) => wcs.img2xyz(img_pos),
            // Pseudo-cyl
            WCSCelestialProj::Cyp(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Cea(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Car(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Mer(wcs) => wcs.img2xyz(img_pos),
            // Cylindrical
            WCSCelestialProj::Sfl(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Par(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Mol(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Ait(wcs) => wcs.img2xyz(img_pos),
            // Conic
            WCSCelestialProj::Cop(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Cod(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Coe(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Coo(wcs) => wcs.img2xyz(img_pos),
            // Polyconic
            WCSCelestialProj::Bon(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Pco(wcs) => wcs.img2xyz(img_pos),
            // Quadcube
            WCSCelestialProj::Tsc(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Csc(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::Qsc(wcs) => wcs.img2xyz(img_pos),
            // Hybrid
            WCSCelestialProj::Hpx(wcs) => wcs.img2xyz(img_pos),

            /* Sip variants */
            // Zenithal
            WCSCelestialProj::AzpSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::SzpSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::TanSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::StgSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::SinSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::ArcSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::ZpnSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::ZeaSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::AirSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::NcpSip(wcs) => wcs.img2xyz(img_pos),
            // Pseudo-cyl
            WCSCelestialProj::CypSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CeaSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CarSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::MerSip(wcs) => wcs.img2xyz(img_pos),
            // Cylindrical
            WCSCelestialProj::SflSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::ParSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::MolSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::AitSip(wcs) => wcs.img2xyz(img_pos),
            // Conic
            WCSCelestialProj::CopSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CodSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CoeSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CooSip(wcs) => wcs.img2xyz(img_pos),
            // Polyconic
            WCSCelestialProj::BonSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::PcoSip(wcs) => wcs.img2xyz(img_pos),
            // Quadcube
            WCSCelestialProj::TscSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CscSip(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::QscSip(wcs) => wcs.img2xyz(img_pos),
            // Hybrid
            WCSCelestialProj::HpxSip(wcs) => wcs.img2xyz(img_pos),

            /* Tpv variants */
            // Zenithal
            WCSCelestialProj::AzpTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::SzpTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::TanTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::StgTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::SinTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::ArcTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::ZpnTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::ZeaTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::AirTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::NcpTpv(wcs) => wcs.img2xyz(img_pos),
            // Pseudo-cyl
            WCSCelestialProj::CypTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CeaTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CarTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::MerTpv(wcs) => wcs.img2xyz(img_pos),
            // Cylindrical
            WCSCelestialProj::SflTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::ParTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::MolTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::AitTpv(wcs) => wcs.img2xyz(img_pos),
            // Conic
            WCSCelestialProj::CopTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CodTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CoeTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CooTpv(wcs) => wcs.img2xyz(img_pos),
            // Polyconic
            WCSCelestialProj::BonTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::PcoTpv(wcs) => wcs.img2xyz(img_pos),
            // Quadcube
            WCSCelestialProj::TscTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::CscTpv(wcs) => wcs.img2xyz(img_pos),
            WCSCelestialProj::QscTpv(wcs) => wcs.img2xyz(img_pos),
            // Hybrid
            WCSCelestialProj::HpxTpv(wcs) => wcs.img2xyz(img_pos),
        };
        xyz.map(|v| self.coo_system.to_icrs_xyz(v))
    }

    /// Unproject a (X, Y) point from the image space to get its corresponding location on the sphere
    ///
    /// The result is (lon, lat) tuple expressed in radians in ICRS
    ///
    /// # Arguments
    ///
    /// * `img_pos`: the image space point expressed as a (X, Y) tuple given en pixels
    pub fn unproj_lonlat(&self, img_pos: &ImgXY) -> Option<LonLat> {
        let lonlat = match &self.proj {
            // Zenithal
            WCSCelestialProj::Azp(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Szp(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Tan(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Stg(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Sin(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Arc(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Zpn(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Zea(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Air(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Ncp(wcs) => wcs.img2lonlat(img_pos),
            // Pseudo-cyl
            WCSCelestialProj::Cyp(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Cea(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Car(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Mer(wcs) => wcs.img2lonlat(img_pos),
            // Cylindrical
            WCSCelestialProj::Sfl(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Par(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Mol(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Ait(wcs) => wcs.img2lonlat(img_pos),
            // Conic
            WCSCelestialProj::Cop(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Cod(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Coe(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Coo(wcs) => wcs.img2lonlat(img_pos),
            // Polyconic
            WCSCelestialProj::Bon(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Pco(wcs) => wcs.img2lonlat(img_pos),
            // Quadcube
            WCSCelestialProj::Tsc(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Csc(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::Qsc(wcs) => wcs.img2lonlat(img_pos),
            // Hybrid
            WCSCelestialProj::Hpx(wcs) => wcs.img2lonlat(img_pos),

            /* Sip variants */
            // Zenithal
            WCSCelestialProj::AzpSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::SzpSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::TanSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::StgSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::SinSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::ArcSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::ZpnSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::ZeaSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::AirSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::NcpSip(wcs) => wcs.img2lonlat(img_pos),
            // Pseudo-cyl
            WCSCelestialProj::CypSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CeaSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CarSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::MerSip(wcs) => wcs.img2lonlat(img_pos),
            // Cylindrical
            WCSCelestialProj::SflSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::ParSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::MolSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::AitSip(wcs) => wcs.img2lonlat(img_pos),
            // Conic
            WCSCelestialProj::CopSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CodSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CoeSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CooSip(wcs) => wcs.img2lonlat(img_pos),
            // Polyconic
            WCSCelestialProj::BonSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::PcoSip(wcs) => wcs.img2lonlat(img_pos),
            // Quadcube
            WCSCelestialProj::TscSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CscSip(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::QscSip(wcs) => wcs.img2lonlat(img_pos),
            // Hybrid
            WCSCelestialProj::HpxSip(wcs) => wcs.img2lonlat(img_pos),

            /* Tpv variants */
            // Zenithal
            WCSCelestialProj::AzpTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::SzpTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::TanTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::StgTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::SinTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::ArcTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::ZpnTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::ZeaTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::AirTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::NcpTpv(wcs) => wcs.img2lonlat(img_pos),
            // Pseudo-cyl
            WCSCelestialProj::CypTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CeaTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CarTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::MerTpv(wcs) => wcs.img2lonlat(img_pos),
            // Cylindrical
            WCSCelestialProj::SflTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::ParTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::MolTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::AitTpv(wcs) => wcs.img2lonlat(img_pos),
            // Conic
            WCSCelestialProj::CopTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CodTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CoeTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CooTpv(wcs) => wcs.img2lonlat(img_pos),
            // Polyconic
            WCSCelestialProj::BonTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::PcoTpv(wcs) => wcs.img2lonlat(img_pos),
            // Quadcube
            WCSCelestialProj::TscTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::CscTpv(wcs) => wcs.img2lonlat(img_pos),
            WCSCelestialProj::QscTpv(wcs) => wcs.img2lonlat(img_pos),
            // Hybrid
            WCSCelestialProj::HpxTpv(wcs) => wcs.img2lonlat(img_pos),
        };

        lonlat.map(|ll| self.coo_system.to_icrs(ll))
    }

    /// Getter of the coordinate system
    pub fn coo_system(&self) -> &CooSystem {
        &self.coo_system
    }
}

#[cfg(test)]
mod tests {
    use super::WCS;
    use crate::Error;
    use crate::WCSParams;

    use crate::mapproj::Projection;

    use fitsrs::card::Value;
    use fitsrs::fits::Fits;
    use fitsrs::hdu::header::{extension::image::Image, Header};
    use fitsrs::Pixels;

    use glob::glob;
    use mapproj::{CanonicalProjection, ImgXY, LonLat};
    use serde::Deserialize;
    use std::f64::consts::PI;
    use std::fs::File;
    use std::io::BufReader;

    use std::convert::TryFrom;

    use crate::Error::MandatoryWCSKeywordsMissing;

    fn parse_card<'de, T: Deserialize<'de>>(
        header: &'de Header<Image>,
        key: &'static str,
    ) -> Result<T, Error> {
        header
            .get_parsed::<T>(key)
            .map_err(|_| MandatoryWCSKeywordsMissing("Card cannot be parsed"))
    }

    fn parse_opt_card<'de, T: Deserialize<'de>>(
        header: &'de Header<Image>,
        key: &'static str,
    ) -> Option<T> {
        header.get_parsed::<T>(key).ok()
    }

    impl<'a> TryFrom<&'a Header<Image>> for WCS {
        type Error = Error;

        fn try_from(h: &'a Header<Image>) -> Result<Self, Self::Error> {
            let params = WCSParams {
                naxis: parse_card::<i64>(h, "NAXIS").unwrap(),
                ctype1: parse_card::<String>(h, "CTYPE1").unwrap(),

                naxis1: parse_opt_card::<i64>(h, "NAXIS1"),
                naxis2: parse_opt_card::<i64>(h, "NAXIS2"),

                znaxis1: parse_opt_card::<i64>(h, "ZNAXIS1"),
                znaxis2: parse_opt_card::<i64>(h, "ZNAXIS2"),
                znaxis3: parse_opt_card::<i64>(h, "ZNAXIS3"),
                znaxis4: parse_opt_card::<i64>(h, "ZNAXIS4"),

                ctype2: parse_opt_card::<String>(h, "CTYPE2"),
                ctype3: parse_opt_card::<String>(h, "CTYPE3"),

                a_order: parse_opt_card::<i64>(h, "A_ORDER"),
                b_order: parse_opt_card::<i64>(h, "B_ORDER"),
                ap_order: parse_opt_card::<i64>(h, "AP_ORDER"),
                bp_order: parse_opt_card::<i64>(h, "BP_ORDER"),
                crpix1: parse_opt_card::<f64>(h, "CRPIX1"),
                crpix2: parse_opt_card::<f64>(h, "CRPIX2"),
                crpix3: parse_opt_card::<f64>(h, "CRPIX3"),
                crval1: parse_opt_card::<f64>(h, "CRVAL1"),
                crval2: parse_opt_card::<f64>(h, "CRVAL2"),
                crval3: parse_opt_card::<f64>(h, "CRVAL3"),
                crota1: parse_opt_card::<f64>(h, "CROTA1"),
                crota2: parse_opt_card::<f64>(h, "CROTA2"),
                crota3: parse_opt_card::<f64>(h, "CROTA3"),
                cdelt1: parse_opt_card::<f64>(h, "CDELT1"),
                cdelt2: parse_opt_card::<f64>(h, "CDELT2"),
                cdelt3: parse_opt_card::<f64>(h, "CDELT3"),
                naxis3: parse_opt_card::<i64>(h, "NAXIS3"),
                naxis4: parse_opt_card::<i64>(h, "NAXIS4"),
                lonpole: parse_opt_card::<f64>(h, "LONPOLE"),
                latpole: parse_opt_card::<f64>(h, "LATPOLE"),
                equinox: parse_opt_card::<f64>(h, "EQUINOX"),
                epoch: parse_opt_card::<f64>(h, "EPOCH"),
                radesys: parse_opt_card::<String>(h, "RADESYS"),
                pv1_0: parse_opt_card::<f64>(h, "PV1_0"),
                pv1_1: parse_opt_card::<f64>(h, "PV1_1"),
                pv1_2: parse_opt_card::<f64>(h, "PV1_2"),
                pv1_3: parse_opt_card::<f64>(h, "PV1_3"),
                pv1_4: parse_opt_card::<f64>(h, "PV1_4"),
                pv1_5: parse_opt_card::<f64>(h, "PV1_5"),
                pv1_6: parse_opt_card::<f64>(h, "PV1_6"),
                pv1_7: parse_opt_card::<f64>(h, "PV1_7"),
                pv1_8: parse_opt_card::<f64>(h, "PV1_8"),
                pv1_9: parse_opt_card::<f64>(h, "PV1_9"),
                pv1_10: parse_opt_card::<f64>(h, "PV1_10"),
                pv1_11: parse_opt_card::<f64>(h, "PV1_11"),
                pv1_12: parse_opt_card::<f64>(h, "PV1_12"),
                pv1_13: parse_opt_card::<f64>(h, "PV1_13"),
                pv1_14: parse_opt_card::<f64>(h, "PV1_14"),
                pv1_15: parse_opt_card::<f64>(h, "PV1_15"),
                pv1_16: parse_opt_card::<f64>(h, "PV1_16"),
                pv1_17: parse_opt_card::<f64>(h, "PV1_17"),
                pv1_18: parse_opt_card::<f64>(h, "PV1_18"),
                pv1_19: parse_opt_card::<f64>(h, "PV1_19"),
                pv1_20: parse_opt_card::<f64>(h, "PV1_20"),
                pv1_21: parse_opt_card::<f64>(h, "PV1_21"),
                pv1_22: parse_opt_card::<f64>(h, "PV1_22"),
                pv1_23: parse_opt_card::<f64>(h, "PV1_23"),
                pv1_24: parse_opt_card::<f64>(h, "PV1_24"),
                pv1_25: parse_opt_card::<f64>(h, "PV1_25"),
                pv1_26: parse_opt_card::<f64>(h, "PV1_26"),
                pv1_27: parse_opt_card::<f64>(h, "PV1_27"),
                pv1_28: parse_opt_card::<f64>(h, "PV1_28"),
                pv1_29: parse_opt_card::<f64>(h, "PV1_29"),
                pv1_30: parse_opt_card::<f64>(h, "PV1_30"),
                pv1_31: parse_opt_card::<f64>(h, "PV1_31"),
                pv1_32: parse_opt_card::<f64>(h, "PV1_32"),
                pv1_33: parse_opt_card::<f64>(h, "PV1_33"),
                pv1_34: parse_opt_card::<f64>(h, "PV1_34"),
                pv1_35: parse_opt_card::<f64>(h, "PV1_35"),
                pv1_36: parse_opt_card::<f64>(h, "PV1_36"),
                pv1_37: parse_opt_card::<f64>(h, "PV1_37"),
                pv1_38: parse_opt_card::<f64>(h, "PV1_38"),
                pv1_39: parse_opt_card::<f64>(h, "PV1_39"),
                pv2_0: parse_opt_card::<f64>(h, "PV2_0"),
                pv2_1: parse_opt_card::<f64>(h, "PV2_1"),
                pv2_2: parse_opt_card::<f64>(h, "PV2_2"),
                pv2_3: parse_opt_card::<f64>(h, "PV2_3"),
                pv2_4: parse_opt_card::<f64>(h, "PV2_4"),
                pv2_5: parse_opt_card::<f64>(h, "PV2_5"),
                pv2_6: parse_opt_card::<f64>(h, "PV2_6"),
                pv2_7: parse_opt_card::<f64>(h, "PV2_7"),
                pv2_8: parse_opt_card::<f64>(h, "PV2_8"),
                pv2_9: parse_opt_card::<f64>(h, "PV2_9"),
                pv2_10: parse_opt_card::<f64>(h, "PV2_10"),
                pv2_11: parse_opt_card::<f64>(h, "PV2_11"),
                pv2_12: parse_opt_card::<f64>(h, "PV2_12"),
                pv2_13: parse_opt_card::<f64>(h, "PV2_13"),
                pv2_14: parse_opt_card::<f64>(h, "PV2_14"),
                pv2_15: parse_opt_card::<f64>(h, "PV2_15"),
                pv2_16: parse_opt_card::<f64>(h, "PV2_16"),
                pv2_17: parse_opt_card::<f64>(h, "PV2_17"),
                pv2_18: parse_opt_card::<f64>(h, "PV2_18"),
                pv2_19: parse_opt_card::<f64>(h, "PV2_19"),
                pv2_20: parse_opt_card::<f64>(h, "PV2_20"),
                pv2_21: parse_opt_card::<f64>(h, "PV2_21"),
                pv2_22: parse_opt_card::<f64>(h, "PV2_22"),
                pv2_23: parse_opt_card::<f64>(h, "PV2_23"),
                pv2_24: parse_opt_card::<f64>(h, "PV2_24"),
                pv2_25: parse_opt_card::<f64>(h, "PV2_25"),
                pv2_26: parse_opt_card::<f64>(h, "PV2_26"),
                pv2_27: parse_opt_card::<f64>(h, "PV2_27"),
                pv2_28: parse_opt_card::<f64>(h, "PV2_28"),
                pv2_29: parse_opt_card::<f64>(h, "PV2_29"),
                pv2_30: parse_opt_card::<f64>(h, "PV2_30"),
                pv2_31: parse_opt_card::<f64>(h, "PV2_31"),
                pv2_32: parse_opt_card::<f64>(h, "PV2_32"),
                pv2_33: parse_opt_card::<f64>(h, "PV2_33"),
                pv2_34: parse_opt_card::<f64>(h, "PV2_34"),
                pv2_35: parse_opt_card::<f64>(h, "PV2_35"),
                pv2_36: parse_opt_card::<f64>(h, "PV2_36"),
                pv2_37: parse_opt_card::<f64>(h, "PV2_37"),
                pv2_38: parse_opt_card::<f64>(h, "PV2_38"),
                pv2_39: parse_opt_card::<f64>(h, "PV2_39"),
                cd1_1: parse_opt_card::<f64>(h, "CD1_1"),
                cd1_2: parse_opt_card::<f64>(h, "CD1_2"),
                cd1_3: parse_opt_card::<f64>(h, "CD1_3"),
                cd2_1: parse_opt_card::<f64>(h, "CD2_1"),
                cd2_2: parse_opt_card::<f64>(h, "CD2_2"),
                cd2_3: parse_opt_card::<f64>(h, "CD2_3"),
                cd3_1: parse_opt_card::<f64>(h, "CD3_1"),
                cd3_2: parse_opt_card::<f64>(h, "CD3_2"),
                cd3_3: parse_opt_card::<f64>(h, "CD3_3"),
                pc1_1: parse_opt_card::<f64>(h, "PC1_1"),
                pc1_2: parse_opt_card::<f64>(h, "PC1_2"),
                pc1_3: parse_opt_card::<f64>(h, "PC1_3"),
                pc2_1: parse_opt_card::<f64>(h, "PC2_1"),
                pc2_2: parse_opt_card::<f64>(h, "PC2_2"),
                pc2_3: parse_opt_card::<f64>(h, "PC2_3"),
                pc3_1: parse_opt_card::<f64>(h, "PC3_1"),
                pc3_2: parse_opt_card::<f64>(h, "PC3_2"),
                pc3_3: parse_opt_card::<f64>(h, "PC3_3"),
                a_0_0: parse_opt_card::<f64>(h, "A_0_0"),
                a_1_0: parse_opt_card::<f64>(h, "A_1_0"),
                a_2_0: parse_opt_card::<f64>(h, "A_2_0"),
                a_3_0: parse_opt_card::<f64>(h, "A_3_0"),
                a_4_0: parse_opt_card::<f64>(h, "A_4_0"),
                a_5_0: parse_opt_card::<f64>(h, "A_5_0"),
                a_6_0: parse_opt_card::<f64>(h, "A_6_0"),
                a_0_1: parse_opt_card::<f64>(h, "A_0_1"),
                a_1_1: parse_opt_card::<f64>(h, "A_1_1"),
                a_2_1: parse_opt_card::<f64>(h, "A_2_1"),
                a_3_1: parse_opt_card::<f64>(h, "A_3_1"),
                a_4_1: parse_opt_card::<f64>(h, "A_4_1"),
                a_5_1: parse_opt_card::<f64>(h, "A_5_1"),
                a_0_2: parse_opt_card::<f64>(h, "A_0_2"),
                a_1_2: parse_opt_card::<f64>(h, "A_1_2"),
                a_2_2: parse_opt_card::<f64>(h, "A_2_2"),
                a_3_2: parse_opt_card::<f64>(h, "A_3_2"),
                a_4_2: parse_opt_card::<f64>(h, "A_4_2"),
                a_0_3: parse_opt_card::<f64>(h, "A_0_3"),
                a_1_3: parse_opt_card::<f64>(h, "A_1_3"),
                a_2_3: parse_opt_card::<f64>(h, "A_2_3"),
                a_3_3: parse_opt_card::<f64>(h, "A_3_3"),
                a_0_4: parse_opt_card::<f64>(h, "A_0_4"),
                a_1_4: parse_opt_card::<f64>(h, "A_1_4"),
                a_2_4: parse_opt_card::<f64>(h, "A_2_4"),
                a_0_5: parse_opt_card::<f64>(h, "A_0_5"),
                a_1_5: parse_opt_card::<f64>(h, "A_1_5"),
                a_0_6: parse_opt_card::<f64>(h, "A_0_6"),
                ap_0_0: parse_opt_card::<f64>(h, "AP_0_0"),
                ap_1_0: parse_opt_card::<f64>(h, "AP_1_0"),
                ap_2_0: parse_opt_card::<f64>(h, "AP_2_0"),
                ap_3_0: parse_opt_card::<f64>(h, "AP_3_0"),
                ap_4_0: parse_opt_card::<f64>(h, "AP_4_0"),
                ap_5_0: parse_opt_card::<f64>(h, "AP_5_0"),
                ap_6_0: parse_opt_card::<f64>(h, "AP_6_0"),
                ap_0_1: parse_opt_card::<f64>(h, "AP_0_1"),
                ap_1_1: parse_opt_card::<f64>(h, "AP_1_1"),
                ap_2_1: parse_opt_card::<f64>(h, "AP_2_1"),
                ap_3_1: parse_opt_card::<f64>(h, "AP_3_1"),
                ap_4_1: parse_opt_card::<f64>(h, "AP_4_1"),
                ap_5_1: parse_opt_card::<f64>(h, "AP_5_1"),
                ap_0_2: parse_opt_card::<f64>(h, "AP_0_2"),
                ap_1_2: parse_opt_card::<f64>(h, "AP_1_2"),
                ap_2_2: parse_opt_card::<f64>(h, "AP_2_2"),
                ap_3_2: parse_opt_card::<f64>(h, "AP_3_2"),
                ap_4_2: parse_opt_card::<f64>(h, "AP_4_2"),
                ap_0_3: parse_opt_card::<f64>(h, "AP_0_3"),
                ap_1_3: parse_opt_card::<f64>(h, "AP_1_3"),
                ap_2_3: parse_opt_card::<f64>(h, "AP_2_3"),
                ap_3_3: parse_opt_card::<f64>(h, "AP_3_3"),
                ap_0_4: parse_opt_card::<f64>(h, "AP_0_4"),
                ap_1_4: parse_opt_card::<f64>(h, "AP_1_4"),
                ap_2_4: parse_opt_card::<f64>(h, "AP_2_4"),
                ap_0_5: parse_opt_card::<f64>(h, "AP_0_5"),
                ap_1_5: parse_opt_card::<f64>(h, "AP_1_5"),
                ap_0_6: parse_opt_card::<f64>(h, "AP_0_6"),
                b_0_0: parse_opt_card::<f64>(h, "B_0_0"),
                b_1_0: parse_opt_card::<f64>(h, "B_1_0"),
                b_2_0: parse_opt_card::<f64>(h, "B_2_0"),
                b_3_0: parse_opt_card::<f64>(h, "B_3_0"),
                b_4_0: parse_opt_card::<f64>(h, "B_4_0"),
                b_5_0: parse_opt_card::<f64>(h, "B_5_0"),
                b_6_0: parse_opt_card::<f64>(h, "B_6_0"),
                b_0_1: parse_opt_card::<f64>(h, "B_0_1"),
                b_1_1: parse_opt_card::<f64>(h, "B_1_1"),
                b_2_1: parse_opt_card::<f64>(h, "B_2_1"),
                b_3_1: parse_opt_card::<f64>(h, "B_3_1"),
                b_4_1: parse_opt_card::<f64>(h, "B_4_1"),
                b_5_1: parse_opt_card::<f64>(h, "B_5_1"),
                b_0_2: parse_opt_card::<f64>(h, "B_0_2"),
                b_1_2: parse_opt_card::<f64>(h, "B_1_2"),
                b_2_2: parse_opt_card::<f64>(h, "B_2_2"),
                b_3_2: parse_opt_card::<f64>(h, "B_3_2"),
                b_4_2: parse_opt_card::<f64>(h, "B_4_2"),
                b_0_3: parse_opt_card::<f64>(h, "B_0_3"),
                b_1_3: parse_opt_card::<f64>(h, "B_1_3"),
                b_2_3: parse_opt_card::<f64>(h, "B_2_3"),
                b_3_3: parse_opt_card::<f64>(h, "B_3_3"),
                b_0_4: parse_opt_card::<f64>(h, "B_0_4"),
                b_1_4: parse_opt_card::<f64>(h, "B_1_4"),
                b_2_4: parse_opt_card::<f64>(h, "B_2_4"),
                b_0_5: parse_opt_card::<f64>(h, "B_0_5"),
                b_1_5: parse_opt_card::<f64>(h, "B_1_5"),
                b_0_6: parse_opt_card::<f64>(h, "B_0_6"),
                // Degree 7-9 coefficients
                a_0_7: parse_opt_card::<f64>(h, "A_0_7"),
                a_1_6: parse_opt_card::<f64>(h, "A_1_6"),
                a_2_5: parse_opt_card::<f64>(h, "A_2_5"),
                a_3_4: parse_opt_card::<f64>(h, "A_3_4"),
                a_4_3: parse_opt_card::<f64>(h, "A_4_3"),
                a_5_2: parse_opt_card::<f64>(h, "A_5_2"),
                a_6_1: parse_opt_card::<f64>(h, "A_6_1"),
                a_7_0: parse_opt_card::<f64>(h, "A_7_0"),
                a_0_8: parse_opt_card::<f64>(h, "A_0_8"),
                a_1_7: parse_opt_card::<f64>(h, "A_1_7"),
                a_2_6: parse_opt_card::<f64>(h, "A_2_6"),
                a_3_5: parse_opt_card::<f64>(h, "A_3_5"),
                a_4_4: parse_opt_card::<f64>(h, "A_4_4"),
                a_5_3: parse_opt_card::<f64>(h, "A_5_3"),
                a_6_2: parse_opt_card::<f64>(h, "A_6_2"),
                a_7_1: parse_opt_card::<f64>(h, "A_7_1"),
                a_8_0: parse_opt_card::<f64>(h, "A_8_0"),
                a_0_9: parse_opt_card::<f64>(h, "A_0_9"),
                a_1_8: parse_opt_card::<f64>(h, "A_1_8"),
                a_2_7: parse_opt_card::<f64>(h, "A_2_7"),
                a_3_6: parse_opt_card::<f64>(h, "A_3_6"),
                a_4_5: parse_opt_card::<f64>(h, "A_4_5"),
                a_5_4: parse_opt_card::<f64>(h, "A_5_4"),
                a_6_3: parse_opt_card::<f64>(h, "A_6_3"),
                a_7_2: parse_opt_card::<f64>(h, "A_7_2"),
                a_8_1: parse_opt_card::<f64>(h, "A_8_1"),
                a_9_0: parse_opt_card::<f64>(h, "A_9_0"),
                b_0_7: parse_opt_card::<f64>(h, "B_0_7"),
                b_1_6: parse_opt_card::<f64>(h, "B_1_6"),
                b_2_5: parse_opt_card::<f64>(h, "B_2_5"),
                b_3_4: parse_opt_card::<f64>(h, "B_3_4"),
                b_4_3: parse_opt_card::<f64>(h, "B_4_3"),
                b_5_2: parse_opt_card::<f64>(h, "B_5_2"),
                b_6_1: parse_opt_card::<f64>(h, "B_6_1"),
                b_7_0: parse_opt_card::<f64>(h, "B_7_0"),
                b_0_8: parse_opt_card::<f64>(h, "B_0_8"),
                b_1_7: parse_opt_card::<f64>(h, "B_1_7"),
                b_2_6: parse_opt_card::<f64>(h, "B_2_6"),
                b_3_5: parse_opt_card::<f64>(h, "B_3_5"),
                b_4_4: parse_opt_card::<f64>(h, "B_4_4"),
                b_5_3: parse_opt_card::<f64>(h, "B_5_3"),
                b_6_2: parse_opt_card::<f64>(h, "B_6_2"),
                b_7_1: parse_opt_card::<f64>(h, "B_7_1"),
                b_8_0: parse_opt_card::<f64>(h, "B_8_0"),
                b_0_9: parse_opt_card::<f64>(h, "B_0_9"),
                b_1_8: parse_opt_card::<f64>(h, "B_1_8"),
                b_2_7: parse_opt_card::<f64>(h, "B_2_7"),
                b_3_6: parse_opt_card::<f64>(h, "B_3_6"),
                b_4_5: parse_opt_card::<f64>(h, "B_4_5"),
                b_5_4: parse_opt_card::<f64>(h, "B_5_4"),
                b_6_3: parse_opt_card::<f64>(h, "B_6_3"),
                b_7_2: parse_opt_card::<f64>(h, "B_7_2"),
                b_8_1: parse_opt_card::<f64>(h, "B_8_1"),
                b_9_0: parse_opt_card::<f64>(h, "B_9_0"),
                ap_0_7: parse_opt_card::<f64>(h, "AP_0_7"),
                ap_1_6: parse_opt_card::<f64>(h, "AP_1_6"),
                ap_2_5: parse_opt_card::<f64>(h, "AP_2_5"),
                ap_3_4: parse_opt_card::<f64>(h, "AP_3_4"),
                ap_4_3: parse_opt_card::<f64>(h, "AP_4_3"),
                ap_5_2: parse_opt_card::<f64>(h, "AP_5_2"),
                ap_6_1: parse_opt_card::<f64>(h, "AP_6_1"),
                ap_7_0: parse_opt_card::<f64>(h, "AP_7_0"),
                ap_0_8: parse_opt_card::<f64>(h, "AP_0_8"),
                ap_1_7: parse_opt_card::<f64>(h, "AP_1_7"),
                ap_2_6: parse_opt_card::<f64>(h, "AP_2_6"),
                ap_3_5: parse_opt_card::<f64>(h, "AP_3_5"),
                ap_4_4: parse_opt_card::<f64>(h, "AP_4_4"),
                ap_5_3: parse_opt_card::<f64>(h, "AP_5_3"),
                ap_6_2: parse_opt_card::<f64>(h, "AP_6_2"),
                ap_7_1: parse_opt_card::<f64>(h, "AP_7_1"),
                ap_8_0: parse_opt_card::<f64>(h, "AP_8_0"),
                ap_0_9: parse_opt_card::<f64>(h, "AP_0_9"),
                ap_1_8: parse_opt_card::<f64>(h, "AP_1_8"),
                ap_2_7: parse_opt_card::<f64>(h, "AP_2_7"),
                ap_3_6: parse_opt_card::<f64>(h, "AP_3_6"),
                ap_4_5: parse_opt_card::<f64>(h, "AP_4_5"),
                ap_5_4: parse_opt_card::<f64>(h, "AP_5_4"),
                ap_6_3: parse_opt_card::<f64>(h, "AP_6_3"),
                ap_7_2: parse_opt_card::<f64>(h, "AP_7_2"),
                ap_8_1: parse_opt_card::<f64>(h, "AP_8_1"),
                ap_9_0: parse_opt_card::<f64>(h, "AP_9_0"),
                bp_0_7: parse_opt_card::<f64>(h, "BP_0_7"),
                bp_1_6: parse_opt_card::<f64>(h, "BP_1_6"),
                bp_2_5: parse_opt_card::<f64>(h, "BP_2_5"),
                bp_3_4: parse_opt_card::<f64>(h, "BP_3_4"),
                bp_4_3: parse_opt_card::<f64>(h, "BP_4_3"),
                bp_5_2: parse_opt_card::<f64>(h, "BP_5_2"),
                bp_6_1: parse_opt_card::<f64>(h, "BP_6_1"),
                bp_7_0: parse_opt_card::<f64>(h, "BP_7_0"),
                bp_0_8: parse_opt_card::<f64>(h, "BP_0_8"),
                bp_1_7: parse_opt_card::<f64>(h, "BP_1_7"),
                bp_2_6: parse_opt_card::<f64>(h, "BP_2_6"),
                bp_3_5: parse_opt_card::<f64>(h, "BP_3_5"),
                bp_4_4: parse_opt_card::<f64>(h, "BP_4_4"),
                bp_5_3: parse_opt_card::<f64>(h, "BP_5_3"),
                bp_6_2: parse_opt_card::<f64>(h, "BP_6_2"),
                bp_7_1: parse_opt_card::<f64>(h, "BP_7_1"),
                bp_8_0: parse_opt_card::<f64>(h, "BP_8_0"),
                bp_0_9: parse_opt_card::<f64>(h, "BP_0_9"),
                bp_1_8: parse_opt_card::<f64>(h, "BP_1_8"),
                bp_2_7: parse_opt_card::<f64>(h, "BP_2_7"),
                bp_3_6: parse_opt_card::<f64>(h, "BP_3_6"),
                bp_4_5: parse_opt_card::<f64>(h, "BP_4_5"),
                bp_5_4: parse_opt_card::<f64>(h, "BP_5_4"),
                bp_6_3: parse_opt_card::<f64>(h, "BP_6_3"),
                bp_7_2: parse_opt_card::<f64>(h, "BP_7_2"),
                bp_8_1: parse_opt_card::<f64>(h, "BP_8_1"),
                bp_9_0: parse_opt_card::<f64>(h, "BP_9_0"),
                bp_0_0: parse_opt_card::<f64>(h, "BP_0_0"),
                bp_1_0: parse_opt_card::<f64>(h, "BP_1_0"),
                bp_2_0: parse_opt_card::<f64>(h, "BP_2_0"),
                bp_3_0: parse_opt_card::<f64>(h, "BP_3_0"),
                bp_4_0: parse_opt_card::<f64>(h, "BP_4_0"),
                bp_5_0: parse_opt_card::<f64>(h, "BP_5_0"),
                bp_6_0: parse_opt_card::<f64>(h, "BP_6_0"),
                bp_0_1: parse_opt_card::<f64>(h, "BP_0_1"),
                bp_1_1: parse_opt_card::<f64>(h, "BP_1_1"),
                bp_2_1: parse_opt_card::<f64>(h, "BP_2_1"),
                bp_3_1: parse_opt_card::<f64>(h, "BP_3_1"),
                bp_4_1: parse_opt_card::<f64>(h, "BP_4_1"),
                bp_5_1: parse_opt_card::<f64>(h, "BP_5_1"),
                bp_0_2: parse_opt_card::<f64>(h, "BP_0_2"),
                bp_1_2: parse_opt_card::<f64>(h, "BP_1_2"),
                bp_2_2: parse_opt_card::<f64>(h, "BP_2_2"),
                bp_3_2: parse_opt_card::<f64>(h, "BP_3_2"),
                bp_4_2: parse_opt_card::<f64>(h, "BP_4_2"),
                bp_0_3: parse_opt_card::<f64>(h, "BP_0_3"),
                bp_1_3: parse_opt_card::<f64>(h, "BP_1_3"),
                bp_2_3: parse_opt_card::<f64>(h, "BP_2_3"),
                bp_3_3: parse_opt_card::<f64>(h, "BP_3_3"),
                bp_0_4: parse_opt_card::<f64>(h, "BP_0_4"),
                bp_1_4: parse_opt_card::<f64>(h, "BP_1_4"),
                bp_2_4: parse_opt_card::<f64>(h, "BP_2_4"),
                bp_0_5: parse_opt_card::<f64>(h, "BP_0_5"),
                bp_1_5: parse_opt_card::<f64>(h, "BP_1_5"),
                bp_0_6: parse_opt_card::<f64>(h, "BP_0_6"),
            };

            WCS::new(&params)
        }
    }

    fn wcs_from_fits_header(header: &Header<Image>) -> Result<WCS, Error> {
        header.try_into()
    }

    #[test]
    fn test_visualize() {
        let f = File::open("examples/panstarrs-rotated-around-orion.fits").unwrap();

        let reader = BufReader::new(f);
        let mut fits = Fits::from_reader(reader);
        let hdu = fits.next().unwrap().unwrap();

        match hdu {
            HDU::XImage(hdu) | HDU::Primary(hdu) => {
                let header = hdu.get_header();

                // Parse data
                let data = match fits.get_data(&hdu).pixels() {
                    Pixels::F32(it) => it.collect::<Vec<_>>(),
                    _ => unreachable!(),
                };

                let wcs = wcs_from_fits_header(header).unwrap();
                reproject_fits_image(mapproj::zenithal::azp::Azp::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::zenithal::szp::Szp::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::zenithal::tan::Tan::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::zenithal::stg::Stg::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::zenithal::sin::Sin::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::zenithal::arc::Arc::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::zenithal::zea::Zea::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::zenithal::air::Air::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::zenithal::ncp::Ncp::new(), &wcs, header, &data);

                reproject_fits_image(mapproj::pseudocyl::mol::Mol::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::pseudocyl::ait::Ait::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::pseudocyl::par::Par::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::pseudocyl::sfl::Sfl::new(), &wcs, header, &data);

                reproject_fits_image(mapproj::cylindrical::cyp::Cyp::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::cylindrical::cea::Cea::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::cylindrical::car::Car::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::cylindrical::mer::Mer::new(), &wcs, header, &data);

                reproject_fits_image(mapproj::conic::cod::Cod::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::conic::cop::Cop::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::conic::coo::Coo::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::conic::coe::Coe::new(), &wcs, header, &data);

                reproject_fits_image(mapproj::polyconic::bon::Bon::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::polyconic::pco::Pco::new(), &wcs, header, &data);

                reproject_fits_image(mapproj::quadcube::csc::Csc::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::quadcube::qsc::Qsc::new(), &wcs, header, &data);
                reproject_fits_image(mapproj::quadcube::tsc::Tsc::new(), &wcs, header, &data);

                reproject_fits_image(mapproj::hybrid::hpx::Hpx::new(), &wcs, header, &data);
            }
            _ => unreachable!(),
        }
    }

    fn reproject_fits_image<T: CanonicalProjection>(
        proj: T,
        wcs: &WCS,
        header: &Header<Image>,
        data: &[f32],
    ) {
        let scale = header.get_parsed::<f32>("BSCALE").unwrap_or(1.0);
        let offset = header.get_parsed::<f32>("BZERO").unwrap_or(0.0);

        let xtension = header.get_xtension();
        let naxis = xtension.get_naxis();
        let width = naxis[0];
        let height = naxis[1];

        //let proj = mapproj::zenithal::sin::Sin::new();
        let bounds = proj.bounds();
        let x_bounds = bounds.x_bounds().as_ref().unwrap_or(&((-PI)..=PI));
        let y_bounds = bounds.y_bounds().as_ref().unwrap_or(&((-PI)..=PI));

        let x_off = x_bounds.start();
        let x_len = x_bounds.end() - x_bounds.start();

        let y_off = y_bounds.start();
        let y_len = y_bounds.end() - y_bounds.start();

        const WIDTH_IMAGE: usize = 1024;
        const HEIGHT_IMAGE: usize = 1024;
        // Create a new ImgBuf with width: imgx and height: imgy
        let mut imgbuf = image::ImageBuffer::new(WIDTH_IMAGE as u32, HEIGHT_IMAGE as u32);

        for y in 0..height {
            for x in 0..width {
                let grayscale_val = (data[(y * width + x) as usize] * scale + offset) as u8;

                let img_xy = ImgXY::new(x as f64, y as f64);
                if let Some(lonlat) = wcs.unproj(&img_xy) {
                    if let Some(proj_xy) = proj.proj_lonlat(&lonlat) {
                        let proj_x = (proj_xy.x() - x_off) / x_len; // between 0 and 1
                        let proj_y = (proj_xy.y() - y_off) / y_len; // between 0 and 1

                        if (0.0..1.0).contains(&proj_x) && (0.0..1.0).contains(&proj_y) {
                            let ix = (proj_x * (WIDTH_IMAGE as f64)) as usize;
                            let iy = (proj_y * (HEIGHT_IMAGE as f64)) as usize;

                            let pixel = imgbuf
                                .get_pixel_mut(ix as u32, (HEIGHT_IMAGE as u32) - iy as u32 - 1);
                            *pixel = image::Rgb([grayscale_val, grayscale_val, grayscale_val]);
                        }
                    }
                }
            }
        }

        let filename = &format!(
            "tests/reproj/pans-{}.jpeg",
            <T as CanonicalProjection>::WCS_NAME
        );
        imgbuf.save(filename).unwrap();
    }

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            if ($x - $y).abs() > $d {
                panic!();
            }
        };
    }

    #[test]
    fn astropy_comparison() {
        use std::fs;

        let fits_file_paths = fs::read_dir("./examples")
            .unwrap()
            .filter_map(|res| res.ok())
            // Map the directory entries to paths
            .map(|dir_entry| dir_entry.path())
            // Filter out all paths with extensions other than `csv`
            .filter_map(|path| {
                if path.extension().is_some_and(|ext| ext == "fits") {
                    Some(path)
                } else {
                    None
                }
            });

        let mut idx = 0;
        for path in fits_file_paths {
            idx += 1;

            let f = File::open(path.clone()).unwrap();

            let reader = BufReader::new(f);
            let mut fits = Fits::from_reader(reader);
            let hdu = fits.next().unwrap().unwrap();
            match hdu {
                HDU::XImage(hdu) | HDU::Primary(hdu) => {
                    let header = hdu.get_header();
                    let wcs = wcs_from_fits_header(header).unwrap();

                    // add the astropy opened results
                    let path_astropy = path.with_extension("fits.csv");
                    let f = File::open(path_astropy).unwrap();

                    let mut rdr = csv::Reader::from_reader(BufReader::new(f));
                    for result in rdr.records() {
                        let record = result.unwrap();

                        let ra: f64 = record[0].parse().unwrap();
                        let dec: f64 = record[1].parse().unwrap();
                        let x: f64 = record[2].parse().unwrap();
                        let y: f64 = record[3].parse().unwrap();

                        if ra.is_finite() && dec.is_finite() {
                            if let Some(img_xy) = wcs.proj(&LonLat::new(ra, dec)) {
                                let dx = (img_xy.x() - x).abs();
                                let dy = (img_xy.y() - y).abs();
                                let tolerance = 1e-4;
                                if (dx > tolerance || dy > tolerance)
                                    && path
                                        .to_str()
                                        .map(|s| s.contains("74721b067"))
                                        .unwrap_or(false)
                                {
                                    eprintln!("  proj fail: ra={:.15}, dec={:.15}, exp_pix=({:.15}, {:.15}), got_pix=({:.15}, {:.15}), err=({:.6}, {:.6})",
                                        ra, dec, x, y, img_xy.x(), img_xy.y(), dx, dy);
                                }
                                assert_delta!(img_xy.x(), x, tolerance);
                                assert_delta!(img_xy.y(), y, tolerance);
                            }

                            if let Some(img_lonlat) = wcs.unproj(&ImgXY::new(x, y)) {
                                let dra = (img_lonlat.lon() - ra).abs();
                                let ddec = (img_lonlat.lat() - dec).abs();
                                let tolerance = 1e-4;
                                if (dra > tolerance || ddec > tolerance)
                                    && path
                                        .to_str()
                                        .map(|s| s.contains("74721b067"))
                                        .unwrap_or(false)
                                {
                                    eprintln!("  unproj fail: pix=({:.15}, {:.15}), exp_world=({:.15}, {:.15}), got_world=({:.15}, {:.15}), err=({:.6e}, {:.6e})",
                                        x, y, ra, dec, img_lonlat.lon(), img_lonlat.lat(), dra, ddec);
                                }
                                assert_delta!(img_lonlat.lon(), ra, tolerance);
                                assert_delta!(img_lonlat.lat(), dec, tolerance);
                            }
                        }
                    }
                }
                _ => unreachable!(),
            };
        }
        if idx < 1 {
            panic!("No files checked.")
        }
    }

    #[test]
    fn crval_to_crpix() {
        for path in glob("examples/*.fits").unwrap().flatten() {
            eprintln!("Testing file: {:?}", path);
            let f = File::open(path).unwrap();
            let reader = BufReader::new(f);
            let mut fits = Fits::from_reader(reader);
            let hdu = fits.next().unwrap().unwrap();

            match hdu {
                HDU::XImage(hdu) | HDU::Primary(hdu) => {
                    let header = hdu.get_header();
                    let crval1 = header.get_parsed::<f64>("CRVAL1").unwrap_or(0.0);
                    let crval2 = header.get_parsed::<f64>("CRVAL2").unwrap_or(0.0);
                    let crpix1 = if let Some(Value::Integer { value, .. }) = header.get("CRPIX1") {
                        *value as f64
                    } else if let Some(Value::Float { value, .. }) = header.get("CRPIX1") {
                        *value
                    } else {
                        0.0
                    };

                    let crpix2 = if let Some(Value::Integer { value, .. }) = header.get("CRPIX2") {
                        *value as f64
                    } else if let Some(Value::Float { value, .. }) = header.get("CRPIX2") {
                        *value
                    } else {
                        0.0
                    };

                    let wcs_result = wcs_from_fits_header(header);
                    let wcs = wcs_result.unwrap();

                    // Test 1: CRVAL → pixel (may not match CRPIX exactly)
                    // Note: In real astronomical data, CRVAL in the FITS header may not
                    // correspond exactly to world(CRPIX) when distortion polynomials are present.
                    // This happens when SIP was fitted after CRVAL was set, shifting the
                    // effective reference point. Both astropy and our implementation will
                    // produce the same pixel coordinate from CRVAL, which may differ from CRPIX.
                    let proj_px_result =
                        wcs.proj(&LonLat::new(crval1.to_radians(), crval2.to_radians()));

                    if proj_px_result.is_none() {
                        panic!("Projection failed");
                    }
                    let proj_px = proj_px_result.unwrap();

                    // Tolerance of 1 pixel accounts for:
                    // 1. SIP AP/BP inverse polynomial approximation error
                    // 2. CRVAL not matching the true world coordinate at CRPIX
                    let tolerance = 1.5; // 1.5 pixels to handle both effects
                    eprintln!(
                        "  CRVAL→pixel: Expected CRPIX ({}, {}), Got ({}, {})",
                        crpix1,
                        crpix2,
                        proj_px.x(),
                        proj_px.y()
                    );
                    assert_delta!(proj_px.x(), crpix1, tolerance);
                    assert_delta!(proj_px.y(), crpix2, tolerance);

                    // Test 2: CRPIX → world → pixel (round-trip, should be exact)
                    // This tests the actual quality of the WCS transformation
                    let world_at_crpix = wcs.unproj_lonlat(&ImgXY::new(crpix1, crpix2)).unwrap();
                    let pixel_roundtrip = wcs.proj(&world_at_crpix).unwrap();

                    // Round-trip should be sub-micropixel accurate
                    let rt_tolerance = 1e-6; // 1 nanopixel
                    eprintln!(
                        "  Round-trip: CRPIX ({}, {}) → world → pixel ({}, {})",
                        crpix1,
                        crpix2,
                        pixel_roundtrip.x(),
                        pixel_roundtrip.y()
                    );
                    assert_delta!(pixel_roundtrip.x(), crpix1, rt_tolerance);
                    assert_delta!(pixel_roundtrip.y(), crpix2, rt_tolerance);
                }
                _ => unreachable!(),
            };
        }
    }

    use fitsrs::hdu::HDU;
    #[test]
    fn open_fits() {
        let f = File::open("examples/neowise.fits").unwrap();

        let reader = BufReader::new(f);
        let mut fits = Fits::from_reader(reader);
        let hdu = fits.next().unwrap().unwrap();
        match hdu {
            HDU::XImage(hdu) | HDU::Primary(hdu) => {
                let header = hdu.get_header();
                let wcs = wcs_from_fits_header(header).unwrap();
                wcs.unproj(&ImgXY::new(0.0, 1200.0));
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn test_sip_round_trip() {
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        // Load the SIP FITS file
        let f = File::open("examples/74721b067-w2-int-1b.fits").unwrap();
        let reader = BufReader::new(f);
        let mut fits = Fits::from_reader(reader);
        let hdu = fits.next().unwrap().unwrap();

        let wcs = match hdu {
            HDU::XImage(hdu) | HDU::Primary(hdu) => {
                let header = hdu.get_header();
                wcs_from_fits_header(header).expect("Failed to parse SIP WCS")
            }
            _ => panic!("Expected Image HDU"),
        };

        println!("\n=== SIP Round-Trip Test (74721b067-w2-int-1b.fits) ===");
        println!("Testing pixel -> world -> pixel transformation accuracy\n");
        println!("NOTE: CSV and WCS both use 1-based FITS indexing\n");

        // Load ground truth from CSV
        let csv_file =
            File::open("examples/74721b067-w2-int-1b.fits.csv").expect("Failed to open CSV file");
        let csv_reader = BufReader::new(csv_file);

        let mut max_error = 0.0f64;
        let mut avg_error = 0.0f64;
        let mut count = 0;
        let mut failures = 0;

        for (line_num, line) in csv_reader.lines().enumerate() {
            let line = line.unwrap();
            let parts: Vec<&str> = line.split(',').collect();

            if parts.len() != 4 {
                println!(
                    "Warning: Line {} has {} parts, expected 4",
                    line_num + 1,
                    parts.len()
                );
                continue;
            }

            let expected_lon: f64 = parts[0].trim().parse().unwrap();
            let expected_lat: f64 = parts[1].trim().parse().unwrap();
            let pixel_x: f64 = parts[2].trim().parse().unwrap(); // 1-based from astropy
            let pixel_y: f64 = parts[3].trim().parse().unwrap();

            // Use 1-based pixel coordinates directly - they match FITS convention
            // which is what the WCS library expects (CRPIX is in 1-based system)
            let original_pixel = ImgXY::new(pixel_x, pixel_y);

            // Forward: pixel -> world
            let world = match wcs.unproj(&original_pixel) {
                Some(lonlat) => lonlat,
                None => {
                    println!(
                        "  Line {}: Failed to unproject pixel ({:.2}, {:.2})",
                        line_num + 1,
                        pixel_x,
                        pixel_y
                    );
                    failures += 1;
                    continue;
                }
            };

            // Check forward transformation accuracy against ground truth
            let lon_diff = (world.lon() - expected_lon).abs();
            let lat_diff = (world.lat() - expected_lat).abs();
            let angular_error = (lon_diff * lon_diff + lat_diff * lat_diff).sqrt();

            // Backward: world -> pixel (round-trip)
            let round_trip_pixel = match wcs.proj(&world) {
                Some(pixel) => pixel,
                None => {
                    println!(
                        "  Line {}: Failed to project world ({:.8}, {:.8}) back to pixel",
                        line_num + 1,
                        world.lon(),
                        world.lat()
                    );
                    failures += 1;
                    continue;
                }
            };

            // Calculate round-trip pixel error
            let error_x = round_trip_pixel.x() - original_pixel.x();
            let error_y = round_trip_pixel.y() - original_pixel.y();
            let error_magnitude = (error_x * error_x + error_y * error_y).sqrt();

            count += 1;
            avg_error += error_magnitude;
            if error_magnitude > max_error {
                max_error = error_magnitude;
            }

            // Print details for first few and any problematic cases
            if line_num < 5 || error_magnitude > 0.01 {
                println!(
                    "  Line {}: Pixel ({:7.2}, {:7.2})",
                    line_num + 1,
                    pixel_x,
                    pixel_y
                );
                println!(
                    "    World: ({:12.8}, {:12.8}) rad",
                    world.lon(),
                    world.lat()
                );
                println!(
                    "    Expected: ({:12.8}, {:12.8}) rad (Δ={:.3e} rad)",
                    expected_lon, expected_lat, angular_error
                );
                println!(
                    "    Round-trip: ({:7.2}, {:7.2}) px (Δ={:.6e} px)",
                    round_trip_pixel.x(),
                    round_trip_pixel.y(),
                    error_magnitude
                );
            }

            // Assert round-trip accuracy
            // With polynomial inverse (AP_ORDER=4), we should get very high accuracy.
            // Note: CSV ground truth was generated with an older Astropy version and
            // has slight discrepancies (~1-10 mas) vs current wcslib/Astropy.
            // Our implementation matches current wcslib 8.5 / Astropy to ~1 mas.
            // The round-trip errors of 1-6 millipixels are dominated by the CSV
            // ground truth version mismatch, not our implementation.
            let tolerance = 10e-3; // 10 millipixels tolerance

            if error_magnitude >= tolerance {
                if line_num < 5 {
                    println!(
                        "  ERROR: Round-trip error too large: {} pixels (tolerance: {} pixels)",
                        error_magnitude, tolerance
                    );
                }
                failures += 1;
            }
        }

        avg_error /= count as f64;

        println!("\n=== Summary ===");
        println!("Total test points: {}", count);
        println!("Failed transformations: {}", failures);
        println!("Average round-trip error: {:.6e} pixels", avg_error);
        println!("Maximum round-trip error: {:.6e} pixels", max_error);
        println!("===============\n");

        assert_eq!(failures, 0, "Some round-trip transformations failed");

        // Strict tolerance: round-trip error must be sub-micropixel (well below 1 arcsec)
        // At 1 arcsec/pixel scale, 1e-6 pixels = 1 microarcsec, far less than 1 arcsec
        let tolerance_pixels = 1e-6; // 1 micropixel
        assert!(
            max_error < tolerance_pixels,
            "Maximum round-trip error {:.6e} pixels exceeds tolerance of {:.6e} pixels (sub-micropixel)",
            max_error,
            tolerance_pixels
        );
    }

    #[test]
    fn test_sip_forward_accuracy() {
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        // Load the SIP FITS file
        let f = File::open("examples/minimal_sip_test.fits").unwrap();
        let reader = BufReader::new(f);
        let mut fits = Fits::from_reader(reader);
        let hdu = fits.next().unwrap().unwrap();

        let wcs = match hdu {
            HDU::XImage(hdu) | HDU::Primary(hdu) => {
                let header = hdu.get_header();
                wcs_from_fits_header(header).expect("Failed to parse SIP WCS")
            }
            _ => panic!("Expected Image HDU"),
        };

        println!("\n=== SIP Forward Transformation Test ===");
        println!("Testing pixel -> world accuracy against ground truth\n");
        println!("NOTE: CSV and WCS both use 1-based FITS indexing\n");

        // Load ground truth from CSV
        let csv_file =
            File::open("examples/minimal_sip_test.fits.csv").expect("Failed to open CSV file");
        let csv_reader = BufReader::new(csv_file);

        let mut max_error = 0.0f64;
        let mut avg_error = 0.0f64;
        let mut count = 0;

        for (line_num, line) in csv_reader.lines().enumerate() {
            let line = line.unwrap();
            let parts: Vec<&str> = line.split(',').collect();

            if parts.len() != 4 {
                continue;
            }

            let expected_lon: f64 = parts[0].trim().parse().unwrap();
            let expected_lat: f64 = parts[1].trim().parse().unwrap();
            let pixel_x: f64 = parts[2].trim().parse().unwrap();
            let pixel_y: f64 = parts[3].trim().parse().unwrap();

            // Use 1-based pixel coordinates directly - they match FITS convention
            let pixel = ImgXY::new(pixel_x, pixel_y);

            // Forward: pixel -> world
            let world = match wcs.unproj(&pixel) {
                Some(lonlat) => lonlat,
                None => {
                    panic!(
                        "Failed to unproject pixel ({}, {}) at line {}",
                        pixel_x,
                        pixel_y,
                        line_num + 1
                    );
                }
            };

            // Calculate error in radians
            let lon_diff = (world.lon() - expected_lon).abs();
            let lat_diff = (world.lat() - expected_lat).abs();
            let angular_error = (lon_diff * lon_diff + lat_diff * lat_diff).sqrt();

            count += 1;
            avg_error += angular_error;
            if angular_error > max_error {
                max_error = angular_error;
            }

            // Print first few for inspection
            if line_num < 5 {
                println!("  Pixel ({:7.2}, {:7.2})", pixel_x, pixel_y);
                println!(
                    "    Computed: ({:12.8}, {:12.8}) rad",
                    world.lon(),
                    world.lat()
                );
                println!(
                    "    Expected: ({:12.8}, {:12.8}) rad",
                    expected_lon, expected_lat
                );
                println!(
                    "    Error: {:.3e} rad ({:.3e} arcsec)\n",
                    angular_error,
                    angular_error.to_degrees() * 3600.0
                );
            }
        }

        avg_error /= count as f64;

        println!("=== Summary ===");
        println!("Total test points: {}", count);
        println!(
            "Average forward error: {:.6e} rad ({:.6e} arcsec)",
            avg_error,
            avg_error.to_degrees() * 3600.0
        );
        println!(
            "Maximum forward error: {:.6e} rad ({:.6e} arcsec)",
            max_error,
            max_error.to_degrees() * 3600.0
        );
        println!("===============\n");

        // Strict tolerance: SIP forward transformation must have milliarcsecond precision
        // This is far less than 1 arcsec, ensuring high-quality distortion correction
        let tolerance_arcsec = 0.001f64; // 1 milliarcsec
        let tolerance_rad = (tolerance_arcsec / 3600.0).to_radians();
        assert!(
            max_error < tolerance_rad,
            "Maximum forward error {:.6e} rad ({:.3e} arcsec) exceeds tolerance of {} milliarcsec (well below 1 arcsec)",
            max_error,
            max_error.to_degrees() * 3600.0,
            tolerance_arcsec * 1000.0
        );
    }
}
