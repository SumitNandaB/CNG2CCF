# Desktop entry point (used by Run_CNG2CCF.bat). Logic lives in cng2ccf_core.py.
# Outputs:
#   OutputCCF/                        - rotated CNG, in microns
#   OutputCCF_OriginalResolution/     - same trees rescaled to the Source resolution
#                                       (only written when the Source is not already in microns)
# Usage: python CNG2CCF_SWC_conversion_Updated.py [--scale 1000] [--no-original-resolution]
from cng2ccf_core import process_swc_files
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--scale", type=float, default=None, help="Force Source/CNG unit ratio; default auto-detect")
ap.add_argument("--no-original-resolution", action="store_true", help="Skip the OriginalResolution copies")
ap.add_argument("--no-refine", action="store_true", help="Use the 3-anchor rotation only (no whole-tree refinement)")
a = ap.parse_args()
process_swc_files("InputSource", "InputCNG", "OutputCCF", scale_override=a.scale,
                  original_res_folder=False if a.no_original_resolution else None,
                  refine=not a.no_refine)
