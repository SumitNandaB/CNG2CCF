"""CNG -> CCF SWC converter: web front end (Streamlit)."""
import io
import os
import tempfile
import zipfile

import streamlit as st

from cng2ccf_core import convert_pair, output_stem, pair_files

st.set_page_config(page_title="CNG → CCF SWC Converter", page_icon="🧠")
st.title("CNG → CCF SWC Converter")
st.write(
    "Rotates NeuroMorpho.Org **CNG** SWC files back into the original (e.g. CCF) space of their "
    "**Source** SWC. Unit differences (nm vs µm, voxel vs µm) are detected automatically. "
    "The download contains **CCF/** in micrometers, plus **CCF_OriginalResolution/** "
    "(coordinates and radii scaled back to the Source resolution) whenever scaling was needed. "
    "Files are paired by name (`X.swc` ↔ `X_CNG.swc` / `X.CNG.swc`)."
)

c1, c2 = st.columns(2)
src_up = c1.file_uploader("Source SWC files", type=["swc"], accept_multiple_files=True)
cng_up = c2.file_uploader("CNG SWC files", type=["swc"], accept_multiple_files=True)

with st.expander("Advanced"):
    refine = st.checkbox("Refine the rotation on all nodes (recommended)", value=True,
                         help="After the 3-anchor fit, a trimmed ICP pass fine-tunes the same "
                              "rotation about the root. Kept only if it measurably improves the fit.")
    manual = st.checkbox("Override automatic scale detection")
    scale = st.number_input("Source units per micron (e.g. 1000 for nm, 1/voxel size for voxels)",
                            min_value=1e-6, value=1000.0, format="%.6g", disabled=not manual)

if st.button("Convert", type="primary", disabled=not (src_up and cng_up)):
    logs, rows = [], []
    with tempfile.TemporaryDirectory() as tmp:
        sdir, cdir, odir, rdir = (os.path.join(tmp, d) for d in ("src", "cng", "out", "out_orig"))
        for d in (sdir, cdir, odir, rdir):
            os.makedirs(d)
        for folder, files in ((sdir, src_up), (cdir, cng_up)):
            for f in files:
                with open(os.path.join(folder, os.path.basename(f.name)), "wb") as fh:
                    fh.write(f.getbuffer())

        pairs, notes = pair_files([os.path.join(sdir, f.name) for f in src_up],
                                  [os.path.join(cdir, f.name) for f in cng_up])
        for n in notes:
            st.warning(n)

        prog = st.progress(0.0)
        for i, (s, c) in enumerate(pairs, 1):
            out = os.path.join(odir, f"{output_stem(c)}.CNG.CCF.swc")
            orig = os.path.join(rdir, f"{output_stem(c)}.CNG.CCF.OriginalResolution.swc")
            logs.append(f"{os.path.basename(s)} + {os.path.basename(c)}")
            try:
                r = convert_pair(s, c, out, scale if manual else None,
                                 log=lambda m: logs.append(m), original_res_path=orig, refine=refine)
                st_ = r["refine"]
                rows.append({"neuron": output_stem(c), "scale": f"1/{r['factor']:g}", "units": r["kind"],
                             "max anchor residual (µm)": round(max(r["anchor_residuals_um"]), 3),
                             "median node distance (µm)": round(st_["median_after_um"], 3) if st_ else "—",
                             "refined by (deg)": round(st_["angle_deg"], 2) if st_ and st_["accepted"] else "—",
                             "warnings": "; ".join(r["warnings"])})
            except Exception as e:
                logs.append(f"  [Error] {e}")
                rows.append({"neuron": output_stem(c), "scale": "—", "units": "FAILED", "warnings": str(e)})
            prog.progress(i / len(pairs))

        outputs, originals = sorted(os.listdir(odir)), sorted(os.listdir(rdir))
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for fn in outputs:
                z.write(os.path.join(odir, fn), f"CCF/{fn}")
            for fn in originals:
                z.write(os.path.join(rdir, fn), f"CCF_OriginalResolution/{fn}")
            z.writestr("conversion_log.txt", "\n".join(logs))

    if rows:
        st.dataframe(rows, use_container_width=True)
    if outputs:
        msg = f"Converted {len(outputs)} of {len(pairs)} pair(s)."
        msg += (f" {len(originals)} also written at the original Source resolution."
                if originals else " No scaling was needed, so no OriginalResolution copies.")
        st.success(msg)
        st.download_button("Download results (.zip)", buf.getvalue(),
                           file_name="CCF_output.zip", mime="application/zip")
    with st.expander("Log"):
        st.code("\n".join(logs))

st.caption("Questions: NeuroMorpho.Org · Sumit Nanda (snanda@mednet.ucla.edu)")
