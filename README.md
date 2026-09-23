# CNG → CCF SWC Converter

Rotates NeuroMorpho.Org CNG SWC files back into the space of their Source SWC.
Unit differences between Source and CNG (nm vs µm, voxel vs µm) are detected
automatically.

Two outputs are produced:
- **OutputCCF/** – `<name>.CNG.CCF.swc`, always in micrometers.
- **OutputCCF_OriginalResolution/** – `<name>.CNG.CCF.OriginalResolution.swc`,
  with X, Y, Z *and* the radius multiplied back by the detected factor, so the
  tree overlays the Source file in its native units. Written only when scaling
  was actually needed; the folder is created on demand and left out otherwise.

## Rotation refinement
The 3-anchor SVD fit gives the orientation; a trimmed-ICP pass then fine-tunes the
*same* rotation about the root using every node, so small anchor-placement errors
do not leave a residual twist. The root stays exactly on the Source root and no
translation or scaling is introduced - it is a pure rotation, so the tree's shape
is untouched. The refined rotation is kept only if it improves the median
node-to-node distance by at least 5% and the correction is under 15 degrees;
otherwise the anchor-only result is used. Disable with `--no-refine`
(command line) or the checkbox under "Advanced" (web app).

## Files
- `cng2ccf_core.py` – all conversion logic (also a command-line tool)
- `app.py` – web front end (Streamlit)
- `CNG2CCF_SWC_conversion_Updated.py` + `Run_CNG2CCF.bat` – desktop batch use, as before
- `requirements.txt` – dependencies for the web app

## Desktop use
Put Source files in `InputSource/`, CNG files in `InputCNG/`, double-click
`Run_CNG2CCF.bat`. Options:

    python CNG2CCF_SWC_conversion_Updated.py --scale 1000           # force the factor
    python CNG2CCF_SWC_conversion_Updated.py --no-original-resolution
    python CNG2CCF_SWC_conversion_Updated.py --no-refine                # anchors only

## Run the web app locally
    pip install -r requirements.txt
    streamlit run app.py

The download is a zip with `CCF/` and, when scaling was needed,
`CCF_OriginalResolution/`, plus a conversion log.

## Put it online (Streamlit Community Cloud, free)
1. Create a GitHub repository and upload these files (app.py, cng2ccf_core.py, requirements.txt).
2. Sign in at https://share.streamlit.io with GitHub.
3. "Create app" → pick the repo, branch `main`, main file `app.py` → Deploy.
4. Share the resulting *.streamlit.app URL. Pushing to GitHub redeploys automatically.
