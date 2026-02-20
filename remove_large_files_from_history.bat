@echo off
REM Remove large files from git history so push to GitHub succeeds.
REM Run from repo root. Backup or clone first if needed.

git filter-branch --force --index-filter "git rm --cached --ignore-unmatch testData/lightTest4DCT.p testData/veryLightDynMod.p testData/lightDynSeqWithMod.p opentps_core/opentps/core/examples/JupyterNotebooks/Data/SimpleRealDoseComputationOptimization_plan.tps opentps_core/opentps/core/examples/JupyterNotebooks/Data/Plan_WaterPhantom_cropped_resampled.tps opentps_core/opentps/core/examples/planOptimization/Plan_WaterPhantom_cropped_resampled.tps" --prune-empty --tag-name-filter cat -- --all

echo.
echo Done. If no errors, run: git push origin master --force
echo (Force push rewrites remote history; only do this on your own repo.)
pause
