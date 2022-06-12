
include("../src/io/loadNii.jl")
include("../src/utils/generateMask.jl")
include("../src/utils/fitPoly3d.jl")
include("../src/utils/set_global_data.jl")
include("../src/utils/erodeMask.jl")
include("../src/unwrap/unwrapLaplacian.jl")
include("../src/inversion/rts.jl")

dllpath = "................../libs/libjulia_qsm.so"


B0 = 3
GYRO = 267.513 

filenameMag_jl = "............../magn_jl.nii"
filenamePhase = "............../phs_tissue_jl.nii"


nii_obj = loadNii(filenameMag_jl)
mag = nii_obj.img
nii_obj = loadNii(filenamePhase)
phas = nii_obj.img


set_nii_obj_info(nii_obj.hdr.dime.datatype)

vsz = nii_obj.hdr.dime.pixdim[2:4]
TE = nii_obj.hdr.dime.pixdim[5]


bdir = [0, 0, 1]
mask0 = generateMask(mag[:,:,:,end], vsz, "-m -n -f 0.5")
mask1 = erodeMask(mask0, 5);

uphas, _ = unwrapLaplacian(phas, mask1, vsz);

for t = 1:size(uphas, 4)
   uphas[:,:,:,t] = uphas[:,:,:,t] ./ (B0 * GYRO * TE[t]);
end

P = fitPoly3d(uphas, 4, mask1, vsz);
fl = uphas - mask1.*P;


x = rts(fl, mask1, vsz, bdir);

saveNii(outfile, x, vsz);

