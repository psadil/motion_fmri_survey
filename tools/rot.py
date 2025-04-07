import nibabel as nb
import numpy as np
from scipy import ndimage

nii = nb.load("../data/sub-2820279_ses-3_task-rest_bold01.nii.gz")

d = np.asarray(nii.get_fdata())
d2 = ndimage.rotate(d, angle=30, axes=(1, 2), reshape=False)
nb.nifti1.Nifti1Image(d2, affine=nii.affine, header=nii.header).to_filename(
    "/Users/psadil/Desktop/rmsd2/shifted.nii.gz"
)
