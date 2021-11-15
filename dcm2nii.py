from glob import glob

import dcmstack

src_dcms = glob(
    'TCGA-DD-A4NJ/03-30-2004-NA-MRI ABDOMEN WWO CONTRAST 74183-95782/17.000000-Ax Vibe Post5 MIN DELAYED-21418/*.dcm')
stacks = dcmstack.parse_and_stack(src_dcms)
stack = stacks.values
nii = stack.to_nifti()
nii.to_filename('output.nii.gz')
