
# D2N2S
## stands for Dicom to Nifti to Struct

This repository extends SPM12 (https://github.com/spm/spm12) to greatly simplify scripting of diffusion MRI processing pipelines in Matlab.

The functions in this repository are intended to make preprocessing faster and prettier and less prone to user error.  

I implemented two main ideas. I thought all the metadata surrounding a DWI (Json files, bvecs, bvals) should be subsumed and associated in a Matlab structure very much like one created using spm_vol. This way, the data can be reference more easily, and less file-handling is needed.  

The second idea was to put Matlab arrays directly into SPM functions, instead of files.   

The best way to use these functions is to use dcm2niix on some dicom sequence folder, use d2n2s on the same folder, preprocess your data, and then use d2n2s_write() to write the preprocessed data.  

If you want to know more about usage, you should probably ask me, since I haven't written extensive instructions on how to use the functions anywhere.  

### outside code
The contents of .../utilities/external/ are all written by others. Their licenses can be found in comments in "choose_output.m", "fftw3.h", and "ringRm.cpp".
