# WAVEFRONT
Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG


## Versioning

We use [SemVer](http://semver.org/) for versioning. 

## Authors

Corresponding Authors: 

    Name: Alireza Chamanzar and Pulkit Grover  
    Institution: Carnegie Mellon University  
    Address: Hamerschlag Hall B200, 5000 Forbes Ave., Pittsburgh PA 15213 United States  
    Emails: achamanz@andrew.cmu.edu  ,  pgrover@andrew.cmu.edu 
   

## License

This project is licensed - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

This work was supported, in part, by grants from the National Science Foundation (NSF), the Chuck Noll Foundation for Brain Injury Research, the Center for Machine Learning and Health at CMU, under the Pittsburgh Health Data Alliance, the Neil and Jo Bushnell Fellowship in Engineering, the Hsu Chang Memorial Fellowship, the CMU Swartz Center for Entrepreneurship Innovation Commercialization Fellows program, and by the Office of the Assistant Secretary of Defense for Health Affairs, through the Defense Medical Research and Development Program under Award No. W81XWH-16-2-0020. Dr. Elmer’s research time was supported by the National Institutes of Health (NIH) through grant 5K23NS097629. Opinions, interpretations, conclusions, and recommendations are those of the authors and are not necessarily endorsed by the Department of Defense. We thank Maysamreza Chamanzar, Shilpa George, Neil Mehta, David Okonkwo, and Praveen Venkatesh for helpful discussions.

## How to cite this work?

You can cite our papers in [1,2], and our WAVEFRONT software as below:

Chamanzar, A., Elmer, J., Shutter, L., Hartings, J. A., and Grover, P., "WAVEFRONT: open source code and software", GitHub. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8210380.svg)](https://doi.org/10.5281/zenodo.8210380)


## Notes

The WAVEFRONT Pipeline contains 6 Preprocessing scripts and 1 Detection Algorithm script. 

This code was written using eeglab version 2020.0 and is dependent on functions from this package; using eeglab of matching version is highly encouraged. 

## File Structure 
Users are free to use their own file organization strategies. However, the code assumes that the files are stored in the following format:

- *Folder of Patient IDs*
   - Sub-ID A
     - Sub-ID A
	    - *All files for Sub-ID A*
   - Sub-ID B
   - Sub-ID C

## Brief Description of Functions for Each Script, in Sequential Order
### StepMII 
- Load raw files (EEG, ECoG, Impedance, ECG, PLETH, and RESP data; as well as Annotation and CSD events) and combines into streamlined set of files
- Bandpass filter EEG and ECoG data between [0.1 - 50] Hz 
- Process ECG, PLETH and RESP data

### StepMI
- Lowpass EEG data [0 - 30] Hz 
- Resample EEG data from 256 Hz to 64 Hz

### StepZero
- Number CSD Events
- Combine EEG and ECoG data into single file 

### StepI
- Extract Delta Band of EEG and ECoG signal [0.5 - 4] Hz
- Normalize and remove outliers from ECoG signal

### StepII
- Normalize and remove outliers from EEG signal

### StepIII
- Filter events 
- Mask temporally isolated portions of EEG data
- Extract Power Envelope of EEG data
- Cross-correlate Power Envelope with First-Derivative Kernel

### Detection Algorithm 
- Load Processed Data from Preprocessing Pipeline 
- Project data into 2D space using electrode locations and spatial conversions
- Smooth projections and create binary images representing relevant data patterns
- Calculate Optical Flows from binary images 
- Score Flows based on speed and orientation
- At temporal frames with sufficient Optical Flow scores, detect Spreading Depolarization
- Compare detection results with ground truths for each Data Session and Threshold Combination 


## References

[1] Chamanzar, A., Elmer, J., Shutter, L. et al. Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG. Commun Med 3, 113 (2023). https://doi.org/10.1038/s43856-023-00344-3

[2] Chamanzar, A., George, S., Venkatesh, P., Chamanzar, M., Shutter, L., Elmer, J., and Grover, P. (2018). An algorithm for automated, noninvasive detection of cortical spreading depolarizations based on EEG simulations. IEEE Transactions on Biomedical Engineering, 66(4), 1115-1126. https://doi.org/10.1109/TBME.2018.2867112
