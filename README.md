# DMDpapers
codes from the two DMD papers. 
[1] A. Alassaf, L. Fan, "Dynamic Mode Decomposition in Various Power System Applications," NAPS 2019. pdf, ppt can be found from USF SPS's group website: power.eng.usf.edu--> publications-> conference articles. 
[2] A. Alassaf and L. Fan, "Randomized Dynamic Mode Decomposition for Oscillation Modal Analysis," accepted, IEEE trans. Power Systems. pdf

The first paper creates a couple of data and DMD is used in a similar way as FFT to find the oscillation frequency and the related phasors. 

The second paper analyzed five real-world ISO New England PMU data, replicated measurements after extracting eigenvalues and the coefficients (phi, b). Furthermore, for frequency measurements at different locations, phi (complex number) of eigenvalues are plots to judge whether a mode is interarea of internal. If a mode is an internal mode, we should see complex phasors at different locations seperated by angles greater than 90 degree. Otherwise, all phasors reside in a cone of an acute angle. 
