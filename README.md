# cyto2func
This repository is for all the amazing people trying to do a raw flow cytometry data conversion to cicuit component compatibility.


To do this we mainly focus on the fluorescence compared to a standard, in FSC_condition_data we correct for observed the variation in FSC. Since this has a covariance with the fluorescence channel this reduces the standard devaition to $\sqrt(1-\rho^2)\approx 0.72$ of the non-corrected data.

In Heat2compt we us this corrected data to find compatibility between gate compindations
