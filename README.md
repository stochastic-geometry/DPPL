# Determinantal Point Process Learning (DPPL)
### Matlab scripts for the paper "Machine Learning meets Stochastic Geometry: Determinantal Subset Selection for Wireless Networks"
#### Authors: Chiranjib Saha and Harpreet S. Dhillon
##### Email for correspondence: csaha@vt.edu

This repository contains the matlab scripts of DPPL proposed in the paper "Machine Learning meets Stochastic Geometry:
Determinantal Subset Selection for Wireless
Networks". 

Run "TrainDPP.m" to generate the results.
Under the folder "GenerateTrainingSet", use generateTrainingSet.m to generate new training sets.

Please cite the following paper if the code is reused. 
```
@article{saha2019load,
  title={Machine Learning meets Stochastic Geometry: {D}eterminantal Subset Selection for Wireless Networks},
  author={Saha, Chiranjib and Dhillon, Harpreet S},
  note={available online: arxiv.org/abs/1905.00504},
  year={2019}
}
```
**Acknowledgements.**
This repository reuses codes from the following sources. 
- ggplab: A Matlab toolbox for geometric programming (link: https://web.stanford.edu/~boyd/ggplab/)
  Note: please download the matlab scripts of geometric programming from this link and put the folder in "GenerateTrainingSets"
- hpaulkeeler/DetPoisson_MATLAB: A github repositoty for determinantal point process
- Codes for the paper --
 Gillenwater, J., Kulesza, A., & Taskar, B. (2012). Near-optimal map inference for determinantal point processes. In Advances in Neural  Information Processing Systems (pp. 2735-2743).
  Link of source code: http://jgillenw.com/dpp-map.html 

