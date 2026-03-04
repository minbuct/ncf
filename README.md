# Hyperellipsoid Fitting with a Novel Cost Function (LMF-NCF)

The synthetic dataset consists of a Gaussian-noise dataset used to evaluate the robustness of the LMF-NCF algorithm. The naming convention for files is as follows:

1) The first character "G" represents "an instance of Gaussian-noise dataset".
2) The next two numbers (separated by "num") represent the noise level (standard deviation σ) and the instance index, respectively.
3) The file "PointCloud_Standard.txt" represents the ground-truth (noise-free) ellipsoid.

Examples:
- "G0.25num4.txt": Represents an instance containing Gaussian noise of σ=0.25, with instance index 4.
- "G0.1num10.txt": Represents an instance containing Gaussian noise of σ=0.1, with instance index 10.

----------------------------------------------------------------------------
The paper is under review. For prevention of copyright disputes, the MATLAB code and correpsonding test dataset will be normally available after this paper is published. Even so, the Matlab code is still attached and it is only for better peer-review. 

If you have any questions, please contact me.
