# RANSAC-Based-Image-Segmentation
Here we provide a modified version MATLAB imlpementation of robust regression based image segmentation.

In this project a robust regression approach is proposed for segmenting an image into background and foreground layers.
Background is defined as the smooth part of the image and foreground as the graphics and texts.
We then use RANSAC algorithm to robustly fit the smooth model to the image, and treat the outliers of that model as foreground pixels.

For further details please look at the following papers:

[1] Shervin Minaee and Yao Wang, "Screen Content Image Segmentation Using Robust Regression and Sparse Decomposition", IEEE Journal on Emerging and Selected Topics in Circuits and Systems, no.99, pp.1-12, 2016.
[2] Shervin Minaee and Yao Wang, "Image Decomposition Using a Robust Regression Approach." arXiv preprint arXiv:1609.03874, 2016.

