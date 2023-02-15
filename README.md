# A conditional mutual information estimator

Conditional mutual information (CMI) is a model-free measure of correlation between two variables $X$ and $Y$ controlling for the effect of a third variable $Z$. In its simplest form, CMI can be expressed as follows

$$ I(X;Y|Z) = \sum_{z\in \mathcal{Z}} \sum_{y\in \mathcal{Y}} \sum_{x\in \mathcal{X}} p_{X,Y,Z}(x,y,z) \log \frac{p_Z(z)p_{X,Y,Z}(x,y,z)}{p_{X,Z}(x,z)p_{Y,Z}(y,z)} $$

In many real-life applications, you may need to control for more than one variable. For this purpose, I implemented a generalized CMI estimator that can calculate CMIs of the form

$$ I(X_1 \dots X_m;Y_1 \dots Y_n|Z_1 \dots Z_o) $$

for discrete random variables. If your dataset contains continuous fields, I suggest you discretize them, e.g. by binning according to quantile.