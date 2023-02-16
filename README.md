![C++](https://img.shields.io/badge/c++-%2300599C.svg?style=flat-square&logo=c%2B%2B&logoColor=white) ![Python](https://img.shields.io/badge/python-3670A0?style=flat-square&logo=python&logoColor=ffdd54)

# A conditional mutual information estimator

Conditional mutual information (CMI) is a **model-free measure of correlation between two variables, $X$ and $Y$, controlling for the effect of a third variable $Z$**. In its simplest form, CMI can be expressed as follows

$$ I(X;Y|Z) = \sum_{z\in \mathcal{Z}} \sum_{y\in \mathcal{Y}} \sum_{x\in \mathcal{X}} p_{X,Y,Z}(x,y,z) \log \frac{p_Z(z)p_{X,Y,Z}(x,y,z)}{p_{X,Z}(x,z)p_{Y,Z}(y,z)} $$

In many real-life applications, one may need to control for more than one variable. For this purpose, I implemented a generalized CMI estimator that can calculate CMIs of the form

$$ I(X_1 \dots X_m;Y_1 \dots Y_n|Z_1 \dots Z_o) $$

for discrete random variables. If your dataset contains continuous fields, I suggest you discretize them, e.g. by binning according to quantile.

If you have **time-series data**, you may want to assess the impact of past states of X on current states of Y while also adjusting for past states of Y, i.e. $I(Y; X_{\mathrm{past}} |Y_{\mathrm{past}})$. This special case of CMI is also known as **transfer entropy**, denoted $T_{X\to Y}$.

## Build and usage

If you're on MacOS/Linux/Unix run ```./setup.sh``` file to build the shared library. If you're on Windows, you can get access to a free hosted Linux runtime through Google Colab.

### Pandas DataFrame support

Here's an example of how the CMI values can be estimated from Pandas DataFrames. This code would need to be contained in the project's root directory. **Ensure all variables are of float type**.

```python
import pandas as pd
import cmipy

data = pd.read_csv('your_path.csv')

cmi, p_val = cmipy.cond_mutual_info(data=data, x='price', y='area', z=['location', 'bathrooms'], p_samples=10000)

```

### Numpy array support

In addition to Pandas DataFrames, NumPy Arrays are also supported out-of-the-box. **Ensure all variables are of float type**.

```python
import numpy as np
import cmipy

x = np.random.randint(0, )

x = np.random.randint(0, 4, size=(10000, 1)).astype(float)
y = np.random.randint(0, 4, size=(10000, 1)).astype(float)
z = np.random.randint(0, 4, size=(10000, 1)).astype(float)

cmi, p_val = cmipy.cond_mutual_info(x=x, y=y, z=z, p_samples=10000, base=2.0)

```