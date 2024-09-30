# X-ray transits with batman

An adaption of [batman](https://github.com/lkreidberg/batman) for generating X-ray transit light curves, using a basic coronal emission model. You can find the paper detailing our model and discussing some of the examples presented here on NASA ADS (*link to be added - paper accepted to AAS Journals and will be published imminently*)

We thank Laura Kreidberg for providing the community with this excellent package, and her ongoing support of it. Her 2015 paper describing batman can be found [here](https://ui.adsabs.harvard.edu/abs/2015PASP..127.1161K/abstract)
 

## Installation

Our adapation of batman is based on version 2.5.1. The highest level versions of dependancies that we have tested our adapted code with are:

* Python 3.12.5
* numpy 1.26.4
* scipy 1.14.1
* matplotlib 3.9.2
* astropy 6.1.3

We recommend you install this code in a separate conda environment (or similar), especially if you also use batman for generating standard limb darkened transit profiles in the optical/near-IR. While the adapted code works as expected for the few optical transits we have generated with it, we have not comprehensively tested this. Therefore, we cannot guarantee there are no unexpected side effects of our changes to the code when using the standard in-built limb darkening laws.

Substantial changes in numpy version 2 current make it incompatable with the current batman code. Based on this and the other above requirements, we suggest the following for creating a new conda enviroment called "xrtransits":

`conda create --name xrtransits python numpy='<2' scipy matplotlib astropy`

Once this environment has been activated, `cd` into the batman-edited directory. The edited package can then be installed by running

`python -m pip install .`

You will need to `cd` out of this directory before using the package.

We also provide several useful functions for implementing our model in the xrt\_fns.py. Copy this script into your working directory. Use of these functions are described below and demonstrated in our various Jupyter notebooks included in this repo.


## Use

Our custom version of batman can in principle be used with any radially-symmetric model of emission, provided the relative brightnesses at 10001 points between x=0 (centre of the stellar disc) and x=1 (edge of the calculated emission region in 2D) are accurated calculated. These values must have been normalised such that the integrated intensity across the full emission region is 1, i.e.

$`\int^{2\pi}_0 \int^1_0 I(x)xdxd\theta = 1`$

To use our radially-symmetric, isothermal coronal model, we have provide code in xrt\_fns.get\_intVals for calculating the appropriate values at each of the 10001 required points, given an emission scale height h, and a number of emission scale heights you wish to calculate the coronal emission out to, Nscale \(the default is 6\). Note that h must have been converted to the batman coordinate system. To convert to and from units of stellar radii, we provide the functions xrt\_fns.hfromH and xrt\_fns.Hfromh. To find scale height from a given coronal temperature, use xrt\_fns.HfromT

```
import xt_fns as xf

#Stellar parameters
Rs = 0.780 #radius (Rsun)
Ms = 0.823 #mass (Msun)
T = 0.34 #kT of corona in keV

#get emission scale height in stellar radii
He = xf.HfromT(T, Rs, Ms)
#get emission scale height in batman coordinates
he = xf.hfromH(He)

#number of emission scale heights to calculate corona out to
Nscale = 6
#integral values
intVals = xf.get_intVals(he, Nscale)
```

The resulting array should be placed in the zeroth element of a list of length 6. This list should be set as the limb darkening parameter of a batman TransitParams object.

```
#set up the object to store the transit parameters
params = batman.TransitParams()

params.u = [0]*6
params.u[0] = intVals
```

batman handles the planet radius and semi-major axis in units of stellar radii. Since we have moved the location of the edge of the photospheric disc in the batman x=0 to x=1 coordinate system, we need to adjust these planetary values else they will be too large in the calculation. We multiply them by a factor Rx, the location of the edge of the stellar disc in batman coordinates.

```
#photosphere edge in batman coords
Rx = xf.calcRx(he, Nscale)

#planetary parameters
RpRs = 0.15641 #Rp/Rs of HD189733b
aRs = 8.863 #a/Rs of HD189733b

params.rp = RpRs * Rx        #planet radius (in units of stellar radii)
params.a = aRs * Rx          #semi-major axis (in units of stellar radii)
```

All other parameters can be set up in the normal way for batman, and there are also no changes to the model setup and normalised flux calculation steps.


## Examples

In this repo, we provide three example Jupyter Notebooks, closely linked to the examples discussed in detail in the paper.

**Examples.ipynb** - covers the examples in Section 3 of the paper. This includes setting up a basic model, and compares the effect of varying absorber radius, coronal scale height, and impact parameter on the transit morphology.

**Test_Detectability.ipynb** - covers the investigation into scaling relations for the effects of changing the various parameters on the overall detectability - Section 4 of the paper.

**Limitation_mitigation.ipynb** - covers examples of how to mitigate some of the limitations of the model discussed in Section 5 of the paper.




