# Gaialaxy: making all-sky Milky Way images

---

### Summary

Gaialaxy is a Modern Fortran code designed to produce an all sky color
image encoding the mean integrated flux coming from all light sources
measured by the fantastic [Gaia
satellite](https://en.wikipedia.org/wiki/Gaia_(spacecraft)) (ESA),
i.e., from more than a billion stars from our Galaxy and its nearby
satellites.


By default, the fluxes are calibrated in W/m^2/sr and the output
colors are in the linear sRGB color space. Other units and color
transformation matrices can be trivially set in the main program file
(gaialaxy.f90).

### Compilation

Please ensure that you have a working installation of **gfortran** and
**gcc** compilers (or alternatives), the
[cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) and
[wcslib](https://www.atnf.csiro.au/people/mcalabre/WCS/wcslib/)
libraries. Editing the provided Makefile is certainly needed to
specify the install location of these libraries.

---

### Gaia Data

The (empty) directory "data/" should be filled with Gaia data under
the form of votable files (in the fits format) and named as
"gaiadata_01.fits", gaiadata_02.fits" etc...

The tables are expected to be of 6 columns (tfields) containing an
identificator of the source [source_id], its galactic coordinates [l]
and [b], the flux in the [G] band, the flux in the [RP] band and,
finally, the flux in the [BP] band. These data can be obtained
from the [ESA/Gaia/DPAC data archives](https://gea.esac.esa.int/archive/).

An example ADQL language query, which retrieves 100 sources in some sky
patch, reads:

```
curl -k -b cookies.txt -X POST  --form PHASE=run 
     --form LANG=ADQL --form REQUEST=doQuery --form FORMAT=fits \
     --form QUERY="SELECT TOP 100 source_id,l,b,phot_g_mean_flux,phot_rp_mean_flux,phot_bp_mean_flux \
     	    FROM gaiaedr3.gaia_source \
	    WHERE phot_g_mean_mag >=10 AND phot_rp_mean_mag >=10 AND phot_bp_mean_mag >=10 \
	    AND l<180.0 and b >=0" "https://gea.esac.esa.int/tap-server/tap/async"
```

Have fun!

---

### Example output

![gaialaxy.jpg](/docs/gaialaxy.jpg)

---

### Acknowledgements

This code uses data which are in the public domain and are provided by
the European Space Agency (ESA) mission
[Gaia](https://www.cosmos.esa.int/gaia), processed by the Gaia Data
Processing and Analysis Consortium
[DPAC](https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding
for the DPAC has been provided by national institutions, in particular
the institutions participating in the Gaia Multilateral Agreement.
