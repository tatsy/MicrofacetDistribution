MicrofacetDistribution
===

> Simple rendering program using microfacet distribution for glossy surfaces.

## Build

```shell
$ git clone https://github.com/tatsy/MicrofacetDistribution.git
$ cd MicrofacetDistribution
$ mkdir build && cd build
$ cmake ..
$ make -j4
```

## Run

```shell
$ ./microfacet \
    --width [ image width ] \
    --height [ image height ] \
    --samples [ sample per pixel ] \
    --samplevis [ false: walter07, true: heitz14 ] \
    --distrib [ beckmann, ggx ] \
    --alphax [ x-axis roughness ] \
    --alphay [ y-axis roughness ] \
    --output [ output file name ]
```

## Results

**Beckmann** (ax = 0.1, ay = 0.1, 2000spp)

| Walter et al. 2007 | Heitz et al. 2014 |
|:------------------:|:-----------------:|
| ![beckmann_walter07_0.1_0.1][1] | ![beckmann_heitz14_0.1_0.1][2] |

  [1]: results/beckmann_walter07_0.1_0.1.png
  [2]: results/beckmann_heitz14_0.1_0.1.png

**GGX** (ax = 0.1, ay = 0.1, 2000spp)

| Walter et al. 2007 | Heitz et al. 2014 |
|:------------------:|:-----------------:|
| ![ggx_walter07_0.1_0.1][3] | ![ggx_heitz14_0.1_0.1][4] |

  [3]: results/ggx_walter07_0.1_0.1.png
  [4]: results/ggx_heitz14_0.1_0.1.png

**Beckmann** (ax = 0.5, ay = 0.1, 2000spp)

| Walter et al. 2007 | Heitz et al. 2014 |
|:------------------:|:-----------------:|
| ![beckmann_walter07_0.5_0.1][5] | ![beckmann_heitz14_0.5_0.1][6] |

  [5]: results/beckmann_walter07_0.5_0.1.png
  [6]: results/beckmann_heitz14_0.5_0.1.png

**GGX** (ax = 0.5, ay = 0.1, 2000spp)

| Walter et al. 2007 | Heitz et al. 2014 |
|:------------------:|:-----------------:|
| ![ggx_walter07_0.5_0.1][7] | ![ggx_heitz14_0.5_0.1][8] |

  [7]: results/ggx_walter07_0.5_0.1.png
  [8]: results/ggx_heitz14_0.5_0.1.png

**Beckmann** (ax = 0.5, ay = 0.1, 2000spp)

| Walter et al. 2007 | Heitz et al. 2014 |
|:------------------:|:-----------------:|
| ![beckmann_walter07_0.5_0.1][9] | ![beckmann_heitz14_0.5_0.1][10] |

  [9]: results/beckmann_walter07_0.5_0.5.png
  [10]: results/beckmann_heitz14_0.5_0.5.png

**GGX** (ax = 0.5, ay = 0.1, 2000spp)

| Walter et al. 2007 | Heitz et al. 2014 |
|:------------------:|:-----------------:|
| ![ggx_walter07_0.5_0.1][11] | ![ggx_heitz14_0.5_0.1][12] |

  [11]: results/ggx_walter07_0.5_0.5.png
  [12]: results/ggx_heitz14_0.5_0.5.png


## Note

In this program, I denote a ray direction towards lights as `wo` and that from an eye as `wi`, which follow the notation in [Walter et al. 2007] but are opposite from the current standard notations.

## References

* Walter et al., "Microfacet Models for Refraction through Rough Surfaces", Eurographics Symposium on Rendering, 2007.
* Heitz et al., "Importance Sampling Microfacet-Based BSDFs with the Distribution of Visible Normals", Eurographics Symposium on Rendering, 2014.