# Coriolis-under-ice-positioning

## Coriolis under-ice positioning algorithm :


The algorithm implemented the “Terrain-following interpolation for under-ice floats” method presented by Kaihe Yamazaki during ADMT 22.

The method is described in Appendix A of Kaihe Yamazaki et al. article ([https://doi.org/10.1029/2019JC015406](https://doi.org/10.1029/2019JC015406)).

The “Terrain-following” method explained in this paper should be understood before reading the documentation and using the code.

The method has been improved by considering in situ data measured by the float to constrain the algorithm : 
- The average drift depth
- The max depth of the profile 
- The grounded flag 


In this branch, several improvements to the initial algorithm have been implemented. 
These improvements are described in more details here: 

- [performance improvements](https://github.com/cabanesc/Coriolis-under-ice-positioning/pull/1)

- [implementation of the geodesic computation](https://github.com/cabanesc/Coriolis-under-ice-positioning/pull/2)

- [algorithm improvements](https://github.com/cabanesc/Coriolis-under-ice-positioning/pull/3)



## Documentation
[How to use the code](https://github.com/euroargodev/Coriolis-under-ice-positioning/blob/main/estimate_profile_locations_V1.0_20220825.pdf)
