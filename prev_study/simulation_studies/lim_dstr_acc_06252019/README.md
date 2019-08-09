##During Creating

### Date : 06/25/2019
In this round of simulaitons, I am attempting to see what impact increasing the number of draws from the estiamted limiting distribution will have on the power of the test.  Due to the large number of people using the clusters right now, I am limiting the scope of my simulations.  Hopefully these simulaitons will let me know if it is worth it to run my final simulations with a large number of draws taken from the limiting distribution.

## Simulation Results

The tables show that there is not too much benefit to using a larger number of draws to estimate the limiting distribution.  While I didn't do any statistical tests to verify this, I feel as though there is no clear benefit to increasing the limiting distribution draws for norm selection.   The number of limiting distribution draws was varied for the l2, lp, max, ZL, and sum of square norms.  For the bonferoni method, while there are different columns, and there were separate simulations run for each column and row, the simulation settings were not altered for the different columns.  This was used to show the variation that would be expected just due to random variation.  400 separate tests were used to estimate power for each simulation setting. There are still some simulations that need to be finished. 
