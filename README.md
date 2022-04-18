# How Much Demand Flexibility Could Have Spared Texas from the 2021 Outage?

[![GitHub commit](https://img.shields.io/github/last-commit/tamu-engineering-research/2021TXPowerOutage)](https://github.com/tamu-engineering-research/2021TXPowerOutage/commits/master) &nbsp;
[![GitHub license](https://img.shields.io/badge/license-MIT-yellow)](https://choosealicense.com/licenses/mit/)


## Suggested Citation 
- Please cite the following paper when you use this data hub:  
`
Wu, Dongqi, Xiangtian Zheng, Ali Menati, Lane Smith, Bainan Xia, Yixing Xu, Chanan Singh, and Le Xie. "How Much Demand Flexibility Could Have Spared
Texas from the 2021 Outage?", Working Paper
`

This repository contains the code and data associated with the above working paper.

## Features
- This repository contains four parts: 
	1) Processed public data during the 2021 Texas Power Outage
	2) Estimated load profile and demand sector share data 
	3) Synethetic network that resembles the ERCOT grid
	4) Source code used to run simulated event reproduction and hypothetical scenarios with flexible demand.


## Navigation
This data hub mainly contains five components: source data, released data, supplementary resources, parser codes, and  quick start tutorials. We navigate this data hub as follows.

- `From EIA (Generation, Demand, Demand Forecasts, and Interchange by BA)` contains public data for ERCOT, MISO and SPP as well as estimated counterfactual solar and wind generation produced by the model from Breakthrough Energy Sciences (BES). All data are in csv format and have been pre-processed to follow the same format. The first column is the name of the RTO region; The second column is timestamp when the data are sampled and the other are data values.
- `case` contains the network model for Matpower. The `case_Mar3_5pm.mat` is a stored Matpower case object that can be directly load by Matlab.
- `code` contains the main Matlab script used to run the various simulation scenarios in the paper.
- `data` contains additional data that are not complete or not official. It contains the raw ERCOT demand profile data, estimated sector percentage data and the population and geography related data to determine the load sector share of all load buses in the network.

## Contact Us
Please contact us if you need further technical support or search for cooperation. Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.\
Email contact: &nbsp; [Le Xie](mailto:le.xie@tamu.edu?subject=[GitHub]%20TX_Outage), &nbsp; [Dongqi Wu](mailto:dqwu@tamu.edu?subject=[GitHub]%20TX_Outage), &nbsp; [Xiangtian Zheng](mailto:zxt0515@tamu.edu?subject=[GitHub]%20TX_Outage).
