# Southern_Ocean_Lead
Code for the manuscript "Antarctic icebergs reorganize ocean circulation during Pleistocene glacials" (Nature, 589(7841), 236-241) by Starr et al., 2021. (https://doi.org/10.1038/s41586-020-03094-7)

Please cite as: 
~~~~
@article{starr2021antarctic,
  title={Antarctic icebergs reorganize ocean circulation during Pleistocene glacials},
  author={Starr, Aidan and Hall, Ian R and Barker, Stephen and Rackow, Thomas and Zhang, Xu and Hemming, Sidney R and van der Lubbe, HJL and Knorr, Gregor and Berke, Melissa A and Bigg, Grant R and others},
  journal={Nature},
  volume={589},
  number={7841},
  pages={236--241},
  year={2021},
  publisher={Nature Publishing Group}
}
~~~~

## /raw
Containing the raw data, also archived at https://doi.pangaea.de/10.1594/PANGAEA.921315

## /analysis
Sub-directory containing the Matlab code for analysis of the raw data. Split into 4 parts (.m files) called by master.m. The dependencies "Ar1Sur.m" and "nexcf.m" are from http://tocsy.pik-potsdam.de/nest.php and are detailed in 1,2,3.

## /outputs
Containing the tables created by scripts in /analysis

## /pyberg_results
Containing some outputs from the Pyberg model runs performed for this manuscript. The full Pyberg code can be found at https://github.com/trackow/pyberg



----
## Useful Contacts regarding code or data
Contact:
- Aidan Starr (lead author) [here](mailto:StarrA1@Cardiff.ac.uk) or [here](mailto:Aidan.M.Starr@gmail.com)
- Ian R. Hall (senior co-author and IODP expedition co-PI) [here](mailto:Hall@Cardiff.ac.uk)
- Thomas Rackow (led the iceberg modeling and author of the pyberg iceberg module) [here](mailto:rackow@awi.de)
- Xu Zhang (for further info on the COSMOS model) [here](mailto:zhang@hotmailcom)
- Addresses for all co-authors are available through the online article

## Selected (code-related) References
References
1. Rehfeld, K., Marwan, N., Heitzig, J., Kurths, J.: Comparison of correlation analysis techniques for irregularly sampled time series, Nonlin. Proc. Geophys., 18(3), 389-404, 2011.
2. Rehfeld, K., Marwan, N., Breitenbach, S., Kurths, J.: Comparison of correlation analysis techniques for irregularly sampled time series, Climate Dynamics, Late Holocene Asian Monsoon Dynamics from small but complex paleoclimate networks, 41(1), 3-19 2013.
3. Rehfeld, K., Kurths, J.: Similarity measures for irregular and age uncertain time series, Clim. Past., 10, 107-122, 2014.
4. Barker, S., Chen, J., Gong, X., Jonkers, L., Knorr, G., & Thornalley, D.: Icebergs not the trigger for North Atlantic cold events, Nature, 520(7547), 333-336. 2015.
5. Rackow, T., Wesche, C., Timmermann, R., Hellmer, H. H., Juricke, S., and Jung, T.: A simulation of small to giant Antarctic iceberg evolution: Differential impact on climatology estimates, J. Geophys. Res. Oceans, 122, 3170–3190, 2017.
