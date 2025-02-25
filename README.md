# Statistical methods for detecting outlier evaluators

The R script is for detecting ‘outlier’ evaluators whose evaluation results tend to be higher or lower than their counterparts. In the first stage, evaluators’ effects are obtained by fitting a regression model. In the second stage, hypothesis tests are performed to detect ‘outlier’ evaluators, where we consider both the power of each hypothesis test and the false discovery rate (FDR) among all tests. A proposed FDR vs. Power decision plot is fitted in the second stage, and based on this plot the evaluator-specific significance levels can be selected by achieving an acceptable tradeoff between FDR and power. See the paper below for details.

Wu Y, Rosner B, Curhan G, Curhan S, Wang M. Analytical methods for detecting outlier evaluators. BMC Medical Research Methodology. 2023;23(1):177. doi: 10.1186/s12874-023-01988-4

