# Statistical methods for detecting outlier evaluators

Statistical methods for hearing studies developed by my group. Please enter the branches. 

The R script is for detecting ‘outlier’ evaluators whose evaluation results tend to be higher or lower than their counterparts. In the first stage, evaluators’ effects are obtained by fitting a regression model. In the second stage, hypothesis tests are performed to detect ‘outlier’ evaluators, where we consider both the power of each hypothesis test and the false discovery rate (FDR) among all tests. A proposed FDR vs. Power decision plot is fitted in the second stage, and based on this plot the evaluator-specific significance levels can be selected by achieving an acceptable tradeoff between FDR and power.
