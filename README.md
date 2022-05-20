# actinobacteria_antibiotics
Storage of R scripts and bash scripts used to prepare and analyze high-throughput screening data
testing R studio

**Errors and Issues**

1. Warnings when running the function regarding the exact "p" values
  - Don't worry about that, since the test that is used is a Wilcoxon rank test, if there are ties between rank values, it is impossible to calculate the exact p value

2. Errors in the line where dataframes regarding spline or OD differences are are being "spread"
  - This shouldn't happen if you are using the newest version of the code. Don't be mad but the column selections were previously hardcoded, which meant that you needed to change it according to the input dataset. This was fixed.

3. Errors in blocks of code were certain dataframes are "spread"
  - This also shouldn't happen, but if it does it is because certain blocks of code have been run in the wrong order. Use the function!!! Don't run the code!!!

4. Missing a dataframe with the name "ttests"... at the end when trying to generate the ggplot graph
  - Older versions of the code might have this section in the wrong order (e.g., generating the dataframe after trying to plot it)

5. Propogating errors that just make the entire thing not work
  - Probably an issue with filtering out the first time point. For older datasets, there are only 8 time points, while for newer datasets there are anywhere between 9-10. If you filter out the first time point for things with only 8 time points, it will break the entire function
