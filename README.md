# Project-IV: Managing uncertainty in a bank's liquidity forecasts


It is essential for a bank to manage its inflows and outflows scrupulously to ensure that it has sufficient funds available to meet its commitments under all given scenarios. Many of these cash flows are extremely hard to control, but there are actions that a bank is able to take in order to impact the levels of certain flows. One key issue to consider is how much cash, or liquid assets a bank should hold in its reserves. It is not optimal to hold large amounts, however there are minimum amounts that are required to be held by the regulator. This is the main motivation for the project: how can a bank keep its liquidity levels as low as optimally possible? 


In a mathematical sense, we can build a model which we can analyse, by imposing stressed conditions on the forecasts of various cash flows, motivated by the current models that the bank has, as well as historical data that is available. It would be impractical to carry out individual exploration of certain scenarios, therefore we will aim to fit a statistical model, from which we can introduce the effects of various scenarios and perform an advanced uncertainty quantification based on various approaches in the analysis of uncertainty in complex models.



Atom bank is the UK’s first app-only bank, based in Durham, with main products including mortgages, savings accounts and business loans. As is the case with all banks, it is essential to have the right type and quantity of funds, in the right place, at the right time in order for these products to function properly. I am grateful to Atom bank for their support with my project, which has involved the understanding of the bank's liquidity forecasting model and the uncertainties around the inputs to the model. 

The outline of this project is as follows: firstly, I will introduce Bayes linear methods and the crucial Bayes linear update equations. I will present these alongside their application in state-of-the-art approaches in the analysis of uncertainty in complex models, such as emulation and history matching. Then, the concept of liquidity will be introduced, and certain liquidity ratios we are aiming to model, as well as Atom bank's necessary considerations regarding these. Finally, I will discuss various methods we can use to model these, and then perform an advanced uncertainty quantification using our earlier methods. There is a substantial use of coding in RStudio, for which the main functions and packages used are documented in this GitHub repository.
